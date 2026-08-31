## ---------------------------------------------------------------------------
##  read_chrom.R -- format-agnostic chromatogram folder reader.
##
##  Implements read_chrom_folder(), which
##    * scans a folder for every extension chromConverter can parse, plus
##      directory-based formats (.D, .raw) that list.files() would miss;
##    * dispatches vendor binaries to chromConverter::read_chroms();
##    * falls back to a dialect-sniffing text reader for .csv/.tsv/.txt/.arw,
##      so plain two-column exports work with no configuration and without
##      chromConverter installed at all;
##    * normalises everything to the (time_min, intensity) + metadata schema
##      the rest of the pipeline expects;
##
##  Requires data.table. chromConverter is optional and only loaded when a
##  vendor format is actually encountered.
## ---------------------------------------------------------------------------


# ---------- Format registry ---------------------------------------------------

## `ext`        file extension, lower case, no dot
## `format_in`  value to pass to chromConverter::read_chroms(format_in = )
## `container`  "file" or "dir" -- .D and Waters .raw are directories, so
##              list.files() never sees them as candidates
## `engine`     "chromconverter" or "text"
## `note`       parsers that need a separate manual install
CHROM_FORMATS <- data.table::data.table(
  ext = c(
    "arw",  "csv",   "tsv",   "txt",  "dat",
    
    "ch",   "uv",    "mwd",   "fid",  "ms",   "dx",  "d",  "sp",  "dad",
    "lcd",  "gcd",   "qgd",   "c0",
    "raw",  "sms",   "chrom", "mdf",
    
    "mzml", "mzxml", "cdf",   "wiff", "asm"
  ),
  format_in = c(
    "waters_arw", "csv", "csv", NA, NA,
    
    "chemstation_ch", "chemstation_uv", "chemstation_uv", "chemstation_fid",
    "chemstation_ms", "agilent_dx",     "agilent_d",      "masshunter_dad",  "masshunter_dad",
    "shimadzu_lcd",   "shimadzu_gcd",   "shimadzu_qgd",   "shimadzu_fid",
    "waters_raw",     "varian_sms",     "chromatotec",    "mdf",
    
    "mzml", "mzxml", "cdf", "other", "asm"
  ),
  container = c(
    rep("file", 5L),
    
    "file", "file", "file", "file", "file", "file", "dir", "file", "file",
    "file", "file", "file", "file",
    "dir",  "file", "file", "file",
    
    "file", "file", "file", "file", "file"
  ),
  engine = c(
    "text", "text", "text", "text", "text",
    rep("chromconverter", 9L),
    rep("chromconverter", 4L),
    rep("chromconverter", 4L),
    rep("chromconverter", 5L)
  ),
  note = c(
    rep(NA_character_, 5L),
    rep(NA_character_, 9L),
    rep(NA_character_, 4L),
    "needs ThermoRawFileParser or rainbow", NA, NA, NA,
    NA, NA, NA, "needs OpenChrom/entab", NA
  )
)

## .txt and .dat are ambiguous: Chromeleon and Shimadzu both export .txt.
## They are handled by the text reader, which sniffs the dialect, but an operator
## who knows the exact format can force chromConverter with format_override.
CHROM_TEXT_EXT <- c("arw", "csv", "tsv", "txt", "dat")

chrom_supported_formats <- function() CHROM_FORMATS[]

.cc_available <- function() requireNamespace("chromConverter", quietly = TRUE)

# ---------- Candidate discovery -----------------------------------------------

#' List every parseable chromatogram entry in a folder.
#'
#' @param path folder to scan.
#' @param recursive descend into subfolders (off by default; note that a
#'   directory-based format like .D is itself a folder, and is always detected
#'   whether or not recursion is on).
#' @param extensions restrict to these extensions (lower case, no dot).
#' @param exclude_dirs exceptions which do not contain chromatograms
#'   do not place any other file into path other than the chromatograms
#' @return data.table(path, name, ext, format_in, engine, container)
find_chrom_files <- function(
    path,
    recursive    = FALSE,
    extensions   = NULL,
    exclude_dirs = c("peakfinder_QC", "peakfinder_inspect")
  ) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)

  reg <- CHROM_FORMATS
  if (!is.null(extensions)) reg <- reg[ext %in% tolower(extensions)]

  ## --- directory-based formats (.D, Waters .raw) -----------------------------
  dir_ext <- reg[container == "dir", ext]
  dirs <- list.dirs(path, recursive = recursive, full.names = TRUE)
  dirs <- setdiff(dirs, path)
  dir_hits <- data.table::data.table()
  if (length(dirs) && length(dir_ext)) {
    de <- tolower(tools::file_ext(dirs))
    ok <- de %in% dir_ext
    if (any(ok)) {
      dir_hits <- data.table::data.table(
        path = dirs[ok], name = basename(dirs[ok]), ext = de[ok],
        container = "dir"
      )
    }
  }

  ## --- ordinary files -------------------------------------------------------
  files <- list.files(path, recursive = recursive, full.names = TRUE,
                      all.files = FALSE, no.. = TRUE)
  files <- files[!dir.exists(files)]
  files <- files[basename(files) != "peakfinder_lookup.csv"] ## exclude lookup file if exists
  
  ## never pick up our own outputs, or files sitting inside a .D / .raw folder
  if (length(exclude_dirs)) {
    drop <- Reduce(`|`, lapply(exclude_dirs, function(dd)
      grepl(paste0("(^|/)", dd, "(/|$)"), files)), FALSE)
    files <- files[!drop]
  }
  if (nrow(dir_hits)) {
    inside <- Reduce(`|`, lapply(dir_hits$path, function(dd)
      startsWith(files, paste0(dd, "/"))), FALSE)
    files <- files[!inside]
  }

  file_hits <- data.table::data.table()
  if (length(files)) {
    fe <- tolower(tools::file_ext(files))
    ok <- fe %in% reg[container == "file", ext]
    if (any(ok)) {
      file_hits <- data.table::data.table(
        path = files[ok], name = basename(files[ok]), ext = fe[ok],
        container = "file"
      )
    }
  }

  out <- data.table::rbindlist(list(dir_hits, file_hits), use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(out)

  out <- merge(out, reg[, .(ext, format_in, engine, note)], by = "ext",
               all.x = TRUE, sort = FALSE)
  data.table::setorder(out, name)
  out[]
}

# ---------- Text reader with dialect sniffing ---------------------------------

## Vendor ASCII exports often disagree on `sep`, `dec`, and header structure.
## This helper finds the longest trailing block of lines that
## parse as a consistent numeric table and takes the dialect that produces it.
.sniff_text <- function(path, n_probe = 400L) {
  lines <- readLines(path, n = n_probe, warn = FALSE, encoding = "latin1")
  ## Keep ORIGINAL line numbering -- skip= is passed to fread in file
  ## coordinates, so blank lines must not be squeezed out here. Only trailing
  ## blanks are dropped, since some exporters end with one.
  while (length(lines) && !nzchar(trimws(lines[length(lines)])))
    lines <- lines[-length(lines)]
  if (length(lines) < 3L) return(NULL)

  seps <- c("\t", ";", ",", "|")
  best <- NULL

  for (sep in seps) for (dec in c(".", ",")) {
    if (identical(sep, dec)) next
    num_re <- if (dec == ".") {
      "^[+-]?(\\d+\\.?\\d*|\\.\\d+)([eEdD][+-]?\\d+)?$"
    } else {
      "^[+-]?(\\d+,?\\d*|,\\d+)([eEdD][+-]?\\d+)?$"
    }

    parts <- strsplit(lines, sep, fixed = TRUE)
    nf <- lengths(parts)
    numeric_row <- vapply(parts, function(p) {
      p <- trimws(p)
      p <- p[nzchar(p)]
      length(p) >= 2L && all(grepl(num_re, p))
    }, logical(1))

    ## longest suffix of numeric rows sharing one field count
    if (!any(numeric_row)) next
    k <- length(lines)
    while (k >= 1L && numeric_row[k]) k <- k - 1L
    start <- k + 1L
    if (start > length(lines)) next
    run     <- start:length(lines)
    tab_nf  <- nf[run]
    mode_nf <- as.integer(names(sort(table(tab_nf), decreasing = TRUE))[1L])
    run     <- run[tab_nf == mode_nf]
    if (length(run) < 3L) next

    cand <- list(sep = sep, dec = dec, start_line = min(run),
                 ncol = mode_nf, nrow = length(run))
    if (is.null(best) || cand$nrow > best$nrow ||
        (cand$nrow == best$nrow && cand$ncol > best$ncol)) best <- cand
  }

  if (is.null(best)) return(NULL)

  ## a header row is the line directly above the data block, if it is not itself
  ## numeric and splits into the same number of fields
  hdr <- NULL
  j   <- best$start_line - 1L
  while (j >= 1L && !nzchar(trimws(lines[j]))) j <- j - 1L   # step over blanks
  if (j >= 1L) {
    h <- trimws(strsplit(lines[j], best$sep, fixed = TRUE)[[1]])
    if (length(h) == best$ncol && any(nzchar(h))) hdr <- h
  }
  best$header <- hdr
  best$preamble <- if (best$start_line > 1L) lines[seq_len(best$start_line - 1L)] else character()
  best
}

#' Read a text chromatogram, sniffing the dialect.
#' @param wavelength if the file has more than two columns, pick the column
#'   whose header parses to the numeric value nearest this; otherwise column 2.
read_text_chrom <- function(path, wavelength = NULL) {
  sn <- .sniff_text(path)
  if (is.null(sn)) stop("Could not find a numeric data block in: ", basename(path))

  d <- data.table::fread(
    path,
    skip         = sn$start_line - 1L,
    sep          = sn$sep,
    dec          = sn$dec,
    header       = FALSE,
    fill         = TRUE,
    showProgress = FALSE,
    na.strings   = c("", "NA", "n.a.", "-"),
    select       = seq_len(sn$ncol),
    encoding     = "Latin-1"
  )
  if (!ncol(d)) stop("Empty data block in: ", basename(path))

  if (!is.null(sn$header) && length(sn$header) == ncol(d)) {
    nm <- make.unique(ifelse(nzchar(sn$header), sn$header, paste0("V", seq_len(ncol(d)))))
    data.table::setnames(d, nm)
  }

  ## every column must be numeric; fread may have read some as character if the
  ## decimal mark guess was wrong for a subset of rows
  for (j in seq_len(ncol(d))) {
    if (!is.numeric(d[[j]])) {
      v <- as.character(d[[j]])
      if (sn$dec == ",") v <- gsub(",", ".", v, fixed = TRUE)
      data.table::set(d, j = j, value = suppressWarnings(as.numeric(v)))
    }
  }

  tcol <- 1L
  if (ncol(d) == 1L) stop("Only one column found in: ", basename(path))

  if (ncol(d) == 2L) {
    icol <- 2L
    wl <- NA_real_
  } else {
    ## multi-wavelength export: choose by header value
    hv <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", names(d))))
    cand <- setdiff(seq_len(ncol(d)), tcol)
    if (!is.null(wavelength) && any(is.finite(hv[cand]))) {
      icol <- cand[which.min(abs(hv[cand] - wavelength))]
    } else {
      icol <- cand[1L]
    }
    wl <- hv[icol]
  }

  out <- data.table::data.table(
    time_min  = as.numeric(d[[tcol]]),
    intensity = as.numeric(d[[icol]])
  )
  data.table::setattr(out, "sniff", sn)
  data.table::setattr(out, "wavelength_from_header", wl)
  data.table::setattr(out, "n_channels", ncol(d) - 1L)
  out
}

# ---------- chromConverter dispatch -------------------------------------------

#' Read one vendor file via chromConverter and reduce it to two columns.
read_vendor_chrom <- function(path, format_in, wavelength = NULL,
                              read_metadata = TRUE, ...) {
  if (!.cc_available()) {
    stop("Reading '", format_in, "' needs the chromConverter package:\n",
         "  install.packages('chromConverter')   # or pak::pak('ethanbass/chromConverter')")
  }

  res <- chromConverter::read_chroms(
    paths         = path,
    format_in     = format_in,
    find_files    = FALSE,
    format_out    = "data.table",
    data_format   = "wide",
    read_metadata = read_metadata,
    progress_bar  = FALSE,
    ...
  )

  x <- if (is.list(res) && !data.table::is.data.table(res)) res[[1L]] else res
  meta <- attributes(x)
  m <- as.matrix(x)

  ## wide format is retention time x wavelength; the time axis lives in the
  ## row names for matrix output and may be the first column for data.table
  rn <- suppressWarnings(as.numeric(rownames(x)))
  if (all(is.finite(rn)) && length(rn) == nrow(m)) {
    tt <- rn
    vals <- m
  } else {
    tt <- suppressWarnings(as.numeric(m[, 1L]))
    vals <- m[, -1L, drop = FALSE]
  }

  hv <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", colnames(vals))))
  j <- if (!is.null(wavelength) && any(is.finite(hv))) {
    which.min(abs(hv - wavelength))
  } else 1L

  out <- data.table::data.table(
    time_min  = as.numeric(tt),
    intensity = as.numeric(vals[, j])
  )
  data.table::setattr(out, "vendor_metadata", meta[setdiff(names(meta), c("dim", "dimnames", "class", "names", "row.names"))])
  data.table::setattr(out, "wavelength_from_header", hv[j])
  data.table::setattr(out, "n_channels", ncol(vals))
  out
}

#' Read one chromatogram of any supported type.
read_chrom_file <- function(
    path,
    format_in  = NULL,
    engine     = NULL,
    wavelength = NULL,
    time_unit  = c("min", "s", "auto"),
    ...
  ) {
  time_unit <- match.arg(time_unit)
  this_ext  <- tolower(tools::file_ext(path))
  reg <- CHROM_FORMATS[ext == this_ext]
  if (is.null(engine))    engine    <- if (nrow(reg)) reg$engine[1L] else "text"
  if (is.null(format_in)) format_in <- if (nrow(reg)) reg$format_in[1L] else NA_character_

  d <- if (identical(engine, "text") || this_ext %in% CHROM_TEXT_EXT) {
    read_text_chrom(path, wavelength = wavelength)
  } else {
    read_vendor_chrom(path, format_in = format_in, wavelength = wavelength, ...)
  }

  d <- d[is.finite(time_min) & is.finite(intensity)]
  data.table::setorder(d, time_min)
  if (nrow(d) < 5L) stop("Fewer than five usable data points in: ", basename(path))

  ## Seconds vs minutes: a chromatogram that "runs" for hundreds of units with
  ## a sub-unit sampling step is almost certainly in seconds.
  span <- diff(range(d$time_min))
  if (time_unit == "s" || (time_unit == "auto" && span > 180)) {
    d[, time_min := time_min / 60]
    data.table::setattr(d, "time_rescaled", TRUE)
  }
  d
}

# ---------- Metadata plug-ins -------------------------------------------------

## A metadata parser takes (path, first_lines) and returns a named list of
## scalar columns. Return NULL for "nothing to add". This is the only piece
## that encodes a local naming convention, so swap it, don't edit the reader.

#' The original metadata convention with ARW files:
#'   should contain:    compound_column_eluent_modifier_temp_10xflow
#'   e.g.               timolol_iA3_MeOH_DEA_25_07 (25 °C and 0.7 mL/min)
arw_metadata_parser <- function(path, first_lines = NULL, ext_re = "arw|arq") {
  filename <- basename(path)
  ff <- regmatches(
    filename,
    regexec(paste0("^(?:.*?)(\\d{3})?nm?(\\d{4,})\\.(", ext_re, ")$"),
            filename, ignore.case = TRUE, perl = TRUE)
  )[[1]]
  if (!length(ff)) return(NULL)

  if (is.null(first_lines)) first_lines <- readLines(path, n = 2L, warn = FALSE)
  if (length(first_lines) < 2L) return(NULL)

  metadata <- trimws(first_lines[2L])
  metadata <- sub('^"', "", sub('"$', "", metadata))

  mm <- regmatches(
    metadata,
    regexec(paste0("^(.+)_([^_]+)_([^_]+)_([^_]+)_",
                   "(-?[0-9]+(?:\\.[0-9]+)?)_([0-9]+(?:\\.[0-9]+)?)$"),
            metadata, perl = TRUE)
  )[[1]]
  if (!length(mm)) return(NULL)

  list(
    wavelength_nm  = as.integer(ff[2L]),
    experimentalID = ff[3L],
    compoundname   = mm[2L],
    columnID       = mm[3L],
    eluentID       = mm[4L],
    modifierID     = mm[5L],
    temp_C         = as.numeric(mm[6L]),
    flow_mL_min    = as.numeric(mm[7L]) / 10,
    metadata_raw   = metadata
  )
}

#' Last-resort parser: take the stem as the experiment id and a trailing
#' 3-digit "###nm" as the wavelength if present. Never fails.
generic_metadata_parser <- function(path, first_lines = NULL) {
  stem <- tools::file_path_sans_ext(basename(path))
  wl <- suppressWarnings(as.integer(
    sub(".*?(\\d{3})\\s*nm.*", "\\1", stem, ignore.case = TRUE)))
  list(
    wavelength_nm  = if (is.finite(wl) && grepl("\\d{3}\\s*nm", stem, ignore.case = TRUE)) wl else NA_integer_,
    experimentalID = stem,
    compoundname   = NA_character_,
    columnID       = NA_character_,
    eluentID       = NA_character_,
    modifierID     = NA_character_,
    temp_C         = NA_real_,
    flow_mL_min    = NA_real_,
    metadata_raw   = NA_character_
  )
}

# ---------- Folder reader -----------------------------------------------------

#' Read every chromatogram in a folder, whatever the format.
#'
#' @param metadata_parsers list of parser functions, tried in order; the first
#'   that returns a non-NULL list wins. Default tries the ARW convention, then
#'   the generic one.
#' @param require_date_folder keep the "folder must start with YYYYMMDD"
#'   rule. Set FALSE to fall back to the folder mtime.
#' @param on_error "stop" (old behaviour) or "skip", which collects failures in
#'   attr(result, "failures") and carries on -- friendlier for a large folder
#'   batch where one file is truncated.
read_chrom_folder <- function(
    path,
    extensions          = NULL,
    recursive           = FALSE,
    wavelength          = NULL,
    time_unit           = c("min", "s", "auto"),
    metadata_parsers    = list(arw_metadata_parser, generic_metadata_parser),
    require_date_folder = TRUE,
    on_error            = c("stop", "skip"),
    hash_files          = TRUE,
    ...
  ) {

  time_unit <- match.arg(time_unit)
  on_error  <- match.arg(on_error)

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  folder_name <- basename(path)

  fm <- regmatches(folder_name,
                   regexec("^(\\d{8})(?:_(.*))?$", folder_name, perl = TRUE))[[1]]
  if (length(fm)) {
    date <- as.Date(fm[2L], format = "%Y%m%d")
    folder_label <- if (length(fm) >= 3L && !is.na(fm[3L])) fm[3L] else ""
  } else if (require_date_folder) {
    stop("Folder name must begin with YYYYMMDD, optionally followed by _junk: ",
         folder_name)
  } else {
    date <- as.Date(file.info(path)$mtime)
    folder_label <- folder_name
  }
  if (is.na(date)) stop("Invalid date for folder: ", folder_name)

  cand <- find_chrom_files(path, recursive = recursive, extensions = extensions)
  if (!nrow(cand)) {
    stop("No parseable chromatogram files found in: ", path, "\n",
         "Recognised extensions: ",
         paste(sort(CHROM_FORMATS$ext), collapse = ", "))
  }

  need_cc <- cand[engine == "chromconverter"]
  if (nrow(need_cc) && !.cc_available()) {
    warning(nrow(need_cc), " file(s) need chromConverter and will be skipped: ",
            paste(unique(need_cc$ext), collapse = ", "))
    cand <- cand[engine != "chromconverter"]
    if (!nrow(cand)) stop("Nothing left to read after skipping vendor formats.")
  }

  failures <- list()

  read_one <- function(k) {
    p <- cand$path[k]
    filename <- cand$name[k]

    d <- read_chrom_file(p, format_in = cand$format_in[k], engine = cand$engine[k],
                         wavelength = wavelength, time_unit = time_unit, ...)

    first_lines <- if (cand$engine[k] == "text") {
      tryCatch(readLines(p, n = 3L, warn = FALSE), error = function(e) NULL)
    } else NULL

    md <- NULL
    for (fn in metadata_parsers) {
      md <- tryCatch(fn(p, first_lines), error = function(e) NULL)
      if (!is.null(md)) break
    }
    if (is.null(md)) md <- generic_metadata_parser(p, first_lines)

    ## a wavelength recovered from the file itself beats one guessed from the name
    wl_hdr <- attr(d, "wavelength_from_header")
    if ((is.null(md$wavelength_nm) || !is.finite(md$wavelength_nm)) &&
        is.finite(wl_hdr)) md$wavelength_nm <- as.integer(round(wl_hdr))
    if (is.null(md$wavelength_nm) || !is.finite(md$wavelength_nm)) {
      md$wavelength_nm <- NA_integer_
    }

    file_hash <- if (isTRUE(hash_files) && requireNamespace("digest", quietly = TRUE) &&
                     cand$container[k] == "file") {
      digest::digest(p, algo = "xxhash64", file = TRUE)
    } else {
      fi <- file.info(p)
      substr(digest_fallback(paste(filename, fi$size, fi$mtime)), 1L, 16L)
    }

    logical_key <- paste(format(date, "%Y%m%d"),
                         md$experimentalID, md$wavelength_nm, sep = "__")
    record_uid <- paste0(logical_key, "__", substr(file_hash, 1L, 12L))

    d[, `:=`(
      record_uid     = record_uid,
      logical_key    = logical_key,
      date           = date,
      folder_label   = folder_label,
      source_folder  = path,
      filename       = filename,
      file_hash      = file_hash,
      extension      = cand$ext[k],
      source_format  = cand$format_in[k],
      n_channels     = attr(d, "n_channels"),
      wavelength_nm  = md$wavelength_nm,
      experimentalID = md$experimentalID,
      compoundname   = md$compoundname,
      columnID       = md$columnID,
      eluentID       = md$eluentID,
      modifierID     = md$modifierID,
      temp_C         = md$temp_C,
      flow_mL_min    = md$flow_mL_min,
      metadata_raw   = md$metadata_raw
    )]

    first_cols <- c("record_uid",    "logical_key",  "date",           "folder_label",
                    "filename",      "file_hash",    "wavelength_nm",  "experimentalID",
                    "compoundname",  "columnID",     "eluentID",       "modifierID",
                    "temp_C",        "flow_mL_min",  "time_min",       "intensity")
    data.table::setcolorder(d, c(first_cols, setdiff(names(d), first_cols)))
    d[]
  }

  pieces <- vector("list", nrow(cand))
  for (k in seq_len(nrow(cand))) {
    z <- tryCatch(read_one(k), error = function(e) e)
    if (inherits(z, "error")) {
      if (on_error == "stop") stop("Failed on ", cand$name[k], ": ", conditionMessage(z))
      failures[[cand$name[k]]] <- conditionMessage(z)
    } else {
      pieces[[k]] <- z
    }
  }

  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) stop("Every candidate file failed to read in: ", path)

  out <- data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
  data.table::setattr(out, "failures", failures)
  data.table::setattr(out, "candidates", cand)
  if (length(failures)) {
    message("read_chrom_folder: skipped ", length(failures), " file(s); see ",
            "attr(x, 'failures').")
  }
  out[]
}

## tiny non-cryptographic fallback so hashing works without digest installed
digest_fallback <- function(s) {
  b <- utf8ToInt(paste(s, collapse = "|"))
  h <- 2166136261
  for (v in b) h <- bitwAnd((bitwXor(h, v) * 16777619), 0xFFFFFFFF)
  sprintf("%08x%08x", h, bitwAnd(h * 2654435761, 0xFFFFFFFF))
}
