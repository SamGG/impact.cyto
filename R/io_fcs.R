#' Convert Navios file to FCS
#'
#' A Navios LMD file is FCS file that contains two segments of data. In the
#' first segment the data are coded in the FCS version 2, data being in log
#' scale. In the second segment, the data are stored in raw mode (20 bits), ie
#' at full resolution. Some information (description, scales) are not reported
#' in the second segment and some conversions are needed. This function extracts
#' the 2nd segment and export it as FCS 3.
#'
#' @param files	        optional character vector with filenames.
#' @param raw_dir       directory where to look for the files.
#' @param pattern       this argument is passed on to \code{dir}, see details.
#' @param fcs_dir       a directory name where FCS files are stored.
#' @param keep_range    logical, TRUE by default. If TRUE, the ranges of the
#'   resulting FCS are scaled to equal to the linear ranges of the first segment
#'   using the $PxE keywords. If FALSE, no range is scaled.
#' @param verbose       integer, verbosity, 0 being minimal.
#'
#' @return Nothing, but reports every steps of processing.
#'
#' @examples
#' \dontrun{
#' # Setup the file locations
#' raw_dir = "raw_files"
#' fcs_dir = "fcs_files"
#' # Then convert
#' convert.navios.FCS(raw_dir, fcs_dir, verbose = 2)
#' }
#'
#' @import flowCore
#'
#' @export
convert.navios.FCS <- function(files=NULL, raw_dir=".", pattern=NULL,
                               fcs_dir, keep_range = TRUE, verbose = 2)
{
  if (is.null(files)) {
    files <- dir(raw_dir, pattern, full.names = TRUE)
    if (length(files) < 1)
      stop(sprintf("No matching files found in \"%s\".", raw_dir))
  } else {
    if (!is.character(files))
      stop("'files' must be a character vector.")
    if (raw_dir != ".")
      files <- file.path(raw_dir, files)
  }

  if (!dir.exists(fcs_dir)) dir.create(fcs_dir)

  for (f in files) {
    if (verbose) message(sprintf("Processing file \"%s\":", basename(f)))

    # 1st segment -------------------------------------------
    ff = read.FCS(filename = f,
                  transformation = FALSE,
                  min.limit = FALSE,
                  truncate_max_range = FALSE,
                  dataset = 1)
    # check version is 2
    if (keyword(ff, "FCSversion") != "2") {
      message(sprintf("1st segment is not FCS v2!", f))
      next()
    }
    # check presence of next segment
    if (is.null(keyword(ff, "$NEXTDATA"))) {
      message("No NEXTDATA keyword. File skipped!")
      next()
    }
    # get parameters
    pd1 = pData(parameters(ff))
    if (verbose > 1) {
      cat("Parameters of 1st segment\n")
      print(pd1)
    }
    # compute the linear range using $PxE
    pd1$linRange = 0
    for (p in rownames(pd1)) {
      kwd = keyword(ff, paste0(p, "E"))
      if (is.null(kwd)) next()
      coef = as.numeric(strsplit(kwd[[1]], ",")[[1]])
      pd1[p, "linRange"] = 10^coef[1] * coef[2]
    }

    # 2nd segment -------------------------------------------
    ff = read.FCS(filename = f,
                  transformation = FALSE,
                  min.limit = FALSE,
                  truncate_max_range = FALSE,
                  dataset = 2)
    # check version is 3
    if (keyword(ff, "FCSversion") != "3") {
      message(sprintf("2nd segment is not FCS v3!", f))
      next()
    }
    # correct SPILLOVER: add name to columns
    ff.spillover = keyword(ff, "$SPILLOVER")[[1]]
    if (is.null(ff.spillover)) {
      message(sprintf("No $SPILLOVER!", f))
    } else {
      col.ids = as.numeric(colnames(ff.spillover))
      colnames(ff.spillover) = colnames(ff)[col.ids]
      keyword(ff) = list("$SPILLOVER" = ff.spillover)
    }
    # get parameters
    pd2 = pData(parameters(ff))
    if (verbose > 1) {
      cat("Parameters of 2nd segment\n")
      print(pd2)
    }
    # match channels by name between two segments
    pd2$id = gsub("(FL\\d+).+", "\\1", pd2$name, perl = TRUE)
    pd1$id = gsub("(FL\\d+).+", "\\1", pd1$name, perl = TRUE)
    pd.mrg = merge(pd2[, c("id", "name", "desc")],
                   pd1[, c("id", "name", "desc", "linRange")],
                   by = "id", sort = FALSE, stringsAsFactors = FALSE)
    # add description from 1st segment
    for (i in seq(nrow(pd.mrg))) {
      pd2[pd2$id == pd.mrg$id[i], "desc"] = pd.mrg[i, "desc.y"]
    }
    # change range according 1st segment (checked vs FlowJo 10.4.2)
    if (keep_range) {
      if (verbose > 1) cat("Scaling ranges\n")
      gain = rep(1, ncol(ff))
      rang = NULL
      for (i in seq(nrow(pd.mrg))) {
        # range from 1st segment
        range1st = pd.mrg$linRange[i]
        if (range1st == 0) next()
        # corresponding channel and range
        j = which(pd2$id == pd.mrg$id[i])
        range2nd = pd2$range[j]
        # gain factor to match true range
        gain[j] = (range1st-1) / (range2nd-1)
        # update range
        pd2$range[j] = range1st
        pd2$maxRange[j] = range1st - 1
        names(range1st) = paste0(rownames(pd2)[j], "R")
        rang = c(rang, range1st)
        names(rang)[]
        keyword(ff) = as.list(range1st)
      }
      exprs(ff) = sweep(exprs(ff), 2, gain, "*")
      # change datatype
      keyword(ff) = list("$DATATYPE" = "F")
    }
      # normalize F/S channel names
    ids = grep("^[FS]S-", pd2$name)
    if (length(ids)) {
      pd2$name = gsub("(^[FS]S)-", "\\1C-", pd2$name)
      pd2$desc[ids] = pd2$name[grep("^[FS]SC-", pd2$name)]
    }
    # remove extra columns
    pd2 = pd2[, colnames(pd2) != "id"]
    # update data columns according parameters
    colnames(ff) = pd2$name
    pData(parameters(ff)) = pd2
    if (verbose) {
      cat("Annotated parameters of 2nd segment\n")
      print(pd2)
    }
    # add comment
    com = keyword(ff, "$COM")
    keyword(ff) = list("$COM" = paste(com, "convert.navios.FCS v0.1",
                                      collapse = ifelse(is.null(com), "", "; ")))
    # write to disk
    write.FCS(ff, filename = file.path(fcs_dir, basename(f)))
  }
}
