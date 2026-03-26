# Helper: accept either an output_dir path or a pre-loaded list
.load_iadhore <- function(x) {
  if (is.character(x)) {
    read_iadhore_output(x)
  } else if (is.list(x) && all(c("genes", "multiplicons", "segments") %in% names(x))) {
    x
  } else {
    stop("`output_dir` must be a path string or a list from read_iadhore_output()")
  }
}

# Count the total number of base positions covered by a set of [start, end]
# intervals (both 0-based, inclusive). Returns an integer.
.count_covered <- function(starts, ends, size) {
  if (length(starts) == 0L) return(0L)
  starts <- pmax(0L, starts)
  ends   <- pmin(size - 1L, ends)
  valid  <- starts <= ends
  starts <- starts[valid]; ends <- ends[valid]
  if (length(starts) == 0L) return(0L)
  covered <- logical(size)
  for (i in seq_along(starts)) {
    covered[(starts[i] + 1L):(ends[i] + 1L)] <- TRUE  # convert to 1-based
  }
  sum(covered)
}


#' Compute colinear portions across genomes
#'
#' For every gene list in every genome, reports how many gene positions (and
#' what percentage) are covered by at least one non-redundant multiplicon
#' segment that is shared with each other genome.  Reimplements the original
#' i-ADHoRe utility \code{colinear_portions.pl} using parsed output tables.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   already returned by \code{\link{read_iadhore_output}}.
#'
#' @return A data frame with one row per (genome, list, other_genome) triple
#'   and columns:
#'   \describe{
#'     \item{genome}{Reference genome}
#'     \item{list}{Gene list (chromosome / scaffold)}
#'     \item{list_size}{Number of genes in the list}
#'     \item{other_genome}{Genome to which colinearity is measured}
#'     \item{colinear_count}{Gene positions colinear with \code{other_genome}}
#'     \item{colinear_pct}{Percentage of the list colinear with \code{other_genome}}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' cp <- colinear_portions("output/arath/")
#' # genome-level summary
#' aggregate(cbind(colinear_count, list_size) ~ genome + other_genome,
#'           data = cp, FUN = sum)
#' }
colinear_portions <- function(output_dir) {
  dat         <- .load_iadhore(output_dir)
  genes        <- dat$genes
  multiplicons <- dat$multiplicons
  segments     <- dat$segments

  # List sizes: max(remapped_coordinate) + 1 per (genome, list).
  # remapped_coordinate is the position after tandem collapsing, matching the
  # Perl API's $genelist->unmapped_size and $gene->coordinate.
  ls_df <- aggregate(remapped_coordinate ~ genome + list, data = genes,
                     FUN = function(x) max(x) + 1L)
  names(ls_df)[3] <- "list_size"

  # Gene-ID → remapped coordinate lookup
  gene_coords <- genes[, c("id", "remapped_coordinate")]

  # Keep only non-redundant multiplicons
  nr_ids  <- multiplicons$id[multiplicons$is_redundant == 0]
  seg_nr  <- segments[segments$multiplicon %in% nr_ids, , drop = FALSE]

  # Attach first / last gene remapped coordinates
  seg_nr <- merge(seg_nr, gene_coords, by.x = "first", by.y = "id",
                  all.x = TRUE, sort = FALSE)
  names(seg_nr)[names(seg_nr) == "remapped_coordinate"] <- "first_coord"
  seg_nr <- merge(seg_nr, gene_coords, by.x = "last",  by.y = "id",
                  all.x = TRUE, sort = FALSE)
  names(seg_nr)[names(seg_nr) == "remapped_coordinate"] <- "last_coord"

  # One row per (multiplicon, genome) pair — distinct genomes in each multiplicon
  mult_genomes <- unique(seg_nr[, c("multiplicon", "genome")])

  # Cross-join to get all (ref_genome, other_genome) pairs per multiplicon
  cross <- merge(seg_nr,
                 mult_genomes,
                 by        = "multiplicon",
                 suffixes  = c("", "_other"),
                 sort      = FALSE)
  cross <- cross[cross$genome != cross$genome_other, , drop = FALSE]

  # Group key: ref genome + list + other genome
  cross$key <- paste(cross$genome, cross$list, cross$genome_other, sep = "\x1F")

  ls_key <- paste(ls_df$genome, ls_df$list, sep = "\x1F")

  keys    <- unique(cross$key)
  results <- vector("list", length(keys))

  for (i in seq_along(keys)) {
    k     <- keys[i]
    rows  <- cross[cross$key == k, ]
    parts <- strsplit(k, "\x1F", fixed = TRUE)[[1]]

    sz_idx <- which(ls_key == paste(parts[1], parts[2], sep = "\x1F"))
    size   <- ls_df$list_size[sz_idx]

    n <- .count_covered(rows$first_coord, rows$last_coord, size)

    results[[i]] <- data.frame(
      genome        = parts[1],
      list          = parts[2],
      list_size     = size,
      other_genome  = parts[3],
      colinear_count = n,
      colinear_pct  = n / size * 100,
      stringsAsFactors = FALSE
    )
  }

  if (length(results) == 0L) {
    return(data.frame(genome = character(), list = character(),
                      list_size = integer(), other_genome = character(),
                      colinear_count = integer(), colinear_pct = numeric(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out[order(out$genome, out$list, out$other_genome), ]
}


#' Compute duplicated / multiplicated portions of each genome
#'
#' For every gene list in every genome, reports how many positions (and what
#' percentage) are covered at each multiplication level (1 = singleton, 2 =
#' found in 2 copies, 3 = found in 3 copies, …).  Each position is assigned
#' the *highest* level observed across all non-redundant multiplicons.
#' Reimplements the original i-ADHoRe utility \code{multiplicated_portions.pl}.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   already returned by \code{\link{read_iadhore_output}}.
#'
#' @return A data frame with one row per (genome, list, level) and columns:
#'   \describe{
#'     \item{genome}{Genome name}
#'     \item{list}{Gene list (chromosome / scaffold)}
#'     \item{list_size}{Number of genes in the list}
#'     \item{level}{Multiplication level (integer, \eqn{\geq 1})}
#'     \item{count}{Number of positions at this level}
#'     \item{pct}{Percentage of the list at this level}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' mp <- multiplicated_portions("output/arath/")
#' # genome-level % duplicated (level >= 2)
#' dup <- mp[mp$level >= 2, ]
#' tot <- aggregate(count ~ genome, data = mp, FUN = sum)
#' dup_sum <- aggregate(count ~ genome, data = dup, FUN = sum)
#' merge(dup_sum, tot, by = "genome", suffixes = c("_dup", "_total"))
#' }
multiplicated_portions <- function(output_dir) {
  dat         <- .load_iadhore(output_dir)
  genes        <- dat$genes
  multiplicons <- dat$multiplicons
  segments     <- dat$segments

  # List sizes (remapped = after tandem collapsing, matching Perl API)
  ls_df <- aggregate(remapped_coordinate ~ genome + list, data = genes,
                     FUN = function(x) max(x) + 1L)
  names(ls_df)[3] <- "list_size"

  # Gene-ID → remapped coordinate lookup
  gene_coords <- genes[, c("id", "remapped_coordinate")]

  # Keep only non-redundant multiplicons
  nr_ids <- multiplicons$id[multiplicons$is_redundant == 0]
  seg_nr <- segments[segments$multiplicon %in% nr_ids, , drop = FALSE]

  # Attach first / last gene remapped coordinates
  seg_nr <- merge(seg_nr, gene_coords, by.x = "first", by.y = "id",
                  all.x = TRUE, sort = FALSE)
  names(seg_nr)[names(seg_nr) == "remapped_coordinate"] <- "first_coord"
  seg_nr <- merge(seg_nr, gene_coords, by.x = "last",  by.y = "id",
                  all.x = TRUE, sort = FALSE)
  names(seg_nr)[names(seg_nr) == "remapped_coordinate"] <- "last_coord"

  # Multiplication level for a genome within a multiplicon = number of
  # segments of that genome in the multiplicon
  lev_df <- aggregate(id ~ multiplicon + genome, data = seg_nr, FUN = length)
  names(lev_df)[3] <- "level_in_genome"
  seg_nr <- merge(seg_nr, lev_df, by = c("multiplicon", "genome"), sort = FALSE)

  # Initialise per-position level arrays (all 1)
  ls_key    <- paste(ls_df$genome, ls_df$list, sep = "\x1F")
  pos_level <- lapply(ls_df$list_size, function(sz) rep(1L, sz))
  names(pos_level) <- ls_key

  # Update positions with their maximum observed level
  seg_key <- paste(seg_nr$genome, seg_nr$list, sep = "\x1F")
  for (i in seq_len(nrow(seg_nr))) {
    k <- seg_key[i]
    if (!k %in% names(pos_level)) next
    sz  <- length(pos_level[[k]])
    s   <- max(1L, seg_nr$first_coord[i] + 1L)  # remapped 0-based → 1-based R index
    e   <- min(sz,  seg_nr$last_coord[i]  + 1L)
    lev <- seg_nr$level_in_genome[i]
    if (s <= e) {
      pos_level[[k]][s:e] <- pmax(pos_level[[k]][s:e], lev)
    }
  }

  # Tabulate levels per list
  results <- vector("list", nrow(ls_df))
  for (i in seq_len(nrow(ls_df))) {
    k    <- ls_key[i]
    sz   <- ls_df$list_size[i]
    tab  <- table(pos_level[[k]])
    results[[i]] <- data.frame(
      genome    = ls_df$genome[i],
      list      = ls_df$list[i],
      list_size = sz,
      level     = as.integer(names(tab)),
      count     = as.integer(tab),
      pct       = as.integer(tab) / sz * 100,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out[order(out$genome, out$list, out$level), ]
}


#' Print a text summary of i-ADHoRe results
#'
#' Prints a human-readable overview of collinearity and duplication statistics
#' to the console, combining \code{\link{colinear_portions}} and
#' \code{\link{multiplicated_portions}}.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   already returned by \code{\link{read_iadhore_output}}.
#' @param digits Number of decimal places for percentages (default 2).
#'
#' @return Invisibly returns a list with elements \code{colinear} and
#'   \code{multiplicated} (the underlying data frames).
#' @export
#' @examples
#' \dontrun{
#' iadhore_summary("output/arath/")
#' }
iadhore_summary <- function(output_dir, digits = 2) {
  dat <- .load_iadhore(output_dir)

  mults <- dat$multiplicons
  n_total      <- nrow(mults)
  n_nonred     <- sum(mults$is_redundant == 0)
  max_level    <- max(mults$level, na.rm = TRUE)

  cat("=== i-ADHoRe run summary ===\n\n")
  cat(sprintf("  Total multiplicons     : %d\n", n_total))
  cat(sprintf("  Non-redundant          : %d\n", n_nonred))
  cat(sprintf("  Maximum level          : %d\n", max_level))
  cat("\n")

  # ---- Collinearity ----
  cp <- colinear_portions(dat)
  if (nrow(cp) > 0) {
    cat("--- Colinear portions (genome level) ---\n\n")
    genome_totals <- aggregate(
      cbind(colinear_count, list_size) ~ genome + other_genome,
      data = cp, FUN = sum)
    genome_totals$pct <- genome_totals$colinear_count /
      genome_totals$list_size * 100

    fmt_cp <- paste0("  %s vs %s : %.", digits, "f%%\n")
    for (g in sort(unique(genome_totals$genome))) {
      sub <- genome_totals[genome_totals$genome == g, ]
      for (j in seq_len(nrow(sub))) {
        cat(sprintf(fmt_cp, sub$genome[j], sub$other_genome[j],
                    round(sub$pct[j], digits)))
      }
    }
    cat("\n")
  }

  # ---- Duplication ----
  mp <- multiplicated_portions(dat)
  cat("--- Duplicated portions (genome level) ---\n\n")
  dup      <- mp[mp$level >= 2, ]
  tot_cnt  <- aggregate(count ~ genome, data = mp,  FUN = sum)
  dup_cnt  <- aggregate(count ~ genome, data = dup, FUN = sum)
  names(dup_cnt)[2] <- "dup_count"
  summ     <- merge(tot_cnt, dup_cnt, by = "genome", all.x = TRUE)
  summ$dup_count[is.na(summ$dup_count)] <- 0L
  summ$pct <- summ$dup_count / summ$count * 100

  fmt_mp <- paste0("  %s : %.", digits, "f%% duplicated\n")
  for (j in seq_len(nrow(summ))) {
    cat(sprintf(fmt_mp, summ$genome[j], round(summ$pct[j], digits)))
  }
  cat("\n")

  invisible(list(colinear = cp, multiplicated = mp))
}
