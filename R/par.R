#' Find Putative Ancestral Regions (PARs)
#'
#' Identifies groups of syntenic segments that likely descend from the same
#' ancestral chromosomal region by applying a hierarchical clustering approach
#' to the i-ADHoRe output.  Follows the method of Jaillon \emph{et al.} (2009)
#' \emph{PNAS} 106:14712.
#'
#' \strong{Algorithm}
#' \enumerate{
#'   \item For every pair of i-ADHoRe segments (one from each genome) count the
#'     number of real anchorpoints they share.
#'   \item Compute a significance score \eqn{S = -\log_{10}(p)} where \eqn{p}
#'     is the Poisson probability of observing \eqn{\ge k} anchorpoints under a
#'     random null (expected proportional to segment sizes).
#'   \item Cluster segments within each genome by the Pearson correlation of
#'     their score profiles (average linkage, minimum correlation \code{r_cutoff}).
#'     Segments with correlated profiles match the same partner regions and
#'     are therefore likely duplicated copies of the same ancestral locus.
#'   \item A cluster pair (one cluster from each genome) whose maximum pairwise
#'     score exceeds \code{-log10(p_cutoff)} is called a PAR.
#' }
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   already returned by \code{\link{read_iadhore_output}}.
#' @param genome_x Character. First genome name.
#' @param genome_y Character. Second genome name (default: same as
#'   \code{genome_x} for intra-genomic PARs).
#' @param p_cutoff Numeric. Poisson p-value threshold for significant segment
#'   pairs (default \code{1e-10}).
#' @param r_cutoff Numeric. Minimum Pearson correlation to merge two segments
#'   into the same cluster (default \code{0.3}).
#' @param real_only Logical. Use only real anchorpoints (default \code{TRUE}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{pairs}}{Data frame of significant segment pairs with PAR
#'     assignment, segment IDs, and genomic coordinate ranges.  Columns:
#'     \code{par_id}, \code{seg_x_id}, \code{list_x}, \code{x_min},
#'     \code{x_max}, \code{seg_y_id}, \code{list_y}, \code{y_min},
#'     \code{y_max}, \code{score}.}
#'   \item{\code{score_matrix}}{Full \eqn{-\log_{10}(p)} matrix (all segments).}
#'   \item{\code{score_matrix_sig}}{Filtered matrix (only significant segments).}
#'   \item{\code{hclust_x}, \code{hclust_y}}{\code{hclust} objects used for
#'     clustering genome\_x and genome\_y segments respectively.}
#'   \item{\code{clusters_x}, \code{clusters_y}}{Named integer vectors mapping
#'     segment ID to PAR cluster assignment.}
#'   \item{\code{seg_bounds_x}, \code{seg_bounds_y}}{Data frames with columns
#'     \code{segment}, \code{list}, \code{seg_min}, \code{seg_max} for the
#'     significant segments.}
#'   \item{\code{r_cutoff}}{The r cutoff used for clustering.}
#'   \item{\code{n_pars}}{Integer. Number of PARs identified.}
#'   \item{\code{genome_x}, \code{genome_y}}{Genome names used.}
#' }
#' @export
#' @examples
#' \dontrun{
#' res  <- read_iadhore_output("output/ath_vvi/")
#' pars <- find_pars(res, genome_x = "ath", genome_y = "vvi")
#' plot_dotplot(res, genome_x = "ath", genome_y = "vvi",
#'              chr_order = "natural", highlight_pars = pars)
#' plot_par_overview(pars)
#' }
find_pars <- function(output_dir,
                       genome_x,
                       genome_y  = genome_x,
                       p_cutoff  = 1e-10,
                       r_cutoff  = 0.3,
                       real_only = TRUE) {

  dat   <- .load_iadhore(output_dir)
  mults <- dat$multiplicons
  ap    <- dat$anchorpoints
  segs  <- dat$segments
  le    <- dat$list_elements
  genes <- dat$genes

  if (real_only) ap <- ap[ap$is_real_anchorpoint != 0, , drop = FALSE]

  # Attach genome/list info to anchorpoints from multiplicons
  ap <- merge(ap,
              mults[, c("id", "genome_x", "list_x", "genome_y", "list_y", "level")],
              by.x = "multiplicon", by.y = "id",
              all.x = TRUE, sort = FALSE)

  # Normalise direction so genome_x/gene_x always correspond to genome_x arg
  intra   <- genome_x == genome_y
  fwd     <- !is.na(ap$genome_x) & ap$genome_x == genome_x & ap$genome_y == genome_y
  rev_dir <- !is.na(ap$genome_x) & ap$genome_x == genome_y & ap$genome_y == genome_x & !intra

  if (any(rev_dir)) {
    ar <- ap[rev_dir, , drop = FALSE]
    for (pair in list(c("gene_x", "gene_y"), c("coord_x", "coord_y"),
                      c("genome_x", "genome_y"), c("list_x", "list_y"))) {
      tmp <- ar[[pair[1L]]]; ar[[pair[1L]]] <- ar[[pair[2L]]]; ar[[pair[2L]]] <- tmp
    }
    ap <- rbind(ap[fwd, , drop = FALSE], ar)
  } else {
    ap <- ap[fwd, , drop = FALSE]
  }

  if (nrow(ap) == 0L) {
    message("No anchorpoints found for ", genome_x, " vs ", genome_y)
    return(invisible(NULL))
  }

  # Segments for each genome
  segs_x <- segs[segs$genome == genome_x, , drop = FALSE]
  segs_y <- if (intra) segs_x else segs[segs$genome == genome_y, , drop = FALSE]

  # Build gene -> segment map: keep lowest-level multiplicon assignment per gene
  .gene_seg_map <- function(seg_sub, le_tbl, mults_df) {
    le_s <- le_tbl[le_tbl$segment %in% seg_sub$id, c("gene", "segment")]
    if (nrow(le_s) == 0L) return(le_s)
    seg_lev <- merge(seg_sub[, c("id", "multiplicon")],
                     mults_df[, c("id", "level")],
                     by.x = "multiplicon", by.y = "id")
    le_s <- merge(le_s, seg_lev[, c("id", "level")],
                  by.x = "segment", by.y = "id", all.x = TRUE)
    le_s <- le_s[order(le_s$gene, le_s$level, na.last = TRUE), ]
    le_s[!duplicated(le_s$gene), c("gene", "segment")]
  }

  gsm_x <- .gene_seg_map(segs_x, le, mults)
  gsm_y <- if (intra) gsm_x else .gene_seg_map(segs_y, le, mults)

  ap$seg_x <- gsm_x$segment[match(ap$gene_x, gsm_x$gene)]
  ap$seg_y <- gsm_y$segment[match(ap$gene_y, gsm_y$gene)]
  ap <- ap[!is.na(ap$seg_x) & !is.na(ap$seg_y), , drop = FALSE]

  if (nrow(ap) == 0L) {
    message("Could not map anchorpoints to segments (check list_elements table)")
    return(invisible(NULL))
  }

  seg_x_ids <- sort(unique(as.integer(segs_x$id)))
  seg_y_ids <- sort(unique(as.integer(segs_y$id)))

  # Anchorpoint count matrix
  ct <- table(factor(ap$seg_x, levels = seg_x_ids),
              factor(ap$seg_y, levels = seg_y_ids))
  counts <- matrix(as.integer(ct),
                   nrow = length(seg_x_ids), ncol = length(seg_y_ids),
                   dimnames = list(as.character(seg_x_ids),
                                   as.character(seg_y_ids)))

  if (intra) counts[lower.tri(counts, diag = TRUE)] <- 0L

  # Gene counts per segment (Poisson null)
  n_x <- tabulate(match(gsm_x$segment, seg_x_ids), nbins = length(seg_x_ids))
  n_y <- if (intra) n_x else
         tabulate(match(gsm_y$segment, seg_y_ids), nbins = length(seg_y_ids))
  n_x[n_x == 0L] <- 1L; n_y[n_y == 0L] <- 1L
  N_total  <- nrow(ap)
  expected <- outer(n_x / sum(n_x), n_y / sum(n_y)) * N_total

  # -log10(Poisson p-value) score matrix
  score_mat <- matrix(0.0, nrow(counts), ncol(counts), dimnames = dimnames(counts))
  pos <- counts > 0L
  if (any(pos)) {
    pp          <- stats::ppois(counts[pos] - 1L, lambda = expected[pos],
                                lower.tail = FALSE)
    pp[pp <= 0] <- .Machine$double.xmin
    score_mat[pos] <- -log10(pp)
  }

  sig_thr <- -log10(p_cutoff)

  has_sig_r <- apply(score_mat, 1L, max) >= sig_thr
  has_sig_c <- apply(score_mat, 2L, max) >= sig_thr
  sm        <- score_mat[has_sig_r, has_sig_c, drop = FALSE]

  if (nrow(sm) == 0L || ncol(sm) == 0L) {
    mx <- round(max(score_mat), 1L)
    message("No segment pairs significant at p_cutoff = ", p_cutoff,
            " (max score = ", mx, ")")
    return(invisible(list(pairs           = .empty_pairs(),
                           score_matrix    = score_mat,
                           score_matrix_sig = sm,
                           n_pars          = 0L,
                           genome_x        = genome_x,
                           genome_y        = genome_y)))
  }

  # Hierarchical clustering — return both cluster assignments AND hclust object
  .clust_segs <- function(mat) {
    if (nrow(mat) < 2L) {
      return(list(clusters = stats::setNames(1L, rownames(mat)), hclust = NULL))
    }
    cr            <- stats::cor(t(mat), method = "pearson")
    cr[is.na(cr)] <- 0
    diss          <- stats::as.dist((1 - cr) / 2)
    hc            <- stats::hclust(diss, method = "average")
    cl            <- stats::cutree(hc, h = (1 - r_cutoff) / 2)
    list(clusters = cl, hclust = hc)
  }

  res_x <- .clust_segs(sm)
  res_y <- .clust_segs(t(sm))
  cl_x  <- res_x$clusters;  hc_x <- res_x$hclust
  cl_y  <- res_y$clusters;  hc_y <- res_y$hclust

  # Segment coordinate bounds for visualisation
  .seg_bounds <- function(seg_ids_chr, gname) {
    le_g <- le[le$segment %in% as.integer(seg_ids_chr), c("segment", "gene")]
    g    <- genes[genes$genome == gname, c("id", "list", "remapped_coordinate")]
    m    <- merge(le_g, g, by.x = "gene", by.y = "id", all.x = TRUE)
    m    <- m[!is.na(m$remapped_coordinate), , drop = FALSE]
    if (nrow(m) == 0L)
      return(data.frame(segment = integer(), list = character(),
                        seg_min = integer(), seg_max = integer(),
                        stringsAsFactors = FALSE))
    do.call(rbind, lapply(sort(unique(m$segment)), function(s) {
      df <- m[m$segment == s, , drop = FALSE]
      data.frame(segment = s, list = df$list[1L],
                 seg_min = min(df$remapped_coordinate),
                 seg_max = max(df$remapped_coordinate),
                 stringsAsFactors = FALSE)
    }))
  }

  bounds_x <- .seg_bounds(rownames(sm), genome_x)
  bounds_y <- .seg_bounds(colnames(sm), if (intra) genome_x else genome_y)

  # Enumerate significant (seg_x, seg_y) pairs and assign PAR IDs
  pairs_rows <- list()
  par_id     <- 1L

  for (cxi in sort(unique(cl_x))) {
    sx <- names(cl_x)[cl_x == cxi]
    for (cyi in sort(unique(cl_y))) {
      sy  <- names(cl_y)[cl_y == cyi]
      sub <- sm[sx, sy, drop = FALSE]
      sig <- which(sub >= sig_thr, arr.ind = TRUE)
      if (nrow(sig) == 0L) next

      for (k in seq_len(nrow(sig))) {
        sxi <- sx[sig[k, 1L]]; syi <- sy[sig[k, 2L]]
        bx  <- bounds_x[bounds_x$segment == as.integer(sxi), , drop = FALSE]
        by_ <- bounds_y[bounds_y$segment == as.integer(syi), , drop = FALSE]
        if (nrow(bx) == 0L || nrow(by_) == 0L) next
        pairs_rows[[length(pairs_rows) + 1L]] <- data.frame(
          par_id   = par_id,
          seg_x_id = as.integer(sxi),
          list_x   = bx$list[1L],
          x_min    = bx$seg_min[1L],
          x_max    = bx$seg_max[1L],
          seg_y_id = as.integer(syi),
          list_y   = by_$list[1L],
          y_min    = by_$seg_min[1L],
          y_max    = by_$seg_max[1L],
          score    = sub[sig[k, 1L], sig[k, 2L]],
          stringsAsFactors = FALSE
        )
      }
      par_id <- par_id + 1L
    }
  }

  n_pars   <- par_id - 1L
  pairs_df <- if (length(pairs_rows) > 0L) do.call(rbind, pairs_rows) else .empty_pairs()

  message("Found ", n_pars, " PARs covering ",
          nrow(pairs_df), " significant segment pairs")

  list(pairs            = pairs_df,
       score_matrix     = score_mat,
       score_matrix_sig = sm,
       hclust_x         = hc_x,
       hclust_y         = hc_y,
       clusters_x       = cl_x,
       clusters_y       = cl_y,
       seg_bounds_x     = bounds_x,
       seg_bounds_y     = bounds_y,
       r_cutoff         = r_cutoff,
       n_pars           = n_pars,
       genome_x         = genome_x,
       genome_y         = genome_y)
}


# Internal: empty pairs data frame with correct column types
.empty_pairs <- function() {
  data.frame(par_id   = integer(), seg_x_id = integer(),
             list_x   = character(), x_min = integer(), x_max = integer(),
             seg_y_id = integer(),
             list_y   = character(), y_min = integer(), y_max = integer(),
             score    = numeric(),
             stringsAsFactors = FALSE)
}
