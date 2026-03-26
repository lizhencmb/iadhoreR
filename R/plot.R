#' Draw a whole-genome synteny dot plot from i-ADHoRe anchorpoints
#'
#' Draws a single dot plot with all chromosomes concatenated along each axis.
#' Each anchorpoint is plotted as a dot at its genome-wide position.
#' Chromosome boundaries are shown as dashed lines with chromosome name labels.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   already returned by \code{\link{read_iadhore_output}}.
#' @param genome_x,genome_y Character. Which genome to place on each axis.
#'   Defaults to the first genome found (for intra-genomic runs both axes are
#'   the same genome).  For two-genome runs supply both names.
#' @param chr_order Character, one of \code{"size"} (default, largest
#'   chromosome first), \code{"alpha"} (alphabetical), or \code{"natural"}
#'   (sort by the numeric part of the name, e.g. chr1 < chr2 < chr10).
#' @param real_only Logical.  If \code{TRUE} (default), only real anchorpoints
#'   (\code{is_real_anchorpoint != 0}) are shown.
#' @param colour_by Character, one of \code{"level"} (default),
#'   \code{"multiplicon"}, or \code{"basecluster"}.  \code{"level"} colours
#'   points by the multiplicon level using a sequential palette (light = low,
#'   dark = high).
#' @param point_size Numeric. Size of points (default \code{0.8}).
#' @param highlight_pars List returned by \code{\link{find_pars}}, or
#'   \code{NULL} (default).  When supplied, each significant segment pair is
#'   drawn as a semi-transparent coloured rectangle behind the anchorpoint
#'   dots, with one colour per PAR.
#'
#' @return Invisibly returns the data frame of points that were plotted.
#' @export
#' @examples
#' \dontrun{
#' plot_dotplot("output/vvi/")
#' plot_dotplot("output/two_sp/", genome_x = "sp1", genome_y = "sp2",
#'              chr_order = "alpha")
#' res  <- read_iadhore_output("output/ath_vvi/")
#' pars <- find_pars(res, "ath", "vvi")
#' plot_dotplot(res, genome_x = "ath", genome_y = "vvi",
#'              chr_order = "natural", highlight_pars = pars)
#' }
plot_dotplot <- function(output_dir,
                         genome_x       = NULL,
                         genome_y       = NULL,
                         chr_order      = c("size", "natural", "alpha"),
                         real_only      = TRUE,
                         colour_by      = c("level", "multiplicon", "basecluster"),
                         point_size     = 0.8,
                         highlight_pars = NULL) {

  chr_order <- match.arg(chr_order)
  colour_by <- match.arg(colour_by)

  dat          <- .load_iadhore(output_dir)
  mults        <- dat$multiplicons
  anchorpoints <- dat$anchorpoints
  genes        <- dat$genes

  # Resolve genomes
  all_genomes <- sort(unique(genes$genome))
  if (is.null(genome_x)) genome_x <- all_genomes[1]
  if (is.null(genome_y)) genome_y <- all_genomes[if (length(all_genomes) >= 2L) 2L else 1L]

  # Helper: ordered chromosome table with genome-wide offsets for one genome
  chr_layout <- function(gname) {
    g  <- genes[genes$genome == gname, ]
    sz <- aggregate(remapped_coordinate ~ list, data = g,
                    FUN = function(x) max(x) + 1L)
    names(sz)[2] <- "size"
    sz <- if (chr_order == "alpha") {
      sz[order(sz$list), ]
    } else if (chr_order == "natural") {
      num <- suppressWarnings(as.integer(gsub("[^0-9]", "", sz$list)))
      sz[order(num, sz$list), ]
    } else {
      sz[order(-sz$size), ]
    }
    sz$offset <- cumsum(c(0L, sz$size[-nrow(sz)]))
    sz
  }

  lx <- chr_layout(genome_x)
  ly <- chr_layout(genome_y)

  total_x <- sum(lx$size)
  total_y <- sum(ly$size)

  # Attach list info to anchorpoints
  ap <- merge(anchorpoints,
              mults[, c("id", "genome_x", "list_x", "genome_y", "list_y", "level")],
              by.x = "multiplicon", by.y = "id",
              all.x = TRUE, sort = FALSE)

  # Keep only anchorpoints that belong to the selected genome pair.
  # i-ADHoRe may store either direction, so accept both orientations and
  # normalise reversed rows so genome_x/list_x always correspond to genome_x arg.
  fwd <- !is.na(ap$genome_x) &
           ap$genome_x == genome_x & ap$genome_y == genome_y
  rev_dir <- !is.na(ap$genome_x) &
               ap$genome_x == genome_y & ap$genome_y == genome_x &
               genome_x != genome_y   # don't double-count intra-genomic
  if (any(rev_dir)) {
    ap_rev <- ap[rev_dir, , drop = FALSE]
    # Swap x/y columns
    tmp_genome          <- ap_rev$genome_x; ap_rev$genome_x <- ap_rev$genome_y;  ap_rev$genome_y <- tmp_genome
    tmp_list            <- ap_rev$list_x;   ap_rev$list_x   <- ap_rev$list_y;    ap_rev$list_y   <- tmp_list
    tmp_gene            <- ap_rev$gene_x;   ap_rev$gene_x   <- ap_rev$gene_y;    ap_rev$gene_y   <- tmp_gene
    tmp_coord           <- ap_rev$coord_x;  ap_rev$coord_x  <- ap_rev$coord_y;   ap_rev$coord_y  <- tmp_coord
    ap <- rbind(ap[fwd, , drop = FALSE], ap_rev)
  } else {
    ap <- ap[fwd, , drop = FALSE]
  }

  if (real_only) ap <- ap[ap$is_real_anchorpoint != 0, , drop = FALSE]

  if (nrow(ap) == 0L) {
    n_ap   <- nrow(anchorpoints)
    n_real <- sum(anchorpoints$is_real_anchorpoint != 0, na.rm = TRUE)
    if (n_ap == 0L) {
      message("No anchorpoints found in output (anchorpoints.txt is empty).")
    } else if (real_only && n_real == 0L) {
      message("No real anchorpoints found (is_real_anchorpoint != 0 in 0 of ",
              n_ap, " rows). Try real_only = FALSE.")
    } else {
      message("No anchorpoints remain after filtering (total: ", n_ap,
              ", real: ", n_real, ").")
    }
    return(invisible(ap))
  }

  # Map local coordinates to genome-wide positions
  off_x <- setNames(lx$offset, lx$list)
  off_y <- setNames(ly$offset, ly$list)
  ap$gx <- ap$coord_x + off_x[ap$list_x]
  ap$gy <- ap$coord_y + off_y[ap$list_y]
  ap    <- ap[!is.na(ap$gx) & !is.na(ap$gy), , drop = FALSE]

  # For intra-genomic runs, mirror each point across the diagonal so the plot
  # is symmetric (i-ADHoRe only reports one direction per anchorpoint)
  if (genome_x == genome_y) {
    mirror  <- ap
    mirror$gx <- ap$gy
    mirror$gy <- ap$gx
    ap <- rbind(ap, mirror)
  }

  # Colour mapping
  if (colour_by == "level") {
    uniq_levels <- sort(unique(ap$level))
    n_lev       <- length(uniq_levels)
    pal <- if (n_lev == 1L) {
      "black"
    } else {
      grDevices::colorRampPalette(c("#4292c6", "#08306b"))(n_lev)  # dark blue range
    }
    ap$.col     <- pal[match(ap$level, uniq_levels)]
  } else {
    colour_ids <- ap[[colour_by]]
    uniq_ids   <- sort(unique(colour_ids))
    pal        <- .dotplot_palette(length(uniq_ids))
    ap$.col    <- pal[match(colour_ids, uniq_ids)]
  }

  # --- Plot ---
  old_par <- graphics::par(mar = c(5, 5, 3, 1), mgp = c(3.5, 0.5, 0))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, total_x), ylim = c(0, total_y))

  # Shaded alternating chromosome bands (subtle)
  for (i in seq_len(nrow(lx))) {
    if (i %% 2 == 0)
      graphics::rect(lx$offset[i], 0, lx$offset[i] + lx$size[i], total_y,
                     col = "#f5f5f5", border = NA)
  }
  for (i in seq_len(nrow(ly))) {
    if (i %% 2 == 0)
      graphics::rect(0, ly$offset[i], total_x, ly$offset[i] + ly$size[i],
                     col = "#f5f5f5", border = NA)
  }

  # Chromosome boundary lines
  for (i in seq_len(nrow(lx))[-1])
    graphics::abline(v = lx$offset[i], col = "grey70", lty = 2, lwd = 0.6)
  for (i in seq_len(nrow(ly))[-1])
    graphics::abline(h = ly$offset[i], col = "grey70", lty = 2, lwd = 0.6)

  # PAR highlight boxes: collect geometry first, draw fills behind dots,
  # then redraw borders on top so they are always visible.
  par_rects <- list()
  if (!is.null(highlight_pars) &&
      is.list(highlight_pars) &&
      !is.null(highlight_pars$pairs) &&
      nrow(highlight_pars$pairs) > 0L) {

    pr <- highlight_pars$pairs
    pr <- pr[pr$list_x %in% lx$list & pr$list_y %in% ly$list, , drop = FALSE]

    if (nrow(pr) > 0L) {
      n_par   <- max(pr$par_id)
      pal_par <- .dotplot_palette(n_par)

      for (i in seq_len(nrow(pr))) {
        col_i <- pal_par[pr$par_id[i]]
        gx1 <- pr$x_min[i] + off_x[pr$list_x[i]]
        gx2 <- pr$x_max[i] + off_x[pr$list_x[i]]
        gy1 <- pr$y_min[i] + off_y[pr$list_y[i]]
        gy2 <- pr$y_max[i] + off_y[pr$list_y[i]]
        par_rects[[length(par_rects) + 1L]] <- list(
          x1 = gx1, x2 = gx2, y1 = gy1, y2 = gy2, col = col_i)
        # Mirror for intra-genomic plots
        if (genome_x == genome_y) {
          par_rects[[length(par_rects) + 1L]] <- list(
            x1 = gy1, x2 = gy2, y1 = gx1, y2 = gx2, col = col_i)
        }
      }
    }
  }

  # Draw fills behind dots
  for (r in par_rects)
    graphics::rect(r$x1, r$y1, r$x2, r$y2,
                   col    = grDevices::adjustcolor(r$col, alpha.f = 0.18),
                   border = NA)

  # Dots
  graphics::points(ap$gx, ap$gy, col = ap$.col, pch = 20, cex = point_size)

  # Draw borders on top of dots
  for (r in par_rects)
    graphics::rect(r$x1, r$y1, r$x2, r$y2,
                   col    = NA,
                   border = grDevices::adjustcolor(r$col, alpha.f = 0.85),
                   lwd    = 1.0)

  graphics::box()

  # Chromosome name labels on axes (at midpoints, suppressing default axis)
  graphics::axis(1, at = lx$offset + lx$size / 2, labels = lx$list,
                 tick = FALSE, las = 2, cex.axis = 0.7, line = 0.5)
  graphics::axis(2, at = ly$offset + ly$size / 2, labels = ly$list,
                 tick = FALSE, las = 2, cex.axis = 0.7, line = 0.5)

  # Axis titles
  title_str <- if (genome_x == genome_y) paste0(genome_x, "  (intra-genomic)")
               else                      paste0(genome_x, " vs ", genome_y)
  graphics::title(main = title_str,
                  xlab = genome_x,
                  ylab = genome_y,
                  cex.main = 0.95)

  invisible(ap[, setdiff(names(ap), c(".col", "gx", "gy"))])
}


#' Plot a single multiplicon collinear block
#'
#' Draws all segments of one multiplicon as horizontal gene tracks.  Each gene
#' is shown as a filled (forward-strand) or open (reverse-strand) box.  Anchor
#' point pairs are connected by coloured semi-transparent bands; all anchors
#' belonging to the same basecluster share a colour.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   from \code{\link{read_iadhore_output}}.
#' @param multiplicon_id Integer. The multiplicon to draw.
#' @param real_only Logical. Use only real anchor points (default \code{TRUE}).
#' @param band_alpha Numeric (0–1). Transparency of the anchor bands
#'   (default \code{0.35}).
#' @param show_gene_ids Logical. If \code{TRUE}, print gene IDs below each
#'   arrow (rotated 90°). Useful for small multiplicons; can be crowded for
#'   large ones (default \code{FALSE}).
#' @param gene_id_cex Numeric. Character size for gene ID labels
#'   (default \code{0.4}).
#'
#' @return Invisibly returns a list with elements \code{segments},
#'   \code{list_elements}, and \code{anchorpoints}.
#' @export
#' @examples
#' \dontrun{
#' plot_multiplicon("output/arath/", multiplicon_id = 1)
#' plot_multiplicon("output/arath/", multiplicon_id = 1, show_gene_ids = TRUE)
#' }
plot_multiplicon <- function(output_dir,
                              multiplicon_id,
                              real_only     = TRUE,
                              band_alpha    = 0.35,
                              show_gene_ids = FALSE,
                              gene_id_cex   = 0.4) {

  dat  <- .load_iadhore(output_dir)
  segs <- dat$segments[dat$segments$multiplicon == multiplicon_id, ]
  if (nrow(segs) == 0L)
    stop("No segments found for multiplicon ", multiplicon_id)
  segs <- segs[order(segs$order), ]
  segs$track <- seq_len(nrow(segs))   # convert 0-based order to 1-based track index

  le <- dat$list_elements[dat$list_elements$segment %in% segs$id, ]
  le <- merge(le, segs[, c("id", "track", "genome", "list")],
              by.x = "segment", by.y = "id", sort = FALSE)

  ap <- dat$anchorpoints[dat$anchorpoints$multiplicon == multiplicon_id, ]
  if (real_only) ap <- ap[ap$is_real_anchorpoint != 0, ]

  n_segs     <- nrow(segs)
  x_range    <- range(le$position, na.rm = TRUE)
  pad        <- max(1L, diff(x_range) * 0.05)
  mult_level <- dat$multiplicons$level[dat$multiplicons$id == multiplicon_id]

  # Collect anchor points from this multiplicon AND all ancestors so that
  # homolog relationships at every level are represented
  all_mult_ids <- multiplicon_id
  mults        <- dat$multiplicons
  cur_id       <- multiplicon_id
  repeat {
    parent_val <- mults$parent[mults$id == cur_id]
    if (length(parent_val) == 0) break
    parent_id  <- suppressWarnings(as.integer(trimws(parent_val)))
    if (is.na(parent_id) || parent_id %in% all_mult_ids) break
    all_mult_ids <- c(all_mult_ids, parent_id)
    cur_id       <- parent_id
  }
  ap_all <- dat$anchorpoints[dat$anchorpoints$multiplicon %in% all_mult_ids, ]
  if (real_only) ap_all <- ap_all[ap_all$is_real_anchorpoint != 0, ]

  # Union-find: build connected components (homolog groups) across all levels
  all_genes <- unique(c(ap_all$gene_x, ap_all$gene_y))
  uf <- new.env(hash = TRUE, parent = emptyenv())
  for (g in all_genes) uf[[g]] <- g

  uf_find <- function(x) {
    while (uf[[x]] != x) {
      uf[[x]] <- uf[[uf[[x]]]]   # path compression
      x <- uf[[x]]
    }
    x
  }
  for (k in seq_len(nrow(ap_all))) {
    rx <- uf_find(ap_all$gene_x[k])
    ry <- uf_find(ap_all$gene_y[k])
    if (rx != ry) uf[[ry]] <- rx
  }
  for (g in all_genes) uf[[g]] <- uf_find(g)   # resolve all roots

  roots      <- vapply(all_genes, function(g) uf[[g]], character(1))
  uniq_roots <- unique(roots)
  pal        <- .dotplot_palette(length(uniq_roots))
  root_col   <- setNames(pal, uniq_roots)
  gene_col   <- setNames(root_col[roots], all_genes)

  bw <- 0.35; bh <- 0.18   # gene box half-width / half-height

  old_par <- graphics::par(mar = c(1, 9, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(xlim = x_range + c(-pad, pad),
                        ylim = c(0.5, n_segs + 0.5),
                        xaxs = "i", yaxs = "i")
  graphics::title(
    main = sprintf("Multiplicon %d  (level %d, %d anchor points)",
                   multiplicon_id, mult_level, nrow(ap_all)),
    cex.main = 0.9
  )

  # ── 1. Anchor bands (drawn first, behind gene boxes) ──────────────────────
  # Use ap_all so parent-level homolog pairs are also drawn
  for (k in seq_len(nrow(ap_all))) {
    rx <- le[le$gene == ap_all$gene_x[k], ]
    ry <- le[le$gene == ap_all$gene_y[k], ]
    if (nrow(rx) == 0L || nrow(ry) == 0L) next

    x1 <- rx$position[1]; y1 <- rx$track[1]
    x2 <- ry$position[1]; y2 <- ry$track[1]

    base_col  <- gene_col[ap_all$gene_x[k]]
    fill_col  <- grDevices::adjustcolor(base_col, alpha.f = band_alpha)
    edge_col  <- grDevices::adjustcolor(base_col, alpha.f = min(1, band_alpha + 0.3))

    # Band edges: match the exact bottom edge of each arrow polygon
    # Forward (+): body spans x-bw to x+bw*0.4
    # Reverse (-): body spans x-bw*0.4 to x+bw
    o1  <- rx$orientation[1]
    o2  <- ry$orientation[1]
    x1l <- if (!is.na(o1) && o1 == "-") x1 - bw * 0.4 else x1 - bw
    x1r <- if (!is.na(o1) && o1 == "-") x1 + bw       else x1 + bw * 0.4
    x2l <- if (!is.na(o2) && o2 == "-") x2 - bw * 0.4 else x2 - bw
    x2r <- if (!is.na(o2) && o2 == "-") x2 + bw       else x2 + bw * 0.4

    y1e <- if (y1 > y2) y1 - bh else y1 + bh
    y2e <- if (y2 > y1) y2 - bh else y2 + bh

    graphics::polygon(
      x      = c(x1l, x1r, x2r, x2l),
      y      = c(y1e, y1e, y2e, y2e),
      col    = fill_col,
      border = edge_col,
      lwd    = 0.6
    )
  }

  # ── 2. Gene tracks (drawn on top of bands) ────────────────────────────────
  for (i in seq_len(n_segs)) {
    y  <- segs$track[i]
    ge <- le[le$segment == segs$id[i], ]
    if (nrow(ge) == 0L) next

    # Chromosome baseline
    x_lo <- min(ge$position, na.rm = TRUE)
    x_hi <- max(ge$position, na.rm = TRUE)
    graphics::segments(x_lo, y, x_hi, y, lwd = 1.5, col = "grey50")

    # Gene arrows: right-pointing = "+" (forward), left-pointing = "-" (reverse)
    # Anchor genes share their basecluster colour; non-anchor genes are grey
    for (j in seq_len(nrow(ge))) {
      xj    <- ge$position[j]
      rev   <- !is.na(ge$orientation[j]) && ge$orientation[j] == "-"
      cfill <- gene_col[ge$gene[j]]
      if (is.na(cfill)) cfill <- "grey80"
      if (rev) {
        px <- c(xj + bw,      xj - bw * 0.4, xj - bw, xj - bw * 0.4, xj + bw)
      } else {
        px <- c(xj - bw,      xj + bw * 0.4, xj + bw, xj + bw * 0.4, xj - bw)
      }
      py <- c(y - bh, y - bh, y, y + bh, y + bh)
      graphics::polygon(px, py,
                        col    = cfill,
                        border = "grey30", lwd = 0.4)
      if (show_gene_ids) {
        graphics::text(xj, y - bh - 0.02, labels = ge$gene[j],
                       srt = 90, adj = c(1, 0.5),
                       cex = gene_id_cex, col = "grey20")
      }
    }

    lbl <- paste0(segs$genome[i], "\n", segs$list[i])
    graphics::mtext(lbl, side = 2, at = y, las = 1, cex = 0.65, line = 0.5)
  }


  invisible(list(segments      = segs,
                 list_elements  = le,
                 anchorpoints   = ap))
}


#' Genome-wide overview of multiplicon segment positions
#'
#' Draws each chromosome as a thick black baseline and shows non-redundant
#' multiplicon segments as coloured rectangles stacked above the baseline, in
#' the style of Fig. 2 of Simillion et al. (2002).
#'
#' Colour scheme (matching the original paper):
#' \itemize{
#'   \item Level 2 (single duplication): light grey
#'   \item Levels 3–4 (two rounds): dark grey
#'   \item Level \eqn{\geq 5} (three or more rounds): distinct colour per
#'     multiplicon
#' }
#'
#' Rectangles are stacked when segments overlap on the same chromosome, so the
#' stack height at any position reflects the local multiplication level.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   from \code{\link{read_iadhore_output}}.
#' @param genome Character.  Genome to plot.  Required when multiple genomes
#'   are present; inferred automatically for single-genome runs.
#' @param chr_order Character, one of \code{"natural"} (default, sort by the
#'   numeric part of the name), \code{"alpha"} (alphabetical), or \code{"size"}
#'   (largest chromosome first, top of plot).
#' @param seg_height Numeric.  Height of each stacked rectangle in gene-index
#'   units (default \code{0.3}).
#' @param band_gap Numeric.  Vertical gap between chromosome bands
#'   (default \code{1.5}).
#'
#' @return Invisibly returns the data frame of segments that were plotted.
#' @export
#' @examples
#' \dontrun{
#' plot_genome_overview("output/arath/")
#' plot_genome_overview("output/multi/", genome = "arath")
#' }
plot_genome_overview <- function(output_dir,
                                  genome     = NULL,
                                  chr_order  = c("natural", "alpha", "size"),
                                  seg_height = 0.3,
                                  band_gap   = 1.5) {

  chr_order <- match.arg(chr_order)

  dat         <- .load_iadhore(output_dir)
  genes        <- dat$genes
  multiplicons <- dat$multiplicons
  segments     <- dat$segments

  # Resolve genome
  all_genomes <- unique(genes$genome)
  if (is.null(genome)) {
    if (length(all_genomes) == 1L) {
      genome <- all_genomes
    } else {
      stop("Multiple genomes detected: ", paste(all_genomes, collapse = ", "),
           ". Please specify `genome`.")
    }
  }
  genes <- genes[genes$genome == genome, ]
  if (nrow(genes) == 0L) stop("No genes found for genome '", genome, "'")

  # Chromosome sizes and sorted order
  chrom_sizes <- aggregate(remapped_coordinate ~ list, data = genes,
                            FUN = function(x) max(x) + 1L)
  names(chrom_sizes)[2] <- "size"
  chroms <- if (chr_order == "alpha") {
    sort(chrom_sizes$list)
  } else if (chr_order == "natural") {
    num <- suppressWarnings(as.integer(gsub("[^0-9]", "", chrom_sizes$list)))
    chrom_sizes$list[order(num, chrom_sizes$list)]
  } else {
    chrom_sizes$list[order(-chrom_sizes$size)]
  }
  n_chr  <- length(chroms)

  # Non-redundant segments on this genome
  nr_ids <- multiplicons$id[multiplicons$is_redundant == 0]
  seg_nr <- segments[segments$multiplicon %in% nr_ids &
                       segments$genome == genome, , drop = FALSE]

  # Resolve first / last gene → remapped coordinate
  gene_coords <- genes[, c("id", "remapped_coordinate")]
  seg_nr <- merge(seg_nr, gene_coords, by.x = "first", by.y = "id",
                  all.x = TRUE, sort = FALSE)
  names(seg_nr)[names(seg_nr) == "remapped_coordinate"] <- "first_coord"
  seg_nr <- merge(seg_nr, gene_coords, by.x = "last",  by.y = "id",
                  all.x = TRUE, sort = FALSE)
  names(seg_nr)[names(seg_nr) == "remapped_coordinate"] <- "last_coord"

  # Ensure first_coord <= last_coord
  swap <- !is.na(seg_nr$first_coord) & !is.na(seg_nr$last_coord) &
            seg_nr$first_coord > seg_nr$last_coord
  tmp                     <- seg_nr$first_coord[swap]
  seg_nr$first_coord[swap] <- seg_nr$last_coord[swap]
  seg_nr$last_coord[swap]  <- tmp

  # Drop rows with missing coordinates
  seg_nr <- seg_nr[!is.na(seg_nr$first_coord) & !is.na(seg_nr$last_coord), ]

  # Attach multiplicon level
  seg_nr <- merge(seg_nr,
                  multiplicons[, c("id", "level")],
                  by.x = "multiplicon", by.y = "id",
                  all.x = TRUE, sort = FALSE)

  # Colour by level: 2 = light grey, 3-4 = dark grey, >=5 = distinct
  high_mults <- sort(unique(seg_nr$multiplicon[!is.na(seg_nr$level) &
                                                   seg_nr$level >= 5]))
  pal_high   <- .dotplot_palette(length(high_mults))
  names(pal_high) <- as.character(high_mults)

  seg_nr$.col <- ifelse(
    is.na(seg_nr$level), "grey80",
    ifelse(seg_nr$level == 2, "grey85",
    ifelse(seg_nr$level <= 4, "grey40",
           pal_high[as.character(seg_nr$multiplicon)])))

  # Greedy interval stacking per chromosome.
  # Sort by level first so level 2 fills the bottom stacks, then level 3, etc.
  seg_nr$.stack <- NA_integer_
  for (ch in chroms) {
    idx <- which(seg_nr$list == ch)
    if (length(idx) == 0L) next
    ord <- order(seg_nr$level[idx], seg_nr$first_coord[idx], na.last = TRUE)
    idx <- idx[ord]
    ends <- numeric(0)
    for (j in seq_along(idx)) {
      fc <- seg_nr$first_coord[idx[j]]
      lc <- seg_nr$last_coord[idx[j]]
      placed <- FALSE
      for (lv in seq_along(ends)) {
        if (fc > ends[lv]) {
          seg_nr$.stack[idx[j]] <- lv
          ends[lv] <- lc
          placed <- TRUE
          break
        }
      }
      if (!placed) {
        ends <- c(ends, lc)
        seg_nr$.stack[idx[j]] <- length(ends)
      }
    }
  }

  max_stack  <- max(seg_nr$.stack, na.rm = TRUE)
  x_max      <- max(chrom_sizes$size, na.rm = TRUE)
  band_h     <- max_stack * seg_height + band_gap

  # y baseline for each chromosome (top chromosome drawn highest)
  y_base <- setNames(
    (n_chr - seq_len(n_chr)) * band_h,
    chroms
  )
  tot_h <- n_chr * band_h

  old_par <- graphics::par(mar = c(3.5, 6, 2.5, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, x_max), ylim = c(-band_gap * 0.5, tot_h))
  graphics::title(
    main = paste("Multiplicon overview \u2013", genome),
    xlab = "Gene index (remapped coordinate)",
    cex.main = 0.9
  )
  graphics::axis(1, cex.axis = 0.8)

  for (ch in chroms) {
    yb  <- y_base[ch]
    csz <- chrom_sizes$size[chrom_sizes$list == ch]

    # Chromosome baseline
    graphics::segments(0, yb, csz, yb, lwd = 2.5, col = "black")
    graphics::mtext(ch, side = 2, at = yb, las = 1, cex = 0.75, line = 0.5)

    # Segments stacked above baseline — draw level 2 first, then 3, 4, …
    idx <- which(seg_nr$list == ch & !is.na(seg_nr$.stack))
    idx <- idx[order(seg_nr$level[idx], na.last = TRUE)]
    for (i in idx) {
      stk <- seg_nr$.stack[i]
      y0  <- yb + (stk - 1L) * seg_height + seg_height * 0.1
      y1  <- yb + stk        * seg_height - seg_height * 0.1
      graphics::rect(seg_nr$first_coord[i], y0,
                     seg_nr$last_coord[i],  y1,
                     col    = seg_nr$.col[i],
                     border = NA)
    }
  }

  # Legend
  leg_cols <- c("grey85", "grey40")
  leg_labs <- c("level 2", "levels 3\u20134")
  if (length(high_mults) > 0L) {
    leg_cols <- c(leg_cols, pal_high[[1L]])
    leg_labs <- c(leg_labs, "level \u22655")
  }
  graphics::legend("topright", legend = leg_labs, fill = leg_cols,
                   border = NA, bty = "n", cex = 0.75)

  invisible(seg_nr[, setdiff(names(seg_nr), c(".col", ".stack"))])
}


#' Draw a clustered segment dot plot with dendrogram (Fig. S1 style)
#'
#' Reproduces the overview plot from Jaillon \emph{et al.} (2009) Fig. S1.
#' All significant segments from both genomes are placed on the two axes,
#' reordered by the same average-linkage hierarchical clustering used in
#' \code{\link{find_pars}}.  Anchorpoint dots are drawn in the main panel;
#' dendrograms are drawn in side strips; PAR blocks are highlighted in gold.
#'
#' @param par_result List returned by \code{\link{find_pars}}.
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   from \code{\link{read_iadhore_output}}.  Required to fetch anchorpoints.
#' @param real_only Logical. Use only real anchorpoints (default \code{TRUE}).
#' @param point_size Numeric. Dot size (default \code{0.4}).
#' @param show_dendro Logical. Draw dendrograms (default \code{TRUE}).
#' @param par_col Character. Colour for PAR highlight boxes (default
#'   \code{"gold"}).
#' @param par_alpha Numeric (0–1). Transparency of PAR boxes (default
#'   \code{0.4}).
#' @param label_cex Numeric. Character size for axis labels (default
#'   \code{0.55}).
#'
#' @return Invisibly returns \code{par_result}.
#' @export
#' @examples
#' \dontrun{
#' res  <- read_iadhore_output("output/ath_vvi/")
#' pars <- find_pars(res, "ath", "vvi")
#' plot_par_overview(pars, res)
#' }
plot_par_overview <- function(par_result,
                               output_dir,
                               real_only  = TRUE,
                               point_size = 0.4,
                               show_dendro = TRUE,
                               par_col    = "gold",
                               par_alpha  = 0.4,
                               label_cex  = 0.55) {

  if (is.null(par_result) || par_result$n_pars == 0L) {
    message("No PARs to plot"); return(invisible(NULL))
  }

  sm       <- par_result$score_matrix_sig   # filtered score matrix (rows=x, cols=y)
  hc_x     <- par_result$hclust_x
  hc_y     <- par_result$hclust_y
  cl_x     <- par_result$clusters_x        # named: seg_id -> cluster
  cl_y     <- par_result$clusters_y
  bx       <- par_result$seg_bounds_x      # segment, list, seg_min, seg_max
  by_      <- par_result$seg_bounds_y
  pairs    <- par_result$pairs
  r_cutoff <- par_result$r_cutoff
  genome_x <- par_result$genome_x
  genome_y <- par_result$genome_y
  intra    <- genome_x == genome_y

  nx <- nrow(sm);  ny <- ncol(sm)

  # ---- Leaf ordering (display position for each row/col of sm) ----
  # hc$order[j] = which row of sm is at display position j
  # We need: for row i, what is its display position?  -> order(hc$order)[i]
  if (!is.null(hc_x)) {
    ord_x <- hc_x$order                   # permutation: display j -> sm row hc$order[j]
    rank_x <- order(ord_x)                # rank_x[i] = display position of sm row i
  } else {
    ord_x <- seq_len(nx); rank_x <- seq_len(nx)
  }
  if (!is.null(hc_y)) {
    ord_y <- hc_y$order
    rank_y <- order(ord_y)
  } else {
    ord_y <- seq_len(ny); rank_y <- seq_len(ny)
  }

  # Segment IDs (character) in sm order
  seg_ids_x <- rownames(sm)   # character, length nx
  seg_ids_y <- colnames(sm)

  # ---- Load & filter anchorpoints ----
  dat   <- .load_iadhore(output_dir)
  ap    <- dat$anchorpoints
  mults <- dat$multiplicons
  if (real_only) ap <- ap[ap$is_real_anchorpoint != 0, , drop = FALSE]
  ap <- merge(ap,
              mults[, c("id", "genome_x", "list_x", "genome_y", "list_y")],
              by.x = "multiplicon", by.y = "id", all.x = TRUE, sort = FALSE)
  fwd     <- !is.na(ap$genome_x) & ap$genome_x == genome_x & ap$genome_y == genome_y
  rev_dir <- !is.na(ap$genome_x) & ap$genome_x == genome_y &
             ap$genome_y == genome_x & !intra
  if (any(rev_dir)) {
    ar <- ap[rev_dir, , drop = FALSE]
    for (p in list(c("gene_x","gene_y"), c("coord_x","coord_y"),
                   c("genome_x","genome_y"), c("list_x","list_y"))) {
      tmp <- ar[[p[1L]]]; ar[[p[1L]]] <- ar[[p[2L]]]; ar[[p[2L]]] <- tmp
    }
    ap <- rbind(ap[fwd, , drop = FALSE], ar)
  } else {
    ap <- ap[fwd, , drop = FALSE]
  }
  if (intra) {   # mirror
    ap_m <- ap
    for (p in list(c("list_x","list_y"), c("coord_x","coord_y"), c("gene_x","gene_y"))) {
      tmp <- ap_m[[p[1L]]]; ap_m[[p[1L]]] <- ap_m[[p[2L]]]; ap_m[[p[2L]]] <- tmp
    }
    ap <- rbind(ap, ap_m)
  }

  # ---- Map anchorpoints -> display coordinates ----
  # seg_id -> row index in sm  (and then -> display position via rank_x/rank_y)
  sid_to_row_x <- stats::setNames(seq_len(nx), seg_ids_x)
  sid_to_row_y <- stats::setNames(seq_len(ny), seg_ids_y)

  # For each anchorpoint, find which x-segment (by coordinate range)
  ap_xi <- integer(nrow(ap))   # sm row index (1-based), 0 = not found
  for (i in seq_len(nrow(bx))) {
    hit <- ap$list_x == bx$list[i] &
           ap$coord_x >= bx$seg_min[i] & ap$coord_x <= bx$seg_max[i]
    ap_xi[hit] <- sid_to_row_x[as.character(bx$segment[i])]
  }
  ap_yi <- integer(nrow(ap))
  for (j in seq_len(nrow(by_))) {
    hit <- ap$list_y == by_$list[j] &
           ap$coord_y >= by_$seg_min[j] & ap$coord_y <= by_$seg_max[j]
    ap_yi[hit] <- sid_to_row_y[as.character(by_$segment[j])]
  }
  keep <- ap_xi > 0L & ap_yi > 0L
  xi_k <- ap_xi[keep]; yi_k <- ap_yi[keep]; ap_k <- ap[keep, , drop = FALSE]

  # Display position of each anchorpoint (segment midpoint ± fractional offset)
  # Coordinate system: x,y in (0.5, n+0.5); segment at rank r centred at r.
  bx_min <- bx$seg_min; bx_max <- bx$seg_max; bx_seg <- bx$segment
  by_min <- by_$seg_min; by_max <- by_$seg_max; by_seg <- by_$segment

  frac_x <- (ap_k$coord_x - bx_min[xi_k]) /
             pmax(1L, bx_max[xi_k] - bx_min[xi_k])
  frac_y <- (ap_k$coord_y - by_min[yi_k]) /
             pmax(1L, by_max[yi_k] - by_min[yi_k])

  # ---- Layout ----
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (show_dendro && (!is.null(hc_x) || !is.null(hc_y))) {
    dw <- 0.13   # fraction of width  for left dendrogram
    dh <- 0.13   # fraction of height for top dendrogram
    graphics::layout(
      matrix(c(0L, 1L, 2L, 3L), nrow = 2L, byrow = TRUE),
      widths  = c(dw, 1 - dw),
      heights = c(dh, 1 - dh)
    )
  } else {
    graphics::layout(matrix(1L, 1L, 1L))
  }

  # ---- Dendrogram helper ----
  # Compute line segments for a hclust tree.
  # Returns list(segs, leaf_rank, max_h):
  #   segs     = data.frame(x0,y0,x1,y1) in (position 1..n, height) space
  #   leaf_rank[i] = display position of leaf i (same as order(hc$order))
  .ds <- function(hc) {
    n    <- length(hc$order)
    lpos <- integer(n); lpos[hc$order] <- seq_len(n)
    npos <- numeric(n - 1L)
    mat  <- matrix(0.0, 3L * (n - 1L), 4L,
                   dimnames = list(NULL, c("x0","y0","x1","y1")))
    row  <- 1L
    .xp <- function(v) if (v < 0L) lpos[-v]      else npos[v]
    .yp <- function(v) if (v < 0L) 0.0           else hc$height[v]
    for (k in seq_len(n - 1L)) {
      a <- hc$merge[k,1L]; b <- hc$merge[k,2L]
      xa <- .xp(a); ya <- .yp(a); xb <- .xp(b); yb <- .yp(b); yk <- hc$height[k]
      npos[k] <- (xa + xb) / 2
      mat[row,  ] <- c(xa, ya, xa, yk)
      mat[row+1,] <- c(xb, yb, xb, yk)
      mat[row+2,] <- c(xa, yk, xb, yk)
      row <- row + 3L
    }
    list(segs = as.data.frame(mat), leaf_rank = lpos, max_h = max(hc$height))
  }

  # Note: the top-left cell (value 0 in the layout matrix) is skipped
  # automatically by graphics::layout() — no plot.new() needed there.

  # Shared margins: the non-shared edges of each dendrogram panel must match
  # the corresponding edge of the main plot so both panels have identical
  # plot-area widths/heights and leaves line up with segments exactly.
  #   top dendro  : left & right must match main plot left & right  (0.5 and 4)
  #   left dendro : top & bottom must match main plot top & bottom  (0.5 and 4)
  mar_main   <- c(4, 0.5, 0.5, 4)   # bottom, left, top, right of main panel
  mar_top_d  <- c(0, mar_main[2L], 1, mar_main[4L])   # top dendro
  mar_left_d <- c(mar_main[1L], 1, mar_main[3L], 0)   # left dendro

  # Both dendrograms and the main plot use the same x/y coordinate space:
  #   x in (0.5, nx+0.5)  —  segment at rank j centred at j
  #   y in (0.5, ny+0.5)  —  segment at rank j centred at j
  # This makes dendrogram leaf positions (integers 1..n) coincide with
  # segment centres in the main panel.

  # ---- Panel 1: top dendrogram (genome_x) ----
  if (show_dendro && !is.null(hc_x)) {
    graphics::par(mar = mar_top_d)
    ds_x    <- .ds(hc_x)
    h_cut_x <- (1 - r_cutoff) / 2
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, nx + 0.5), ylim = c(0, ds_x$max_h * 1.15))
    graphics::segments(ds_x$segs$x0, ds_x$segs$y0,
                       ds_x$segs$x1, ds_x$segs$y1, col = "grey20", lwd = 0.7)
    graphics::abline(h = h_cut_x, lty = 2, col = "grey50", lwd = 0.6)
    graphics::mtext(genome_x, side = 3, line = 0, cex = 0.75, font = 2)
  } else if (show_dendro) {
    graphics::par(mar = mar_top_d); graphics::plot.new()
  }

  # ---- Panel 2: left dendrogram (genome_y, rotated 90°) ----
  if (show_dendro && !is.null(hc_y)) {
    graphics::par(mar = mar_left_d)
    ds_y    <- .ds(hc_y)
    h_cut_y <- (1 - r_cutoff) / 2
    graphics::plot.new()
    # x = depth (root left, leaves right), y = leaf position matching main ylim
    graphics::plot.window(xlim = c(ds_y$max_h * 1.15, 0), ylim = c(0.5, ny + 0.5))
    graphics::segments(ds_y$segs$y0, ds_y$segs$x0,
                       ds_y$segs$y1, ds_y$segs$x1, col = "grey20", lwd = 0.7)
    graphics::abline(v = h_cut_y, lty = 2, col = "grey50", lwd = 0.6)
    graphics::mtext(genome_y, side = 2, line = 0, cex = 0.75, font = 2, las = 0)
  } else if (show_dendro) {
    graphics::par(mar = mar_left_d); graphics::plot.new()
  }

  # ---- Panel 3 (or 1 without dendro): main dot plot ----
  # Coordinate system: x in (0.5, nx+0.5), y in (0.5, ny+0.5).
  # Segment at rank j is centred at j and spans [j-0.5, j+0.5].
  # Anchorpoint at fraction f within a segment → display = rank + f - 0.5.
  px <- rank_x[xi_k] + frac_x - 0.5
  py <- rank_y[yi_k] + frac_y - 0.5

  graphics::par(mar = mar_main)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, nx + 0.5), ylim = c(0.5, ny + 0.5))

  # Alternating-shade columns and rows
  for (i in seq_len(nx)) {
    if (i %% 2 == 0)
      graphics::rect(i - 0.5, 0.5, i + 0.5, ny + 0.5, col = "#f7f7f7", border = NA)
  }
  for (j in seq_len(ny)) {
    if (j %% 2 == 0)
      graphics::rect(0.5, j - 0.5, nx + 0.5, j + 0.5, col = "#f7f7f7", border = NA)
  }

  # PAR highlight fills (drawn first, behind dots)
  par_boxes <- if (nrow(pairs) > 0L) {
    lapply(sort(unique(pairs$par_id)), function(pid) {
      pr  <- pairs[pairs$par_id == pid, , drop = FALSE]
      rxi <- rank_x[sid_to_row_x[as.character(pr$seg_x_id)]]
      ryi <- rank_y[sid_to_row_y[as.character(pr$seg_y_id)]]
      list(pid = pid,
           x1  = min(rxi) - 0.5, x2 = max(rxi) + 0.5,
           y1  = min(ryi) - 0.5, y2 = max(ryi) + 0.5)
    })
  } else list()

  for (b in par_boxes)
    graphics::rect(b$x1, b$y1, b$x2, b$y2,
                   col    = grDevices::adjustcolor(par_col, alpha.f = par_alpha),
                   border = NA)

  # Segment separator lines at boundaries between segments (thin grey)
  for (i in seq_len(nx - 1L))
    graphics::abline(v = i + 0.5, col = "grey80", lty = 1, lwd = 0.3)
  for (j in seq_len(ny - 1L))
    graphics::abline(h = j + 0.5, col = "grey80", lty = 1, lwd = 0.3)

  # Cluster boundary lines (green, thicker) — drawn at cluster group edges
  .cluster_breaks <- function(clusters, ord) {
    cl_ordered <- clusters[ord]
    which(diff(cl_ordered) != 0)
  }
  if (!is.null(hc_x) && length(unique(cl_x)) > 1L) {
    brk_x <- .cluster_breaks(cl_x, ord_x)
    for (b in brk_x)
      graphics::abline(v = b, col = "green3", lty = 1, lwd = 1.0)
  }
  if (!is.null(hc_y) && length(unique(cl_y)) > 1L) {
    brk_y <- .cluster_breaks(cl_y, ord_y)
    for (b in brk_y)
      graphics::abline(h = b, col = "green3", lty = 1, lwd = 1.0)
  }

  # Dots
  if (length(px) > 0L)
    graphics::points(px, py, pch = 20, cex = point_size, col = "black")

  # PAR borders and ID labels drawn on top of dots so they are always visible
  for (b in par_boxes) {
    graphics::rect(b$x1, b$y1, b$x2, b$y2,
                   col    = NA,
                   border = grDevices::adjustcolor(par_col, alpha.f = 0.9),
                   lwd    = 1.2)
    graphics::text(b$x1 + 0.05, b$y2 - 0.05,
                   labels = sprintf("PAR%02d", b$pid),
                   adj = c(0, 1), cex = 0.35, col = "grey10", font = 2)
  }

  graphics::box()

  # Axes: genome-level chromosome labels (one label per chromosome, centred)
  .chr_axis <- function(bounds, rank_fn, total, side) {
    chr_uniq <- unique(bounds$list)
    for (ch in chr_uniq) {
      rows   <- bounds[bounds$list == ch, , drop = FALSE]
      ranks  <- rank_fn(as.character(rows$segment))
      mid    <- (min(ranks) + max(ranks)) / 2
      graphics::mtext(ch, side = side, at = mid,
                      las = if (side %in% c(1L,3L)) 2L else 2L,
                      cex = label_cex, line = 0.3)
    }
  }
  .chr_axis(bx,  function(ids) rank_x[sid_to_row_x[ids]], nx, side = 1L)
  .chr_axis(by_, function(ids) rank_y[sid_to_row_y[ids]], ny, side = 4L)

  invisible(par_result)
}


#' Draw fine-resolution dot plots for individual PARs (Fig. S2 style)
#'
#' For each requested PAR, draws a dot plot with all member segments
#' concatenated along each axis, separated by green lines, and labelled with
#' \code{chromosome:start-stop}.  The layout mirrors Fig. S2 of Tang \emph{et
#' al.} (2009) \emph{PNAS} 106:14712.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   already returned by \code{\link{read_iadhore_output}}.
#' @param par_result List returned by \code{\link{find_pars}}.
#' @param par_ids Integer vector of PAR IDs to draw (default: all).
#' @param real_only Logical. Plot only real anchorpoints (default \code{TRUE}).
#' @param point_size Numeric. Dot size (default \code{0.8}).
#' @param label_cex Numeric. Character size for segment axis labels
#'   (default \code{0.45}).
#'
#' @return Invisibly returns \code{par_result}.
#' @export
#' @examples
#' \dontrun{
#' res  <- read_iadhore_output("output/ath_vvi/")
#' pars <- find_pars(res, "ath", "vvi")
#' plot_par(res, pars)                 # all PARs, one per page
#' plot_par(res, pars, par_ids = 1:6)  # first 6 only
#' }
plot_par <- function(output_dir,
                      par_result,
                      par_ids    = NULL,
                      real_only  = TRUE,
                      point_size = 0.8,
                      label_cex  = 0.55) {

  if (is.null(par_result) || par_result$n_pars == 0L ||
      nrow(par_result$pairs) == 0L) {
    message("No PARs to plot"); return(invisible(NULL))
  }

  pairs    <- par_result$pairs
  genome_x <- par_result$genome_x
  genome_y <- par_result$genome_y
  intra    <- genome_x == genome_y

  dat   <- .load_iadhore(output_dir)
  ap    <- dat$anchorpoints
  mults <- dat$multiplicons

  if (real_only) ap <- ap[ap$is_real_anchorpoint != 0, , drop = FALSE]

  # Attach genome / list info from multiplicons
  ap <- merge(ap,
              mults[, c("id", "genome_x", "list_x", "genome_y", "list_y")],
              by.x = "multiplicon", by.y = "id",
              all.x = TRUE, sort = FALSE)

  # Normalise direction so genome_x cols always correspond to the genome_x arg
  fwd     <- !is.na(ap$genome_x) & ap$genome_x == genome_x & ap$genome_y == genome_y
  rev_dir <- !is.na(ap$genome_x) & ap$genome_x == genome_y &
             ap$genome_y == genome_x & !intra
  if (any(rev_dir)) {
    ar <- ap[rev_dir, , drop = FALSE]
    for (p in list(c("gene_x","gene_y"), c("coord_x","coord_y"),
                   c("genome_x","genome_y"), c("list_x","list_y"))) {
      tmp <- ar[[p[1L]]]; ar[[p[1L]]] <- ar[[p[2L]]]; ar[[p[2L]]] <- tmp
    }
    ap <- rbind(ap[fwd, , drop = FALSE], ar)
  } else {
    ap <- ap[fwd, , drop = FALSE]
  }

  # For intra-genomic add mirrored anchorpoints (both directions)
  if (intra) {
    ap_m <- ap
    for (p in list(c("list_x","list_y"), c("coord_x","coord_y"), c("gene_x","gene_y"))) {
      tmp <- ap_m[[p[1L]]]; ap_m[[p[1L]]] <- ap_m[[p[2L]]]; ap_m[[p[2L]]] <- tmp
    }
    ap <- rbind(ap, ap_m)
  }

  # Which PARs to draw
  all_ids <- sort(unique(pairs$par_id))
  if (is.null(par_ids)) par_ids <- all_ids
  par_ids <- intersect(as.integer(par_ids), all_ids)
  n_plot  <- length(par_ids)
  if (n_plot == 0L) { message("No valid PAR IDs"); return(invisible(NULL)) }

  old_par <- graphics::par(mfrow = c(1L, 1L),
                             mar   = c(6, 0.5, 2.5, 8))
  on.exit(graphics::par(old_par), add = TRUE)

  GAP <- 5L   # gene units of padding between concatenated segments

  # Natural sort: sort data frame by numeric part of a column then by full string
  .nat_sort <- function(df, chr_col, pos_col) {
    num <- suppressWarnings(as.integer(gsub("[^0-9]", "", df[[chr_col]])))
    df[order(num, df[[chr_col]], df[[pos_col]]), ]
  }

  for (pid in par_ids) {
    pr <- pairs[pairs$par_id == pid, , drop = FALSE]

    # Unique segments for this PAR, sorted naturally
    xs <- unique(pr[, c("list_x", "x_min", "x_max")])
    ys <- unique(pr[, c("list_y", "y_min", "y_max")])
    xs <- .nat_sort(xs, "list_x", "x_min")
    ys <- .nat_sort(ys, "list_y", "y_min")
    nx <- nrow(xs); ny <- nrow(ys)

    # Widths / heights and cumulative offsets (segments separated by GAP)
    xw <- xs$x_max - xs$x_min + 1L
    yh <- ys$y_max - ys$y_min + 1L
    xo <- cumsum(c(0L, head(xw, -1L) + GAP))
    yo <- cumsum(c(0L, head(yh, -1L) + GAP))
    tot_x <- tail(xo, 1L) + tail(xw, 1L)
    tot_y <- tail(yo, 1L) + tail(yh, 1L)

    # Assign each anchorpoint to an x-segment index (0 = not in this PAR)
    xi_vec <- integer(nrow(ap))
    for (i in seq_len(nx)) {
      hit <- ap$list_x == xs$list_x[i] &
             ap$coord_x >= xs$x_min[i] &
             ap$coord_x <= xs$x_max[i]
      xi_vec[hit] <- i
    }
    yi_vec <- integer(nrow(ap))
    for (j in seq_len(ny)) {
      hit <- ap$list_y == ys$list_y[j] &
             ap$coord_y >= ys$y_min[j] &
             ap$coord_y <= ys$y_max[j]
      yi_vec[hit] <- j
    }

    keep  <- xi_vec > 0L & yi_vec > 0L
    xi_k  <- xi_vec[keep]
    yi_k  <- yi_vec[keep]
    ap_k  <- ap[keep, , drop = FALSE]
    px    <- xo[xi_k] + (ap_k$coord_x - xs$x_min[xi_k])
    py    <- yo[yi_k] + (ap_k$coord_y - ys$y_min[yi_k])

    # ---- Draw panel ----
    mg <- GAP
    graphics::plot.new()
    graphics::plot.window(xlim = c(-mg, tot_x + mg), ylim = c(-mg, tot_y + mg))

    # White background with border
    graphics::rect(0, 0, tot_x, tot_y, col = "white", border = "black", lwd = 0.8)

    # Green vertical separators between x-segments
    if (nx > 1L) {
      for (i in seq(2L, nx)) {
        xv <- xo[i] - GAP / 2
        graphics::segments(xv, 0, xv, tot_y, col = "green3", lwd = 0.7)
      }
    }
    # Green horizontal separators between y-segments
    if (ny > 1L) {
      for (j in seq(2L, ny)) {
        yh_line <- yo[j] - GAP / 2
        graphics::segments(0, yh_line, tot_x, yh_line, col = "green3", lwd = 0.7)
      }
    }

    # Anchorpoint dots
    if (length(px) > 0L)
      graphics::points(px, py, pch = 20, cex = point_size, col = "black")

    # Title
    graphics::title(main = sprintf("PAR%02d", pid),
                    col.main = "#005bb5", cex.main = 0.9, font.main = 2,
                    line = 0.4)

    # X-axis segment labels (below, rotated 90°)
    graphics::axis(1,
                   at     = xo + xw / 2,
                   labels = paste0(xs$list_x, ":", xs$x_min, "-", xs$x_max),
                   tick   = FALSE, las = 2, cex.axis = label_cex, line = 0.3)

    # Y-axis segment labels (right side, horizontal)
    graphics::axis(4,
                   at     = yo + yh / 2,
                   labels = paste0(ys$list_y, ":", ys$y_min, "-", ys$y_max),
                   tick   = FALSE, las = 2, cex.axis = label_cex, line = 0.3)
  }

  invisible(par_result)
}


# Internal: qualitative colour palette.
# For n <= 16 uses a fixed high-contrast set; for larger n uses golden-angle
# hue spacing in HCL space, which keeps consecutive colours maximally apart.
.dotplot_palette <- function(n) {
  if (n == 0L) return(character(0))
  base_colours <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
  )
  if (n <= length(base_colours)) {
    return(base_colours[seq_len(n)])
  }
  # Golden-angle spacing: each successive hue is ~137.5° away, so no two
  # nearby indices share a similar hue regardless of n
  hues <- (seq(0, n - 1L) * 137.508) %% 360
  grDevices::hcl(h = hues, c = 70, l = 58)
}
