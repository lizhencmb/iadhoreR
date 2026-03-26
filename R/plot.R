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
#'
#' @return Invisibly returns the data frame of points that were plotted.
#' @export
#' @examples
#' \dontrun{
#' plot_dotplot("output/vvi/")
#' plot_dotplot("output/two_sp/", genome_x = "sp1", genome_y = "sp2",
#'              chr_order = "alpha")
#' }
plot_dotplot <- function(output_dir,
                         genome_x   = NULL,
                         genome_y   = NULL,
                         chr_order  = c("size", "natural", "alpha"),
                         real_only  = TRUE,
                         colour_by  = c("level", "multiplicon", "basecluster"),
                         point_size = 0.8) {

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

  # Dots
  graphics::points(ap$gx, ap$gy, col = ap$.col, pch = 20, cex = point_size)

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
