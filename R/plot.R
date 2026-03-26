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
#'   (\code{is_real_anchorpoint == 1}) are shown.
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

  # Keep only anchorpoints that belong to the selected genome pair
  ap <- ap[!is.na(ap$genome_x) &
             ap$genome_x == genome_x & ap$genome_y == genome_y, , drop = FALSE]

  if (real_only) ap <- ap[ap$is_real_anchorpoint == 1, , drop = FALSE]

  if (nrow(ap) == 0L) {
    n_ap   <- nrow(anchorpoints)
    n_real <- sum(anchorpoints$is_real_anchorpoint == 1, na.rm = TRUE)
    if (n_ap == 0L) {
      message("No anchorpoints found in output (anchorpoints.txt is empty).")
    } else if (real_only && n_real == 0L) {
      message("No real anchorpoints found (is_real_anchorpoint == 1 in 0 of ",
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
    pal         <- grDevices::colorRampPalette(
                     c("#c6dbef", "#084594"))(n_lev)   # light→dark blue
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
#' Draws all segments of one multiplicon as horizontal tracks with gene boxes
#' and anchor-point connections as diagonal lines, in the style of Fig. 1 /
#' Fig. 3 of Simillion et al. (2002).
#'
#' The x-axis uses the profile position from \code{list_elements} (already
#' gap-aligned across segments).  Filled boxes = forward-strand genes,
#' open boxes = reverse-strand genes.  Anchor lines are coloured by
#' orientation agreement between the two connected genes.
#'
#' @param output_dir Path to the i-ADHoRe output directory, **or** a list
#'   from \code{\link{read_iadhore_output}}.
#' @param multiplicon_id Integer. The multiplicon to draw.
#' @param real_only Logical. Use only real anchor points (default \code{TRUE}).
#' @param col_same Colour for anchor lines where both genes share the same
#'   strand orientation (default \code{"steelblue"}).
#' @param col_opposite Colour for anchor lines where orientations differ
#'   (default \code{"firebrick"}).
#' @param lwd_ap Line width for anchor lines (default \code{1}).
#'
#' @return Invisibly returns a list with elements \code{segments},
#'   \code{list_elements}, and \code{anchorpoints}.
#' @export
#' @examples
#' \dontrun{
#' plot_multiplicon("output/arath/", multiplicon_id = 1)
#' }
plot_multiplicon <- function(output_dir,
                              multiplicon_id,
                              real_only    = TRUE,
                              col_same     = "steelblue",
                              col_opposite = "firebrick",
                              lwd_ap       = 1) {

  dat  <- .load_iadhore(output_dir)
  segs <- dat$segments[dat$segments$multiplicon == multiplicon_id, ]
  if (nrow(segs) == 0L)
    stop("No segments found for multiplicon ", multiplicon_id)
  segs <- segs[order(segs$order), ]

  # Profile positions for genes in this multiplicon's segments
  le <- dat$list_elements[dat$list_elements$segment %in% segs$id, ]
  le <- merge(le, segs[, c("id", "order", "genome", "list")],
              by.x = "segment", by.y = "id", sort = FALSE)

  ap <- dat$anchorpoints[dat$anchorpoints$multiplicon == multiplicon_id, ]
  if (real_only) ap <- ap[ap$is_real_anchorpoint == 1, ]

  n_segs  <- nrow(segs)
  x_range <- range(le$position, na.rm = TRUE)
  pad     <- max(1L, diff(x_range) * 0.03)

  mult_level <- dat$multiplicons$level[dat$multiplicons$id == multiplicon_id]

  old_par <- graphics::par(mar = c(3.5, 9, 2.5, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(xlim = x_range + c(-pad, pad),
                        ylim = c(0.4, n_segs + 0.6))
  graphics::title(
    main = sprintf("Multiplicon %d  (level %d, %d anchor points)",
                   multiplicon_id, mult_level, nrow(ap)),
    xlab = "Profile position",
    cex.main = 0.9
  )
  graphics::axis(1, cex.axis = 0.8)

  for (i in seq_len(n_segs)) {
    y  <- segs$order[i]
    ge <- le[le$segment == segs$id[i], ]
    if (nrow(ge) == 0L) next

    x_lo <- min(ge$position, na.rm = TRUE)
    x_hi <- max(ge$position, na.rm = TRUE)

    # Baseline
    graphics::segments(x_lo, y, x_hi, y, lwd = 2.5, col = "black")

    # Gene boxes: filled = "+", open = "-"
    bw <- 0.35; bh <- 0.18
    for (j in seq_len(nrow(ge))) {
      xj    <- ge$position[j]
      cfill <- if (!is.na(ge$orientation[j]) && ge$orientation[j] == "-")
                 "white" else "black"
      graphics::rect(xj - bw, y - bh, xj + bw, y + bh,
                     col = cfill, border = "black", lwd = 0.4)
    }

    # Track label
    lbl <- paste0(segs$genome[i], "\n", segs$list[i])
    graphics::mtext(lbl, side = 2, at = y, las = 1, cex = 0.65, line = 0.5)
  }

  # Anchor lines
  for (k in seq_len(nrow(ap))) {
    rx <- le[le$gene == ap$gene_x[k], ]
    ry <- le[le$gene == ap$gene_y[k], ]
    if (nrow(rx) == 0L || nrow(ry) == 0L) next
    ox <- rx$orientation[1]; oy <- ry$orientation[1]
    col_line <- if (!is.na(ox) && !is.na(oy) && ox != oy)
                  col_opposite else col_same
    graphics::segments(rx$position[1], rx$order[1],
                       ry$position[1], ry$order[1],
                       col = col_line, lwd = lwd_ap)
  }

  graphics::legend("topright",
                   legend = c("same orientation", "opposite orientation"),
                   col    = c(col_same, col_opposite),
                   lwd    = 2, bty = "n", cex = 0.7)

  invisible(list(segments     = segs,
                 list_elements = le,
                 anchorpoints  = ap))
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

  # Greedy interval stacking per chromosome
  seg_nr$.stack <- NA_integer_
  for (ch in chroms) {
    idx <- which(seg_nr$list == ch)
    if (length(idx) == 0L) next
    ord <- order(seg_nr$first_coord[idx])
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

    # Segments stacked above baseline
    idx <- which(seg_nr$list == ch & !is.na(seg_nr$.stack))
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


# Internal: simple colour palette
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
  # Fall back to grDevices rainbow for large n
  grDevices::rainbow(n, s = 0.8, v = 0.85)
}
