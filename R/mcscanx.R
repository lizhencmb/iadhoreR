#' Convert i-ADHoRe results to MCScanX collinearity format
#'
#' Reads i-ADHoRe output tables and writes a `.collinearity` file in the
#' format produced by MCScanX. Each i-ADHoRe basecluster becomes one
#' alignment block.
#'
#' **Score** is set to `number_of_anchorpoints` (a proxy; MCScanX uses a
#' Smith-Waterman score that i-ADHoRe does not compute).
#'
#' **E-value** in the block header is `baseclusters.random_probability`.
#' Per-pair e-values on each gene-pair line are taken from the original BLAST
#' file when `blast_file` is supplied, otherwise the basecluster e-value is
#' used for every pair.
#'
#' @param output_dir Path to the i-ADHoRe output directory
#' @param file Path to write the `.collinearity` file
#' @param blast_file Optional path to the original BLAST tabular file
#'   (outfmt 6). When provided, the best e-value for each gene pair is looked
#'   up and written on the gene-pair lines. Default: NULL
#' @param include_redundant If FALSE (default), baseclusters belonging to
#'   redundant multiplicons are skipped, giving a non-overlapping block set
#'   comparable to MCScanX output.
#'
#' @return Invisibly returns the path to the written file.
#' @export
#'
#' @examples
#' \dontrun{
#' results_dir <- "iadhore_output/"
#' write_mcscanx_collinearity(results_dir, file = "iadhore.collinearity")
#'
#' # With per-pair e-values from the original BLAST file
#' write_mcscanx_collinearity(results_dir, file = "iadhore.collinearity",
#'                            blast_file = "all_vs_all.blast")
#' }
write_mcscanx_collinearity <- function(output_dir, file,
                                       blast_file = NULL,
                                       include_redundant = FALSE) {

  multiplicons <- read_multiplicons(output_dir)
  baseclusters <- read_baseclusters(output_dir)
  anchorpoints <- read_anchorpoints(output_dir)

  # optionally build a per-pair e-value lookup from the BLAST file
  pair_evalue <- NULL
  if (!is.null(blast_file)) {
    col_names <- c("query", "subject", "pident", "length", "mismatch",
                   "gapopen", "qstart", "qend", "sstart", "send",
                   "evalue", "bitscore")
    blast <- utils::read.table(blast_file, header = FALSE, sep = "\t",
                               stringsAsFactors = FALSE, col.names = col_names)
    blast <- blast[order(blast$evalue), ]
    blast <- blast[!duplicated(blast[, c("query", "subject")]), ]
    pair_evalue <- setNames(blast$evalue,
                            paste(blast$query, blast$subject, sep = "\t"))
    message("Loaded ", nrow(blast), " gene pairs from BLAST file.")
  }

  # optionally exclude baseclusters from redundant multiplicons
  if (!include_redundant) {
    non_redundant_mp <- multiplicons$id[multiplicons$is_redundant == 0]
    baseclusters <- baseclusters[baseclusters$multiplicon %in% non_redundant_mp, ]
  }

  lines <- character(0)
  align_id <- 0L

  for (i in seq_len(nrow(baseclusters))) {
    bc <- baseclusters[i, ]
    mp <- multiplicons[multiplicons$id == bc$multiplicon, ]

    chr_x  <- paste0(mp$genome_x, "_", mp$list_x)
    chr_y  <- paste0(mp$genome_y, "_", mp$list_y)
    orient <- if (bc$orientation == "+") "plus" else "minus"

    aps <- anchorpoints[anchorpoints$basecluster == bc$id, ]
    aps <- aps[order(aps$coord_x), ]

    header <- sprintf(
      "## Alignment %d: score=%d e_value=%.2e N=%d %s&%s %s",
      align_id,
      bc$number_of_anchorpoints,
      bc$random_probability,
      bc$number_of_anchorpoints,
      chr_x, chr_y, orient
    )
    lines <- c(lines, header)

    for (j in seq_len(nrow(aps))) {
      ap <- aps[j, ]

      if (!is.null(pair_evalue)) {
        key_fwd <- paste(ap$gene_x, ap$gene_y, sep = "\t")
        key_rev <- paste(ap$gene_y, ap$gene_x, sep = "\t")
        ev <- if (!is.na(pair_evalue[key_fwd])) pair_evalue[[key_fwd]] else
              if (!is.na(pair_evalue[key_rev])) pair_evalue[[key_rev]] else
              bc$random_probability
      } else {
        ev <- bc$random_probability
      }

      lines <- c(lines,
        sprintf("  %d-%3d:\t%s\t%s\t%.2e",
                align_id, j - 1L, ap$gene_x, ap$gene_y, ev))
    }

    align_id <- align_id + 1L
  }

  writeLines(lines, file)
  message("Wrote ", align_id, " alignment blocks to ", file)
  invisible(file)
}


#' Convert i-ADHoRe gene table to MCScanX GFF format
#'
#' Writes a four-column tab-separated file (`chromosome`, `gene`, `start`,
#' `stop`) from i-ADHoRe's `genes.txt`, as required by MCScanX and its
#' downstream visualisation tools.
#'
#' Chromosome names are constructed as `genome_list`
#' (e.g. `ath_Chr1`), matching the names in the `.collinearity` file produced
#' by [write_mcscanx_collinearity()].
#'
#' **Note on coordinates:** i-ADHoRe's `genes.txt` stores gene-order
#' coordinates (1, 2, 3, ...), not base-pair positions. These are sufficient
#' for collinearity analysis and dot-plot visualisation. If you need true bp
#' coordinates for other MCScanX tools, provide the original GFF/GTF file via
#' the `gff_file` argument.
#'
#' @param output_dir Path to the i-ADHoRe output directory
#' @param file Path to write the MCScanX GFF file
#' @param gff_file Optional path to the original GFF3/GTF annotation file.
#'   When provided, actual start/stop base-pair coordinates replace the
#'   gene-order coordinates. `feature_type` and `id_attribute` must also be
#'   supplied. Default: NULL
#' @param feature_type GFF feature type used to extract coordinates (e.g.
#'   `"mRNA"`). Only used when `gff_file` is provided.
#' @param id_attribute GFF attribute key for gene identifiers (e.g. `"ID"`).
#'   Only used when `gff_file` is provided.
#'
#' @return Invisibly returns the path to the written file.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using gene-order coordinates (default)
#' write_mcscanx_gff("iadhore_output/", file = "iadhore.gff")
#'
#' # Using real bp coordinates from the original GFF
#' write_mcscanx_gff("iadhore_output/", file = "iadhore.gff",
#'                   gff_file      = "annotation.gff3",
#'                   feature_type  = "mRNA",
#'                   id_attribute  = "ID")
#' }
write_mcscanx_gff <- function(output_dir, file,
                               gff_file = NULL,
                               feature_type = "mRNA",
                               id_attribute = "ID") {

  genes <- read_genes(output_dir)
  genes$chr <- paste0(genes$genome, "_", genes$list)

  if (!is.null(gff_file)) {
    # read bp coordinates from original GFF
    raw <- utils::read.table(gff_file, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "",
                             comment.char = "#",
                             col.names = c("seqid", "source", "type", "start",
                                           "end", "score", "strand", "phase",
                                           "attributes"))
    raw <- raw[raw$type == feature_type, ]

    is_gtf <- grepl('"', raw$attributes[1])
    if (is_gtf) {
      pattern <- paste0(id_attribute, '\\s+"([^"]+)"')
      m <- regexpr(pattern, raw$attributes)
      raw$gene_id <- regmatches(raw$attributes, m)
      raw$gene_id <- sub(paste0(id_attribute, '\\s+"'), "", raw$gene_id)
      raw$gene_id <- sub('"$', "", raw$gene_id)
    } else {
      pattern <- paste0("(?:^|;)", id_attribute, "=([^;]+)")
      m <- regexpr(pattern, raw$attributes, perl = TRUE)
      raw$gene_id <- regmatches(raw$attributes, m)
      raw$gene_id <- sub(paste0("^;?", id_attribute, "="), "", raw$gene_id)
    }

    coord_map <- setNames(
      data.frame(start = raw$start, stop = raw$end, stringsAsFactors = FALSE),
      c("start", "stop")
    )
    rownames(coord_map) <- raw$gene_id

    matched <- genes$id %in% rownames(coord_map)
    if (!all(matched)) {
      warning(sum(!matched), " gene(s) in i-ADHoRe output not found in GFF; ",
              "gene-order coordinates used for those genes.")
    }

    start <- ifelse(matched, coord_map[genes$id, "start"], genes$coordinate)
    stop  <- ifelse(matched, coord_map[genes$id, "stop"],  genes$coordinate)
  } else {
    start <- genes$coordinate
    stop  <- genes$coordinate
  }

  out <- data.frame(
    chr   = genes$chr,
    gene  = genes$id,
    start = start,
    stop  = stop,
    stringsAsFactors = FALSE
  )

  utils::write.table(out, file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  message("Wrote ", nrow(out), " genes to ", file,
          if (is.null(gff_file)) " (gene-order coordinates)" else
          " (bp coordinates from GFF)")
  invisible(file)
}
