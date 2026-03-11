#' Get path to an MCL tool binary
#' @param tool One of "mcl", "mcxload", "mcxdump"
#' @keywords internal
get_mcl_bin <- function(tool) {
  path <- Sys.which(tool)
  if (path == "") {
    stop(
      "'", tool, "' not found in PATH.\n",
      "Install the MCL suite via conda:\n",
      "  conda install -c bioconda mcl"
    )
  }
  path
}


#' Cluster BLAST/DIAMOND results into gene families using MCL
#'
#' Replicates the `mclblastline` pipeline using the individual MCL tools:
#' `mcxload` to build the graph, `mcl` to cluster, `mcxdump` to export.
#' Writes a two-column gene-to-family table ready for i-ADHoRe
#' (`table_type=family`).
#'
#' @param blast_file Path to BLAST/DIAMOND tabular output (outfmt 6)
#' @param output_file Path to write the gene-to-family table
#' @param gene_list_files Character vector of gene list file paths (as returned
#'   by [parse_gff()]). When provided, any gene present in the gene lists but
#'   absent from the BLAST results is added as a singleton family. Required
#'   for `table_type=family` in i-ADHoRe. Default: NULL
#' @param work_dir Directory for intermediate MCL files. Defaults to the same
#'   directory as `output_file`
#' @param inflation MCL inflation parameter (cluster granularity). Higher =
#'   smaller, tighter clusters. Suggested: 1.5-5.0. Default: 3.0
#' @param evalue_cutoff Maximum e-value to include a hit. Default: 1e-5
#' @param threads Number of threads for MCL. Default: 1
#' @param verbose Print MCL tool output. Default: FALSE
#'
#' @return Invisibly returns a data frame with columns `gene` and `family`
#' @export
#'
#' @examples
#' \dontrun{
#' # Single species
#' blast_file <- "vvi.blast"  # output from run_diamond()
#' families <- blast_to_families(blast_file, output_file = "families.txt",
#'                               gene_list_files = unname(vvi_gene_lists))
#'
#' # Multiple species: combine gene list files from all species
#' all_gene_lists <- c(unname(ath_gene_lists), unname(vvi_gene_lists))
#' families <- blast_to_families("all_vs_all.blast",
#'                               output_file  = "families.txt",
#'                               gene_list_files = all_gene_lists)
#'
#' write_iadhore_config(
#'   genomes     = list(ath = ath_gene_lists, vvi = vvi_gene_lists),
#'   blast_table = "families.txt",
#'   table_type  = "family",
#'   output_path = "output/",
#'   file        = "config.ini"
#' )
#' }
blast_to_families <- function(blast_file, output_file, gene_list_files = NULL,
                              work_dir = NULL, inflation = 3.0,
                              evalue_cutoff = 1e-5, threads = 1,
                              verbose = FALSE) {

  if (!file.exists(blast_file)) {
    stop("BLAST file not found: ", blast_file)
  }

  blast_file  <- normalizePath(blast_file,  mustWork = TRUE)
  output_file <- normalizePath(output_file, mustWork = FALSE)

  if (is.null(work_dir)) work_dir <- dirname(output_file)
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
  work_dir <- normalizePath(work_dir, mustWork = TRUE)

  stdout_dest <- if (verbose) "" else FALSE
  stderr_dest <- if (verbose) "" else FALSE

  # Step 1: parse BLAST outfmt 6 -> ABC format (gene1, gene2, -log10(evalue))
  message("Preparing ABC input from BLAST results...")
  col_names <- c("query", "subject", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore")
  hits <- utils::read.table(blast_file, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, col.names = col_names)

  hits <- hits[hits$query != hits$subject, ]
  hits <- hits[hits$evalue <= evalue_cutoff, ]
  hits$evalue[hits$evalue == 0] <- .Machine$double.xmin  # avoid log(0)

  # keep best (lowest) e-value per ordered pair
  hits <- hits[order(hits$evalue), ]
  hits <- hits[!duplicated(hits[, c("query", "subject")]), ]

  abc <- data.frame(
    gene1 = hits$query,
    gene2 = hits$subject,
    score = -log10(hits$evalue),
    stringsAsFactors = FALSE
  )
  abc_file <- file.path(work_dir, "blast.abc")
  utils::write.table(abc, abc_file, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  message("  ", nrow(abc), " pairs written")

  # Step 2: mcxload - convert ABC to MCL matrix, mirror edges
  mci_file <- file.path(work_dir, "blast.mci")
  tab_file <- file.path(work_dir, "blast.tab")
  message("Running mcxload...")
  load_status <- system2(
    get_mcl_bin("mcxload"),
    args   = c("-abc",          shQuote(abc_file),
               "--stream-mirror",
               "-o",            shQuote(mci_file),
               "-write-tab",    shQuote(tab_file)),
    stdout = stdout_dest,
    stderr = stderr_dest
  )
  if (load_status != 0) stop("mcxload failed with status ", load_status)

  # Step 3: mcl - cluster
  clusters_mci <- file.path(work_dir, "clusters.mci")
  message("Running MCL (inflation = ", inflation, ")...")
  mcl_status <- system2(
    get_mcl_bin("mcl"),
    args   = c(shQuote(mci_file),
               "-I",  inflation,
               "-te", threads,
               "-o",  shQuote(clusters_mci)),
    stdout = stdout_dest,
    stderr = stderr_dest
  )
  if (mcl_status != 0) stop("mcl failed with status ", mcl_status)

  # Step 4: mcxdump - export clusters as tab-separated gene lists
  dump_file <- file.path(work_dir, "clusters.txt")
  message("Running mcxdump...")
  dump_status <- system2(
    get_mcl_bin("mcxdump"),
    args   = c("-icl",   shQuote(clusters_mci),
               "-tabr",  shQuote(tab_file),
               "-o",     shQuote(dump_file)),
    stdout = stdout_dest,
    stderr = stderr_dest
  )
  if (dump_status != 0) stop("mcxdump failed with status ", dump_status)

  # Step 5: parse clusters (one cluster per line, genes tab-separated)
  # and write gene-to-family table for i-ADHoRe
  message("Writing gene-to-family table...")
  cluster_lines <- readLines(dump_file)
  families <- do.call(rbind, lapply(seq_along(cluster_lines), function(i) {
    genes <- strsplit(cluster_lines[i], "\t")[[1]]
    data.frame(
      gene   = genes,
      family = paste0("FAM", formatC(i, width = 6, flag = "0")),
      stringsAsFactors = FALSE
    )
  }))

  # Step 6: add singletons for genes not present in any BLAST hit
  if (!is.null(gene_list_files)) {
    all_genes <- unlist(lapply(gene_list_files, function(f) {
      sub("[+-]$", "", readLines(f))
    }))
    missing <- setdiff(all_genes, families$gene)
    if (length(missing) > 0) {
      next_id <- length(cluster_lines) + seq_along(missing)
      singletons <- data.frame(
        gene   = missing,
        family = paste0("FAM", formatC(next_id, width = 6, flag = "0")),
        stringsAsFactors = FALSE
      )
      families <- rbind(families, singletons)
      message("  Added ", length(missing), " singleton families")
    }
  }

  utils::write.table(families, output_file, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("Done: ", nrow(families), " genes in ",
          length(unique(families$family)), " families -> ", output_file)
  invisible(families)
}
