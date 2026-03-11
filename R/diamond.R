#' Get DIAMOND executable path
#' @keywords internal
get_diamond_bin <- function() {
  path <- Sys.which("diamond")
  if (path == "") {
    stop(
      "DIAMOND not found in PATH.\n",
      "Install it via conda:\n",
      "  conda install -c bioconda diamond"
    )
  }
  path
}


#' Run an all-vs-all DIAMOND blastp search
#'
#' Builds a DIAMOND database from one or more protein FASTA files and runs an
#' all-vs-all blastp search. Output is in BLAST tabular format (outfmt 6),
#' ready for [blast_to_families()] or [parse_blast()].
#'
#' When multiple FASTA files are supplied (e.g. one per species), they are
#' concatenated before building the database so that cross-species hits are
#' included in the output.
#'
#' @param fasta_file Path to a protein FASTA file, or a character vector of
#'   paths (one per species) for multi-species searches
#' @param output_file Path to write the BLAST tabular results (outfmt 6)
#' @param work_dir Directory for the intermediate database files. Defaults to
#'   the same directory as `output_file`
#' @param evalue E-value threshold for reporting hits. Default: 1e-5
#' @param threads Number of CPU threads. Default: 1
#' @param more_sensitive Use `--more-sensitive` mode. Default: FALSE
#' @param verbose Print DIAMOND output. Default: TRUE
#'
#' @return Invisibly returns the path to the output file
#' @export
#'
#' @examples
#' \dontrun{
#' # Single species
#' fasta <- system.file("extdata", "proteome.all_transcripts.vvi.fasta",
#'                      package = "iadhoreR")
#' run_diamond(fasta, output_file = "vvi.blast")
#'
#' # Multiple species (all-vs-all across species)
#' fastas <- c(
#'   system.file("extdata", "proteome.spp1.fasta", package = "iadhoreR"),
#'   system.file("extdata", "proteome.spp2.fasta", package = "iadhoreR")
#' )
#' run_diamond(fastas, output_file = "all_vs_all.blast")
#' }
run_diamond <- function(fasta_file, output_file, work_dir = NULL,
                        evalue = 1e-5, threads = 1,
                        more_sensitive = FALSE, verbose = TRUE) {

  missing_files <- fasta_file[!file.exists(fasta_file)]
  if (length(missing_files) > 0) {
    stop("FASTA file(s) not found: ", paste(missing_files, collapse = ", "))
  }

  fasta_file  <- normalizePath(fasta_file,  mustWork = TRUE)
  output_file <- normalizePath(output_file, mustWork = FALSE)

  if (is.null(work_dir)) work_dir <- dirname(output_file)
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
  work_dir <- normalizePath(work_dir, mustWork = TRUE)

  diamond_bin <- get_diamond_bin()
  db_path     <- file.path(work_dir, "diamond_db")

  stdout_dest <- if (verbose) "" else FALSE
  stderr_dest <- if (verbose) "" else FALSE

  # concatenate multiple FASTAs into a single temp file
  if (length(fasta_file) > 1) {
    message("Combining ", length(fasta_file), " FASTA files...")
    combined_fasta <- file.path(work_dir, "combined.fasta")
    combined_lines <- unlist(lapply(fasta_file, readLines))
    writeLines(combined_lines, combined_fasta)
    query_fasta <- combined_fasta
  } else {
    query_fasta <- fasta_file
  }

  # build database
  message("Building DIAMOND database...")
  db_status <- system2(
    diamond_bin,
    args   = c("makedb",
               "--in", shQuote(query_fasta),
               "--db", shQuote(db_path),
               "--threads", threads),
    stdout = stdout_dest,
    stderr = stderr_dest
  )
  if (db_status != 0) stop("DIAMOND makedb failed with status ", db_status)

  # run all-vs-all blastp
  message("Running all-vs-all DIAMOND blastp...")
  sensitivity <- if (more_sensitive) "--more-sensitive" else "--sensitive"
  blast_status <- system2(
    diamond_bin,
    args   = c("blastp",
               "--db",     shQuote(db_path),
               "--query",  shQuote(query_fasta),
               "--out",    shQuote(output_file),
               "--outfmt", "6",
               "--evalue", evalue,
               "--threads", threads,
               sensitivity),
    stdout = stdout_dest,
    stderr = stderr_dest
  )
  if (blast_status != 0) stop("DIAMOND blastp failed with status ", blast_status)

  message("DIAMOND search complete: ", output_file)
  invisible(output_file)
}
