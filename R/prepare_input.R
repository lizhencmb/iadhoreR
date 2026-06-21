#' Parse a BLAST tabular result file for use with i-ADHoRe
#'
#' Reads a BLAST outfmt 6 file, removes self-hits, applies an e-value cutoff,
#' and writes a two-column tab-separated file of homologous gene pairs suitable
#' for the `blast_table` setting in i-ADHoRe.
#'
#' Note: i-ADHoRe automatically adds the reciprocal pair (B, A) for every
#' pair (A, B) in the table, so duplicate/reciprocal hits do not need to be
#' kept explicitly.
#'
#' @param blast_file Path to the BLAST tabular output file (outfmt 6)
#' @param output_file Path to write the filtered gene-pair table
#' @param evalue_cutoff Maximum e-value to retain a hit. Default: 1e-5
#' @param min_pident Minimum percent identity (0-100) to retain a hit.
#'   Default: 0 (no filter)
#'
#' @return Invisibly returns a data frame with columns `query` and `subject`
#' @export
#'
#' @examples
#' \dontrun{
#' blast_file <- "vvi.blast"  # output from run_diamond()
#' parse_blast(blast_file, output_file = "blast_table.txt", evalue_cutoff = 1e-5)
#' }
parse_blast <- function(blast_file, output_file, evalue_cutoff = 1e-5,
                        min_pident = 0) {
  col_names <- c("query", "subject", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore")

  hits <- utils::read.table(blast_file, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, col.names = col_names)

  # remove self-hits
  hits <- hits[hits$query != hits$subject, ]

  # apply filters
  hits <- hits[hits$evalue <= evalue_cutoff, ]
  if (min_pident > 0) {
    hits <- hits[hits$pident >= min_pident, ]
  }

  # keep only the gene pair columns, drop duplicates
  pairs <- unique(hits[, c("query", "subject")])

  utils::write.table(pairs, file = output_file, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)

  message("Wrote ", nrow(pairs), " gene pairs to ", output_file)
  invisible(pairs)
}


#' Read gene IDs from a FASTA file (part of header before first space)
#' @keywords internal
read_fasta_ids <- function(fasta_file) {
  lines <- readLines(fasta_file)
  headers <- lines[startsWith(lines, ">")]
  sub("^>([^ \t]+).*", "\\1", headers)
}


#' Extract a named attribute from GFF attribute strings
#'
#' Returns a vector the same length as `attributes`, with NA for rows that do
#' not carry the requested key. The key must match a whole attribute name (so
#' asking for "ID" does not accidentally match "protein_id").
#' @keywords internal
extract_gff_attribute <- function(attributes, key, is_gtf) {
  if (is_gtf) {
    pattern <- paste0("(?:^|;|\\s)", key, '\\s+"([^"]+)"')
  } else {
    pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  }
  m <- regexpr(pattern, attributes, perl = TRUE)
  vals <- rep(NA_character_, length(attributes))
  hit  <- m > 0
  if (any(hit)) {
    matched <- regmatches(attributes, m)
    # strip the key prefix/suffix to leave just the value
    if (is_gtf) {
      matched <- sub(paste0('^(?:[;\\s]*)', key, '\\s+"'), "", matched, perl = TRUE)
      matched <- sub('"$', "", matched)
    } else {
      matched <- sub(paste0("^;?", key, "="), "", matched)
    }
    vals[hit] <- matched
  }
  vals
}


#' Discover all attribute keys present in a vector of GFF attribute strings
#' @keywords internal
discover_gff_keys <- function(attributes, is_gtf) {
  if (is_gtf) {
    pat <- '[A-Za-z_][A-Za-z0-9_]*(?=\\s+")'
  } else {
    pat <- '[A-Za-z_][A-Za-z0-9_]*(?==)'
  }
  unique(unlist(regmatches(attributes, gregexpr(pat, attributes, perl = TRUE))))
}


#' Find which GFF attribute matches a non-primary field of the FASTA headers
#'
#' Used as a fallback when no GFF attribute equals the primary FASTA IDs. Splits
#' each FASTA header into all whitespace tokens and `|`/`:`/`,`-delimited
#' subfields, then reports the (feature_type, attribute) whose values best
#' overlap that expanded set. Returns NULL when nothing overlaps.
#' @keywords internal
diagnose_fasta_mismatch <- function(gff, is_gtf, attr_keys, feat_types,
                                    fasta_file) {
  lines   <- readLines(fasta_file)
  headers <- sub("^>", "", lines[startsWith(lines, ">")])
  tokens  <- unlist(strsplit(headers, "[ \t]+"))
  subflds <- unlist(strsplit(tokens, "[|:,]"))
  token_set <- unique(c(tokens, subflds))
  token_set <- token_set[nzchar(token_set)]

  best <- NULL
  for (ft in feat_types) {
    sub_gff <- gff$attributes[gff$type == ft]
    if (length(sub_gff) == 0) next
    for (ak in attr_keys) {
      ids <- extract_gff_attribute(sub_gff, ak, is_gtf)
      uid <- unique(ids[!is.na(ids)])
      if (length(uid) == 0) next
      n   <- length(intersect(uid, token_set))
      if (n > 0 && (is.null(best) || n > best$n)) {
        best <- list(feature_type = ft, id_attribute = ak, n = n)
      }
    }
  }
  best
}


#' Recommend parse_gff() parameters by matching a GFF file against a FASTA file
#'
#' Inspects all combinations of feature type and attribute key in a GFF file,
#' finds the combination whose IDs best overlap with the sequence IDs in a
#' protein FASTA file, and returns the recommended parameters for [parse_gff()].
#'
#' If the number of matched genes differs between the two files, the
#' intersection is used and a warning is issued.
#'
#' @param gff_file Path to the GFF3 or GTF file
#' @param fasta_file Path to the protein FASTA file
#'
#' @return A named list with elements `feature_type`, `id_attribute`,
#'   `chromosomes` (all sequences containing matched genes), `n_fasta`
#'   (IDs in FASTA), `n_gff` (IDs in best GFF combination), and `n_matched`
#'   (overlap). Pass the list elements directly to [parse_gff()].
#' @importFrom utils head
#' @export
#'
#' @examples
#' \dontrun{
#' gff   <- system.file("extdata",
#'            "annotation.selected_transcript.all_features.vvi.gff3",
#'            package = "iadhoreR")
#' fasta <- system.file("extdata", "proteome.all_transcripts.vvi.fasta",
#'            package = "iadhoreR")
#' rec   <- recommend_parse_gff(gff, fasta)
#'
#' gene_lists <- parse_gff(gff,
#'   output_dir    = "vvin_lists",
#'   genome_name   = "vvin",
#'   feature_type  = rec$feature_type,
#'   id_attribute  = rec$id_attribute,
#'   chromosomes   = rec$chromosomes)
#' }
recommend_parse_gff <- function(gff_file, fasta_file) {

  fasta_ids <- read_fasta_ids(fasta_file)
  n_fasta   <- length(fasta_ids)
  cat("FASTA: ", n_fasta, " sequences\n", sep = "")
  cat("  e.g.: ", paste(head(fasta_ids, 3), collapse = ", "), "\n\n")

  # read full GFF once
  gff <- utils::read.table(gff_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "",
                           comment.char = "#",
                           col.names = c("seqid", "source", "type", "start",
                                         "end", "score", "strand", "phase",
                                         "attributes"))

  is_gtf <- grepl('"', gff$attributes[1])
  cat("GFF format: ", if (is_gtf) "GTF" else "GFF3", "\n\n", sep = "")

  # candidate feature types (by frequency) and attribute keys. Discover every
  # key actually present in the file rather than assuming a fixed set: the ID
  # that matches the FASTA can live in non-standard keys such as protein_id,
  # Source_ID, proteinId or transcriptId depending on the annotation source.
  feat_types  <- names(sort(table(gff$type), decreasing = TRUE))
  preferred   <- c("ID", "Name", "transcript_id", "gene_id", "protein_id",
                   "Source_ID", "proteinId", "transcriptId", "Parent")
  discovered  <- discover_gff_keys(gff$attributes, is_gtf)
  attr_keys   <- unique(c(intersect(preferred, discovered),
                          setdiff(discovered, preferred)))

  # score every combination by UNIQUE matched IDs
  # prefer: (1) most unique matches, (2) fewest duplicates per unique ID
  fasta_set <- as.character(fasta_ids)
  results   <- vector("list", length(feat_types) * length(attr_keys))
  idx       <- 0L
  for (ft in feat_types) {
    sub_gff <- gff[gff$type == ft, ]
    if (nrow(sub_gff) == 0) next
    for (ak in attr_keys) {
      ids <- extract_gff_attribute(sub_gff$attributes, ak, is_gtf)
      if (all(is.na(ids))) next
      unique_ids <- unique(ids[!is.na(ids)])
      n_unique_matched <- length(intersect(unique_ids, fasta_set))
      if (n_unique_matched == 0) next
      # duplication ratio: avg features per unique ID (1.0 = no duplicates)
      dup_ratio <- length(ids[!is.na(ids)]) / length(unique_ids)
      idx <- idx + 1L
      results[[idx]] <- list(
        feature_type     = ft,
        id_attribute     = ak,
        n_gff_unique     = length(unique_ids),
        n_matched        = n_unique_matched,
        dup_ratio        = dup_ratio,
        gff_sub          = sub_gff,
        ids              = ids
      )
    }
  }
  results <- results[seq_len(idx)]

  # decide whether the primary IDs were found in any attribute at all
  best_n <- if (length(results) == 0) 0L else max(sapply(results, `[[`, "n_matched"))

  if (best_n < 0.5 * n_fasta) {
    # The FASTA primary IDs are not (mostly) stored as a GFF attribute value.
    # Look for where they hide so the message is actionable instead of cryptic.
    diag <- diagnose_fasta_mismatch(gff, is_gtf, attr_keys, feat_types, fasta_file)
    msg <- paste0(
      "The primary FASTA IDs do not match any GFF attribute value",
      if (best_n > 0) paste0(" well (best: ", best_n, "/", n_fasta, ").") else ".",
      "\n  FASTA ID example : ", paste(head(fasta_ids, 2), collapse = ", "))
    if (!is.null(diag)) {
      msg <- paste0(msg,
        "\n  However GFF attribute '", diag$id_attribute, "' (feature '",
        diag$feature_type, "') matches ", diag$n, " of the FASTA header",
        " sub-fields.\n  Your FASTA primary IDs differ from the IDs stored in",
        " the GFF, so gene lists built from this GFF will NOT match a BLAST",
        " table built from this FASTA.\n  Reconcile the IDs first (rename the",
        " FASTA headers, or rebuild the BLAST table using the same IDs as the",
        " GFF), then re-run.")
    } else {
      msg <- paste0(msg,
        "\n  No GFF attribute matched any FASTA header field. Check that the",
        " GFF and FASTA correspond to the same annotation, and run",
        " inspect_gff() to see the available attribute keys.")
    }
    stop(msg, call. = FALSE)
  }

  # pick best: most unique matches; break ties by lowest duplication ratio
  n_matched_vec <- sapply(results, `[[`, "n_matched")
  dup_vec       <- sapply(results, `[[`, "dup_ratio")
  best_score    <- max(n_matched_vec)
  candidates    <- which(n_matched_vec == best_score)
  best          <- results[[candidates[which.min(dup_vec[candidates])]]]

  cat("Best match: feature_type = \"", best$feature_type, "\", ",
      "id_attribute = \"", best$id_attribute, "\"\n", sep = "")
  cat("  GFF unique IDs : ", best$n_gff_unique, "\n", sep = "")
  cat("  FASTA          : ", n_fasta,            "\n", sep = "")
  cat("  Matched        : ", best$n_matched,     "\n", sep = "")

  # warn if counts differ
  unique_gff_ids <- unique(best$ids[!is.na(best$ids)])
  only_gff       <- setdiff(unique_gff_ids, fasta_set)
  only_fasta     <- setdiff(fasta_set, unique_gff_ids)

  if (length(only_gff) > 0) {
    warning(length(only_gff), " gene(s) found in GFF but NOT in FASTA -- ",
            "these will be excluded.\n",
            "  e.g.: ", paste(head(only_gff, 3), collapse = ", "),
            call. = FALSE)
  }
  if (length(only_fasta) > 0) {
    warning(length(only_fasta), " sequence(s) found in FASTA but NOT in GFF -- ",
            "these will be excluded.\n",
            "  e.g.: ", paste(head(only_fasta, 3), collapse = ", "),
            call. = FALSE)
  }

  # chromosomes that contain at least one matched gene
  matched_mask <- !is.na(best$ids) & best$ids %in% fasta_set
  chromosomes  <- sort(unique(best$gff_sub$seqid[matched_mask]))
  cat("  Chromosomes with matched genes (", length(chromosomes), "):\n", sep = "")
  cat("  ", paste(chromosomes, collapse = ", "), "\n\n", sep = "")

  cat("Recommended call:\n")
  chr_str <- paste0('c("', paste(chromosomes, collapse = '", "'), '")')
  cat('  parse_gff(gff_file,\n')
  cat('    output_dir   = "...",\n')
  cat('    genome_name  = "...",\n')
  cat('    feature_type = "', best$feature_type, '",\n', sep = "")
  cat('    id_attribute = "', best$id_attribute,  '",\n', sep = "")
  cat('    chromosomes  = ', chr_str, ')\n', sep = "")

  invisible(list(
    feature_type = best$feature_type,
    id_attribute = best$id_attribute,
    chromosomes  = chromosomes,
    n_fasta      = n_fasta,
    n_gff        = best$n_gff_unique,
    n_matched    = best$n_matched
  ))
}


#' Inspect a GFF file to guide parse_gff() parameter choices
#'
#' Summarises the feature types and attribute keys present in a GFF file.
#' Use this on any new GFF file to determine the correct `feature_type` and
#' `id_attribute` values to pass to [parse_gff()].
#'
#' @param gff_file Path to the GFF/GTF file
#' @param n_features Number of features to sample per type when showing
#'   example IDs. Default: 3
#'
#' @return Invisibly returns a list with elements `feature_counts` (a named
#'   integer vector) and `attributes` (a named list of attribute key counts
#'   per feature type)
#' @export
#'
#' @examples
#' \dontrun{
#' gff_file <- system.file("extdata", "annotation.selected_transcript.all_features.vvi.gff3",
#'                         package = "iadhoreR")
#' inspect_gff(gff_file)
#' }
inspect_gff <- function(gff_file, n_features = 3) {
  gff <- utils::read.table(gff_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "",
                           comment.char = "#",
                           col.names = c("seqid", "source", "type", "start",
                                         "end", "score", "strand", "phase",
                                         "attributes"))

  # feature type counts
  feat_counts <- sort(table(gff$type), decreasing = TRUE)
  cat("Feature types:\n")
  print(feat_counts)

  # detect format: GTF uses key "value"; GFF3 uses key=value
  sample_attr <- gff$attributes[!is.na(gff$attributes)][1]
  is_gtf <- grepl('"', sample_attr)

  # parse attribute keys for each feature type
  attr_info <- lapply(names(feat_counts), function(ftype) {
    rows <- gff[gff$type == ftype, ]
    if (nrow(rows) == 0) return(NULL)

    if (is_gtf) {
      keys <- regmatches(rows$attributes,
                         gregexpr('[A-Za-z_][A-Za-z0-9_]*(?=\\s+")',
                                  rows$attributes, perl = TRUE))
    } else {
      keys <- regmatches(rows$attributes,
                         gregexpr('[A-Za-z_][A-Za-z0-9_]*(?==)',
                                  rows$attributes, perl = TRUE))
    }
    key_counts <- sort(table(unlist(keys)), decreasing = TRUE)

    # show example IDs for the most likely ID attributes
    id_candidates <- intersect(c("ID", "Name", "gene_id", "transcript_id",
                                 "Parent"), names(key_counts))
    examples <- lapply(id_candidates, function(k) {
      if (is_gtf) {
        pattern <- paste0(k, '\\s+"([^"]+)"')
      } else {
        pattern <- paste0("(?:^|;)", k, "=([^;]+)")
      }
      vals <- regmatches(rows$attributes[seq_len(min(n_features, nrow(rows)))],
                         regexpr(pattern, rows$attributes[seq_len(
                           min(n_features, nrow(rows)))]))
      sub(paste0(".*", k, '[="]+"?'), "", vals)
    })
    names(examples) <- id_candidates
    list(key_counts = key_counts, examples = examples)
  })
  names(attr_info) <- names(feat_counts)

  cat("\nFormat detected:", if (is_gtf) "GTF" else "GFF3", "\n")
  cat("\nAttribute keys per feature type (showing ID candidates with examples):\n")
  for (ftype in names(attr_info)) {
    info <- attr_info[[ftype]]
    if (is.null(info)) next
    cat("\n  [", ftype, "]\n", sep = "")
    cat("  All keys:", paste(names(info$key_counts), collapse = ", "), "\n")
    for (k in names(info$examples)) {
      ex <- info$examples[[k]]
      if (length(ex) > 0)
        cat("  ", k, " e.g.: ", paste(ex, collapse = ", "), "\n", sep = "")
    }
  }

  cat("\nChromosomes/sequences:", paste(sort(unique(gff$seqid)), collapse = ", "), "\n")

  invisible(list(feature_counts = feat_counts, attributes = attr_info,
                 format = if (is_gtf) "GTF" else "GFF3"))
}


#' Parse a GFF/GTF file to create i-ADHoRe gene list files
#'
#' Extracts features of a given type from a GFF3 or GTF file, sorts genes by
#' chromosomal position, and writes one gene list file per chromosome. Each
#' file contains one gene identifier per line with its strand orientation
#' appended (e.g. `AT1G01010+`), as required by i-ADHoRe.
#'
#' Use [inspect_gff()] first on any new file to determine the correct values
#' for `feature_type` and `id_attribute`.
#'
#' The returned named vector can be passed directly as one genome entry in the
#' `genomes` argument of [write_iadhore_config()].
#'
#' @param gff_file Path to the GFF3 or GTF file
#' @param output_dir Directory where gene list files will be written (created
#'   if it does not exist)
#' @param genome_name Genome identifier used as the file name prefix and
#'   returned as the name of the result vector. Default: "genome"
#' @param feature_type GFF feature type to extract gene identifiers from.
#'   Default: "mRNA"
#' @param id_attribute Attribute key to use as the gene identifier. For GFF3
#'   files this is typically "ID" or "Name"; for GTF files "transcript_id" or
#'   "gene_id". Default: "ID"
#' @param chromosomes Character vector of chromosomes to include. If NULL
#'   (default), all sequences present in the GFF are included.
#'
#' @return A named character vector mapping chromosome identifiers to gene list
#'   file paths, suitable for use in [write_iadhore_config()]
#' @export
#'
#' @examples
#' \dontrun{
#' gff_file <- system.file("extdata", "annotation.selected_transcript.all_features.vvi.gff3",
#'                         package = "iadhoreR")
#'
#' # inspect first
#' inspect_gff(gff_file)
#'
#' # then parse
#' gene_lists <- parse_gff(gff_file, output_dir = "vvin_lists",
#'                         genome_name = "vvin", feature_type = "mRNA",
#'                         id_attribute = "ID")
#' }
parse_gff <- function(gff_file, output_dir, genome_name = "genome",
                      feature_type = "mRNA", id_attribute = "ID",
                      chromosomes = NULL) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  gff <- utils::read.table(gff_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "",
                           comment.char = "#",
                           col.names = c("seqid", "source", "type", "start",
                                         "end", "score", "strand", "phase",
                                         "attributes"))

  gff <- gff[gff$type == feature_type, ]
  if (nrow(gff) == 0) {
    stop("No features of type '", feature_type, "' found in ", gff_file,
         "\nRun inspect_gff() to see available feature types.")
  }

  # detect GTF vs GFF3 and extract the requested attribute
  sample_attr <- gff$attributes[1]
  is_gtf <- grepl('"', sample_attr)

  gff$gene_id <- extract_gff_attribute(gff$attributes, id_attribute, is_gtf)

  missing_id <- is.na(gff$gene_id) | gff$gene_id == ""
  if (any(missing_id)) {
    stop(sum(missing_id), " features have no '", id_attribute, "' attribute.",
         "\nRun inspect_gff() to see available attribute keys.")
  }

  # i-ADHoRe gene lists require a known orientation; drop features whose strand
  # is unknown (".") rather than aborting the whole file.
  bad_strand <- !(gff$strand %in% c("+", "-"))
  if (any(bad_strand)) {
    warning(sum(bad_strand), " feature(s) have no strand ('.') and were ",
            "excluded.", call. = FALSE)
    gff <- gff[!bad_strand, ]
  }

  if (!is.null(chromosomes)) {
    gff <- gff[gff$seqid %in% chromosomes, ]
    not_found <- setdiff(chromosomes, unique(gff$seqid))
    if (length(not_found) > 0) {
      warning("Chromosomes not found in GFF: ", paste(not_found, collapse = ", "))
    }
  }

  chr_list  <- split(gff, gff$seqid)
  file_paths <- vapply(names(chr_list), function(chr) {
    chr_genes <- chr_list[[chr]][order(chr_list[[chr]]$start), ]
    # multi-segment feature types (e.g. CDS/exon) repeat the same gene id on
    # consecutive rows; keep one entry per gene at its first (lowest) position.
    chr_genes <- chr_genes[!duplicated(chr_genes$gene_id), ]
    out_file  <- file.path(output_dir, paste0(genome_name, "_", chr, ".lst"))
    write_gene_list(chr_genes$gene_id, chr_genes$strand, out_file)
    out_file
  }, character(1))

  names(file_paths) <- names(chr_list)
  message("Wrote ", length(file_paths), " gene list files to ", output_dir)
  invisible(file_paths)
}
