# Regression tests for the GFF reading / ID-matching logic in prepare_input.R.
# Each block documents the concrete bug it guards against.

# ---- helpers ---------------------------------------------------------------

# Write a vector of lines to a temp file and return its path.
write_tmp <- function(lines, ext) {
  f <- tempfile(fileext = ext)
  writeLines(lines, f)
  f
}

# A small GFF3 modelled on the real failure cases:
#  - the FASTA-matching ID lives in a non-standard key (protein_id),
#  - CDS features repeat the same gene over several rows,
#  - one feature has an unknown strand (".").
sample_gff3 <- function() {
  write_tmp(c(
    "##gff-version 3",
    "chr1\tsrc\tgene\t100\t900\t.\t+\t.\tID=gene:G1;Name=geneone",
    "chr1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;transcript_id=T1",
    "chr1\tsrc\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=transcript:T1;protein_id=P1",
    "chr1\tsrc\tCDS\t400\t600\t.\t+\t1\tID=cds1;Parent=transcript:T1;protein_id=P1",
    "chr1\tsrc\tCDS\t700\t900\t.\t+\t2\tID=cds1;Parent=transcript:T1;protein_id=P1",
    "chr2\tsrc\tgene\t2000\t2500\t.\t-\t.\tID=gene:G2;Name=genetwo",
    "chr2\tsrc\tCDS\t2000\t2500\t.\t-\t0\tID=cds2;Parent=transcript:T2;protein_id=P2",
    "chr2\tsrc\tmRNA\t3000\t3100\t.\t.\t.\tID=transcript:T3;Parent=gene:G3;transcript_id=T3",
    "chr2\tsrc\tCDS\t3000\t3100\t.\t.\t0\tID=cds3;Parent=transcript:T3;protein_id=P3"
  ), ".gff")
}

# ---- extract_gff_attribute -------------------------------------------------

test_that("extract_gff_attribute returns one value per row, NA where absent", {
  # Bug guarded: regmatches() returned only matching rows and was assigned back
  # full-length, recycling values onto the wrong features.
  attrs <- c("ID=A;protein_id=PA", "protein_id=PB;ID=B", "protein_id=PC")
  expect_equal(extract_gff_attribute(attrs, "ID", FALSE), c("A", "B", NA))
  expect_equal(extract_gff_attribute(attrs, "protein_id", FALSE),
               c("PA", "PB", "PC"))
})

test_that("extract_gff_attribute matches whole keys only", {
  # Bug guarded: a request for "ID" must not also match "protein_id".
  attrs <- c("protein_id=PA", "ID=A;protein_id=PA")
  expect_equal(extract_gff_attribute(attrs, "ID", FALSE), c(NA, "A"))
})

test_that("extract_gff_attribute handles GTF quoting", {
  attrs <- c('gene_id "G1"; transcript_id "T1";',
             'transcript_id "T2"; gene_id "G2";',
             'gene_id "G3";')
  expect_equal(extract_gff_attribute(attrs, "transcript_id", TRUE),
               c("T1", "T2", NA))
  expect_equal(extract_gff_attribute(attrs, "gene_id", TRUE),
               c("G1", "G2", "G3"))
})

test_that("discover_gff_keys finds every attribute key", {
  attrs <- c("ID=A;protein_id=PA;Parent=x", "Name=n;ID=B")
  expect_setequal(discover_gff_keys(attrs, FALSE),
                  c("ID", "protein_id", "Parent", "Name"))
})

# ---- parse_gff -------------------------------------------------------------

test_that("parse_gff collapses multi-row features to one entry per gene", {
  # Bug guarded: CDS features wrote the gene id once per CDS segment.
  gff <- sample_gff3()
  out <- tempfile("gl")
  files <- suppressWarnings(
    parse_gff(gff, output_dir = out, genome_name = "g",
              feature_type = "CDS", id_attribute = "protein_id")
  )
  chr1 <- readLines(files[["chr1"]])
  expect_equal(chr1, "P1+")          # three CDS rows -> one entry
})

test_that("parse_gff drops unknown-strand features with a warning", {
  # Bug guarded: a "." strand aborted the whole run via write_gene_list().
  gff <- sample_gff3()
  out <- tempfile("gl")
  expect_warning(
    files <- parse_gff(gff, output_dir = out, genome_name = "g",
                       feature_type = "CDS", id_attribute = "protein_id"),
    "no strand"
  )
  # P3 (strand ".") is excluded; only P2 remains on chr2.
  expect_equal(readLines(files[["chr2"]]), "P2-")
})

test_that("parse_gff orders genes by position", {
  gff <- write_tmp(c(
    "chrX\tsrc\tmRNA\t500\t600\t.\t+\t.\tID=B",
    "chrX\tsrc\tmRNA\t100\t200\t.\t-\t.\tID=A"
  ), ".gff")
  out <- tempfile("gl")
  files <- parse_gff(gff, output_dir = out, genome_name = "g",
                     feature_type = "mRNA", id_attribute = "ID")
  expect_equal(readLines(files[["chrX"]]), c("A-", "B+"))
})

# ---- recommend_parse_gff ---------------------------------------------------

test_that("recommend_parse_gff finds the ID in a non-standard attribute key", {
  # Bug guarded: hardcoded key list missed protein_id, so the real match was
  # never considered and a near-empty wrong attribute was returned silently.
  gff   <- sample_gff3()
  fasta <- write_tmp(c(">P1", "MAAA", ">P2", "MBBB"), ".pep")
  rec <- suppressWarnings(
    capture.output(r <- recommend_parse_gff(gff, fasta))
  )
  expect_equal(r$feature_type, "CDS")
  expect_equal(r$id_attribute, "protein_id")
  expect_equal(r$n_matched, 2)
})

test_that("recommend_parse_gff errors helpfully on a true ID mismatch", {
  # Bug guarded: cryptic "No usable attribute found". Now it names the GFF
  # attribute that corresponds to a FASTA header sub-field.
  gff   <- sample_gff3()
  # primary id (X1/X2) is in no attribute; second token matches protein_id.
  fasta <- write_tmp(c(">X1 P1", "MAAA", ">X2 P2", "MBBB"), ".pep")
  expect_error(
    suppressWarnings(capture.output(recommend_parse_gff(gff, fasta))),
    "protein_id"
  )
})
