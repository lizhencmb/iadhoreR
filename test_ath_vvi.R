# =============================================================================
# iadhoreR — Arabidopsis thaliana + Vitis vinifera two-species test
# =============================================================================
# Run from the repo root with the iadhoreR conda env active:
#   conda activate iadhoreR
#   Rscript test_ath_vvi.R
# or from inside an R session:
#   devtools::load_all(); source("test_ath_vvi.R")
# =============================================================================

devtools::load_all(".")

# ── Paths ──────────────────────────────────────────────────────────────────────
ext  <- system.file("extdata", package = "iadhoreR")
work <- "test_output/ath_vvi"
dir.create(work, recursive = TRUE, showWarnings = FALSE)

ath_gff   <- file.path(ext, "annotation.selected_transcript.all_features.ath.gff3")
ath_fasta <- file.path(ext, "proteome.selected_transcript.ath.fasta")
vvi_gff   <- file.path(ext, "annotation.selected_transcript.all_features.vvi.gff3")
vvi_fasta <- file.path(ext, "proteome.selected_transcript.vvi.fasta")

cat("\n=== iadhoreR Arabidopsis + Vitis two-species test ===\n")
cat("Working directory:", work, "\n\n")

# ── Step 0: check tools ────────────────────────────────────────────────────────
cat("--- Step 0: check_tools() ---\n")
check_tools()

# ── Step 0b: recommend GFF parameters for each species ────────────────────────
cat("\n--- Step 0b: recommend_parse_gff() ---\n")
rec_ath <- recommend_parse_gff(ath_gff, ath_fasta)
rec_vvi <- recommend_parse_gff(vvi_gff, vvi_fasta)

# ── Step 1: parse GFF → gene lists ────────────────────────────────────────────
cat("\n--- Step 1: parse_gff() ---\n")

# ath: nuclear chromosomes only (exclude ChrC and ChrM)
ath_lists <- parse_gff(
  gff_file     = ath_gff,
  output_dir   = file.path(work, "ath_lists"),
  genome_name  = "ath",
  feature_type = rec_ath$feature_type,
  id_attribute = rec_ath$id_attribute,
  chromosomes  = paste0("Chr", 1:5)
)

# vvi: chromosomes 1–19
vvi_lists <- parse_gff(
  gff_file     = vvi_gff,
  output_dir   = file.path(work, "vvi_lists"),
  genome_name  = "vvi",
  feature_type = rec_vvi$feature_type,
  id_attribute = rec_vvi$id_attribute,
  chromosomes  = paste0("chr", 1:19)
)

cat("ath gene list files:\n"); print(ath_lists)
cat("vvi gene list files:\n"); print(vvi_lists)

# ── Step 2: all-vs-all DIAMOND (both species together) ────────────────────────
cat("\n--- Step 2: run_diamond() ---\n")
blast_file <- file.path(work, "ath_vvi.blast")
run_diamond(
  fasta_file  = c(ath_fasta, vvi_fasta),
  output_file = blast_file,
  threads     = 4
)

# ── Step 3: cluster into gene families via MCL ────────────────────────────────
cat("\n--- Step 3: blast_to_families() ---\n")
families_file <- file.path(work, "families.txt")
blast_to_families(
  blast_file      = blast_file,
  output_file     = families_file,
  gene_list_files = c(unname(ath_lists), unname(vvi_lists))
)

# ── Step 4: write config ───────────────────────────────────────────────────────
cat("\n--- Step 4: write_iadhore_config() ---\n")
config_file <- file.path(work, "config.ini")
output_path <- file.path(work, "output")
write_iadhore_config(
  genomes     = list(ath = ath_lists, vvi = vvi_lists),
  blast_table = families_file,
  table_type  = "family",
  output_path = output_path,
  file        = config_file
)
cat("Config written to:", config_file, "\n")

# ── Step 5: run i-ADHoRe ──────────────────────────────────────────────────────
cat("\n--- Step 5: run_iadhore() ---\n")
run_iadhore(config_file)

# ── Step 6: read output ────────────────────────────────────────────────────────
cat("\n--- Step 6: read_iadhore_output() ---\n")
results <- read_iadhore_output(output_path)
cat("Tables loaded:", paste(names(results), collapse = ", "), "\n\n")
cat("multiplicons (first 6):\n");  print(head(results$multiplicons))
cat("\nanchorpoints (first 6):\n"); print(head(results$anchorpoints))

# ── Step 7: statistics ─────────────────────────────────────────────────────────
cat("\n--- Step 7: iadhore_summary() ---\n")
iadhore_summary(results)

cat("\n--- colinear_portions() ---\n")
cp <- colinear_portions(results)
print(cp)

cat("\n--- multiplicated_portions() ---\n")
mp <- multiplicated_portions(results)
print(mp)

cat("\n--- homolog_groups() ---\n")
hg <- homolog_groups(results)
cat("Total genes in homolog groups:", nrow(hg), "\n")
cat("Group size distribution:\n"); print(table(hg$group_size))
cat("Largest groups:\n")
print(head(hg[order(-hg$group_size, hg$homolog_group), ]))

# ── Step 7b: PAR analysis ─────────────────────────────────────────────────────
cat("\n--- Step 7b: find_pars() ath vs vvi ---\n")
pars_av <- find_pars(results, genome_x = "ath", genome_y = "vvi")
cat("PAR pairs (first 6):\n"); print(head(pars_av$pairs))

cat("\n--- Step 7c: find_pars() ath intra ---\n")
pars_ath <- find_pars(results, genome_x = "ath")

# ── Step 8: visualisations ────────────────────────────────────────────────────
cat("\n--- Step 8: visualisations ---\n")

ap_tab <- table(results$anchorpoints$is_real_anchorpoint)
cat("Anchorpoint counts by is_real_anchorpoint:\n"); print(ap_tab)

pdf(file.path(work, "plots.pdf"), width = 10, height = 8)

# 8a. Inter-genomic dot plot: ath vs vvi  (with PAR highlights)
cat("  plot_dotplot() ath vs vvi ...\n")
plot_dotplot(results,
             genome_x       = "ath", genome_y = "vvi",
             chr_order      = "natural",
             highlight_pars = pars_av)

# 8b. Intra-genomic dot plot: ath vs ath  (with PAR highlights)
cat("  plot_dotplot() ath intra ...\n")
plot_dotplot(results,
             genome_x       = "ath", genome_y = "ath",
             chr_order      = "natural",
             highlight_pars = pars_ath)

# 8c. Intra-genomic dot plot: vvi vs vvi
cat("  plot_dotplot() vvi intra ...\n")
plot_dotplot(results,
             genome_x  = "vvi", genome_y = "vvi",
             chr_order = "natural")

# 8d. Genome overview: ath
cat("  plot_genome_overview() ath ...\n")
plot_genome_overview(results, genome = "ath", chr_order = "natural")

# 8e. Genome overview: vvi
cat("  plot_genome_overview() vvi ...\n")
plot_genome_overview(results, genome = "vvi", chr_order = "natural")

# 8f. Track diagram for the first non-redundant multiplicon
first_mult <- results$multiplicons$id[results$multiplicons$is_redundant == 0][1]
cat(sprintf("  plot_multiplicon() for multiplicon %d ...\n", first_mult))
plot_multiplicon(results, multiplicon_id = first_mult)

dev.off()

# ── Step 9: PAR dot plots (Fig. S2 style) ─────────────────────────────────────
cat("\n--- Step 9: plot_par() ---\n")
pdf(file.path(work, "par_plots_ath_vvi.pdf"), width = 12, height = 10)
plot_par(results, pars_av)
dev.off()
cat("PAR plots written to:", file.path(work, "par_plots_ath_vvi.pdf"), "\n")

pdf(file.path(work, "par_plots_ath_intra.pdf"), width = 12, height = 10)
plot_par(results, pars_ath)
dev.off()
cat("PAR plots (ath intra) written to:", file.path(work, "par_plots_ath_intra.pdf"), "\n")

cat("\n=== All steps completed ===\n")
cat("Output:  ", output_path, "\n")
cat("Plots:   ", file.path(work, "plots.pdf"), "\n\n")
