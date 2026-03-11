# iadhoreR

An R interface to [i-ADHoRe 3.0](https://github.com/VIB-PSB/i-ADHoRe) for
detecting collinear (syntenic) regions within and between genomes.

## What it does

i-ADHoRe identifies conserved gene order across chromosomes — evidence of
ancient whole-genome duplications or shared ancestry between species. iadhoreR
handles the full workflow from raw annotation files to parsed results:

```
GFF + FASTA  →  parse_gff()
                run_diamond()  →  blast_to_families()
                                  write_iadhore_config()
                                  run_iadhore()
                                  read_iadhore_output()
```

## Installation

Follow these steps in order.

### Step 1 — Install R

Download and install R (≥ 4.0) from https://cran.r-project.org.

[RStudio](https://posit.co/download/rstudio-desktop/) is recommended as an
IDE but not required.

### Step 2 — Install conda

If you do not already have conda, install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) (lightweight) or
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) (faster
solver, recommended for bioinformatics).

### Step 3 — Install the external tools

All three tools — **i-ADHoRe**, **DIAMOND**, and **MCL** — can be installed
via conda:

**Option A: dedicated environment (recommended)**

```bash
conda env create -f https://raw.githubusercontent.com/lizhencmb/iadhoreR/main/inst/conda/environment.yml
conda activate iadhoreR
```

**Option B: install into an existing environment**

```bash
conda activate my_existing_env
conda install -c lizhencmb -c bioconda -c conda-forge i-adhore diamond mcl
```

> **Important:** always activate the conda environment *before* opening R or
> RStudio, so that the tools are on the PATH. On macOS/Linux, open a terminal,
> activate the environment, then launch R or RStudio from that same terminal:
>
> ```bash
> conda activate iadhoreR
> open -a RStudio   # macOS
> rstudio &         # Linux
> ```

### Step 4 — Install iadhoreR

With the conda environment active, open R and run:

```r
install.packages("remotes")
remotes::install_github("lizhencmb/iadhoreR", build_vignettes = TRUE)
```

### Step 5 — Verify

```r
library(iadhoreR)
check_tools()
#> External tool status:
#>   i-adhore     [OK]     /path/to/conda/envs/iadhoreR/bin/i-adhore
#>   diamond      [OK]     /path/to/conda/envs/iadhoreR/bin/diamond
#>   mcl          [OK]     /path/to/conda/envs/iadhoreR/bin/mcl
#>   mcxload      [OK]     /path/to/conda/envs/iadhoreR/bin/mcxload
#>   mcxdump      [OK]     /path/to/conda/envs/iadhoreR/bin/mcxdump
#> All tools found. You are ready to use iadhoreR.
```

If any tool shows `[MISSING]`, run `setup_instructions()` in R for
troubleshooting guidance.

---

## Windows users (WSL2)

i-ADHoRe and MCL do not have native Windows builds. The recommended approach
is **Windows Subsystem for Linux 2 (WSL2)**, which runs a full Linux
environment inside Windows and is fully supported.

**Set up WSL2:**

1. Open PowerShell as Administrator and run:
   ```powershell
   wsl --install
   ```
   This installs WSL2 with Ubuntu. Restart your computer when prompted.

2. Open the Ubuntu app and install Miniconda inside WSL2:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   # follow the prompts, then restart the shell
   ```

3. Install the tools and R inside WSL2:
   ```bash
   conda env create -f https://raw.githubusercontent.com/lizhencmb/iadhoreR/main/inst/conda/environment.yml
   conda activate iadhoreR
   conda install -c conda-forge r-base
   ```

4. Inside WSL2 R, install iadhoreR:
   ```r
   install.packages("remotes")
   remotes::install_github("lizhencmb/iadhoreR", build_vignettes = TRUE)
   ```

**Using RStudio on Windows with WSL2:**

Install [RStudio Desktop](https://posit.co/download/rstudio-desktop/) on
Windows. RStudio automatically detects WSL2 and can use R installed inside it
— see the
[Posit WSL2 guide](https://support.posit.co/hc/en-us/articles/360049776974)
for setup instructions.

---

## Quick start

### Step 0 — Check your GFF format

Before parsing, inspect your annotation file to find the right feature type
and attribute key for your data:

```r
library(iadhoreR)

# See all feature types and attribute keys in the file
inspect_gff("species1.gff3")

# Auto-match GFF IDs against your protein FASTA and get recommended parameters
recommend_parse_gff("species1.gff3", "species1_proteins.fasta")
#> Best match: feature_type = "mRNA", id_attribute = "ID"
#>   Matched: 27655 / 27655
```

Use the `feature_type` and `id_attribute` values returned in the examples below.

---

### Example 1 — Single species (intra-genomic duplications)

```r
work <- "my_analysis"

# 1. Parse GFF into gene lists (one file per chromosome)
sp1_lists <- parse_gff("species1.gff3",
                       output_dir   = file.path(work, "sp1_lists"),
                       genome_name  = "sp1",
                       feature_type = "mRNA",
                       id_attribute = "ID")

# 2. All-vs-all protein similarity search
run_diamond("species1_proteins.fasta",
            output_file = file.path(work, "sp1.blast"),
            threads = 8)

# 3. Cluster into gene families
blast_to_families(
  blast_file      = file.path(work, "sp1.blast"),
  output_file     = file.path(work, "families.txt"),
  gene_list_files = unname(sp1_lists)
)

# 4. Write config and run i-ADHoRe
write_iadhore_config(
  genomes     = list(sp1 = sp1_lists),
  blast_table = file.path(work, "families.txt"),
  table_type  = "family",
  output_path = file.path(work, "output"),
  file        = file.path(work, "config.ini")
)
run_iadhore(file.path(work, "config.ini"))

# 5. Read results
results <- read_iadhore_output(file.path(work, "output"))
head(results$multiplicons)   # syntenic regions
head(results$anchorpoints)   # homologous gene pairs
```

---

### Example 2 — Two species (inter-genomic synteny)

```r
work <- "my_analysis"

# 1. Parse GFF for each species
sp1_lists <- parse_gff("species1.gff3",
                       output_dir   = file.path(work, "sp1_lists"),
                       genome_name  = "sp1",
                       feature_type = "mRNA",
                       id_attribute = "ID")
sp2_lists <- parse_gff("species2.gff3",
                       output_dir   = file.path(work, "sp2_lists"),
                       genome_name  = "sp2",
                       feature_type = "mRNA",
                       id_attribute = "ID")

# 2. All-vs-all search across both species (pass both FASTAs)
run_diamond(c("species1_proteins.fasta", "species2_proteins.fasta"),
            output_file = file.path(work, "all_vs_all.blast"),
            threads = 8)

# 3. Cluster into gene families (include both species' gene lists)
blast_to_families(
  blast_file      = file.path(work, "all_vs_all.blast"),
  output_file     = file.path(work, "families.txt"),
  gene_list_files = c(unname(sp1_lists), unname(sp2_lists))
)

# 4. Write config and run i-ADHoRe
write_iadhore_config(
  genomes     = list(sp1 = sp1_lists, sp2 = sp2_lists),
  blast_table = file.path(work, "families.txt"),
  table_type  = "family",
  output_path = file.path(work, "output"),
  file        = file.path(work, "config.ini")
)
run_iadhore(file.path(work, "config.ini"))

# 5. Read results
results <- read_iadhore_output(file.path(work, "output"))
head(results$multiplicons)   # syntenic regions
head(results$anchorpoints)   # homologous gene pairs
```

## Full tutorial

A step-by-step vignette using bundled *Arabidopsis* and *Vitis* example data
is available after installation:

```r
vignette("iadhoreR-tutorial", package = "iadhoreR")
```

## Key functions

| Function | Description |
|----------|-------------|
| `check_tools()` | Verify all external tools are on PATH |
| `setup_instructions()` | Print conda installation commands |
| `inspect_gff()` | Explore feature types and attributes in a GFF file |
| `recommend_parse_gff()` | Auto-detect best GFF parameters for your FASTA |
| `parse_gff()` | Create i-ADHoRe gene list files from a GFF |
| `run_diamond()` | All-vs-all protein similarity search |
| `blast_to_families()` | Cluster BLAST results into gene families via MCL |
| `parse_blast()` | Filter BLAST results into a gene-pair table |
| `write_iadhore_config()` | Write i-ADHoRe configuration file |
| `run_iadhore()` | Run i-ADHoRe |
| `read_iadhore_output()` | Read all output tables into a named list |

## License

MIT
