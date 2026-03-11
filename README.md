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

## Requirements

iadhoreR calls three external tools that must be installed via
[conda](https://docs.conda.io):

```bash
conda install -c bioconda -c conda-forge i-adhore diamond mcl
```

Or use the bundled environment file for a self-contained setup:

```bash
# Get the environment file path from R
Rscript -e "system.file('conda', 'environment.yml', package='iadhoreR')"

# Then create the environment
conda env create -f /path/to/environment.yml
conda activate iadhoreR
```

> **Important:** activate the conda environment *before* launching R so that
> the tools are on the PATH.

## Installation

```r
# install.packages("remotes")
remotes::install_github("your-org/iadhoreR")
```

After installation, verify your setup:

```r
library(iadhoreR)
check_tools()
```

## Quick start

```r
library(iadhoreR)

# Paths to your input files
gff1  <- "species1.gff3"
gff2  <- "species2.gff3"
fasta1 <- "species1_proteins.fasta"
fasta2 <- "species2_proteins.fasta"
work  <- "my_analysis"

# 1. Parse GFF files into gene lists
sp1_lists <- parse_gff(gff1, output_dir = file.path(work, "sp1_lists"),
                       genome_name = "sp1", feature_type = "mRNA",
                       id_attribute = "ID")
sp2_lists <- parse_gff(gff2, output_dir = file.path(work, "sp2_lists"),
                       genome_name = "sp2", feature_type = "mRNA",
                       id_attribute = "ID")

# 2. All-vs-all protein similarity search
run_diamond(c(fasta1, fasta2),
            output_file = file.path(work, "all_vs_all.blast"),
            threads = 8)

# 3. Cluster into gene families
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
run_iadhore(file.path(work, "config.ini"), threads = 8)

# 5. Read results
results <- read_iadhore_output(file.path(work, "output"))
head(results$multiplicons)   # syntenic regions
head(results$anchorpoints)   # homologous gene pairs
```

Not sure which GFF features and attributes to use? Start with:

```r
inspect_gff("species1.gff3")
recommend_parse_gff("species1.gff3", "species1_proteins.fasta")
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
