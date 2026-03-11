#' Check that all required external tools are available
#'
#' Looks up i-ADHoRe, DIAMOND, and the MCL suite on the system PATH and prints
#' a status table. If any tool is missing, prints setup instructions.
#'
#' All three tools must be on the PATH before using [run_iadhore()],
#' [run_diamond()], or [blast_to_families()]. The easiest way to install them
#' is via conda. A ready-made environment file is bundled with the package:
#'
#' ```
#' # Option A — create a dedicated environment (recommended)
#' conda env create -f <path from setup_instructions()>
#' conda activate iadhoreR
#'
#' # Option B — install into an existing environment
#' conda install -c bioconda -c conda-forge i-adhore diamond mcl
#' ```
#'
#' Call [setup_instructions()] to print the exact commands for your system.
#'
#' @return Invisibly returns a named logical vector (TRUE = found).
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' check_tools()
check_tools <- function() {
  tools <- c("i-adhore", "diamond", "mcl", "mcxload", "mcxdump")
  paths <- Sys.which(tools)
  found <- nchar(paths) > 0

  cat("External tool status:\n")
  for (i in seq_along(tools)) {
    status <- if (found[i]) paste0("[OK]     ", paths[i]) else "[MISSING]"
    cat(sprintf("  %-12s %s\n", tools[i], status))
  }

  if (!all(found)) {
    cat("\nOne or more tools are missing. Run setup_instructions() for help.\n")
  } else {
    cat("\nAll tools found. You are ready to use iadhoreR.\n")
  }

  invisible(setNames(found, tools))
}


#' Print conda setup instructions for iadhoreR
#'
#' Prints the commands needed to create a conda environment with all tools
#' required by iadhoreR (i-ADHoRe, DIAMOND, MCL). Also prints the path to the
#' bundled `environment.yml` file that can be used directly with
#' `conda env create`.
#'
#' @return Invisibly returns the path to the bundled `environment.yml`.
#' @export
#'
#' @examples
#' setup_instructions()
setup_instructions <- function() {
  env_yml <- system.file("conda", "environment.yml",
                         package = "iadhoreR", mustWork = TRUE)

  cat("iadhoreR requires i-ADHoRe, DIAMOND, and MCL on your PATH.\n")
  cat("Install them via conda (https://docs.conda.io):\n\n")

  cat("  # Option A: create a dedicated environment (recommended)\n")
  cat("  conda env create -f", env_yml, "\n")
  cat("  conda activate iadhoreR\n\n")

  cat("  # Option B: install into your current environment\n")
  cat("  conda install -c bioconda -c conda-forge i-adhore diamond mcl\n\n")

  cat("After activating the environment, launch R from that same shell\n")
  cat("so it inherits the correct PATH.\n\n")
  cat("  # If using RStudio, launch it from the activated shell:\n")
  cat("  conda activate iadhoreR\n")
  cat("  open -a RStudio   # macOS\n")
  cat("  rstudio &         # Linux\n\n")

  cat("Verify the installation with:\n")
  cat("  library(iadhoreR)\n")
  cat("  check_tools()\n")

  invisible(env_yml)
}


#' Get i-ADHoRe executable path
#' @keywords internal
get_iadhore_bin <- function() {
  path <- Sys.which("i-adhore")
  if (path == "") {
    stop(
      "i-ADHoRe not found in PATH.\n",
      "Install it via conda:\n",
      "  conda install -c bioconda i-adhore"
    )
  }
  path
}

#' Run i-ADHoRe analysis
#'
#' @param config_file Path to i-ADHoRe configuration file. All paths inside
#'   the config (gene lists, blast table, output path) must be relative to the
#'   directory containing this file.
#' @param verbose Print i-ADHoRe output (default: TRUE)
#'
#' @return Invisibly returns the exit status (0 for success)
#' @export
#'
#' @examples
#' \dontrun{
#' run_iadhore("config.ini")
#' }
run_iadhore <- function(config_file, verbose = TRUE) {

  config_file <- normalizePath(config_file, mustWork = TRUE)
  config_dir  <- dirname(config_file)
  config_name <- basename(config_file)

  # get the correct executable file
  iadhore_bin <- get_iadhore_bin()

  # i-ADHoRe resolves all paths in the config relative to the working
  # directory, so we must run it from the config file's directory
  old_wd <- setwd(config_dir)
  on.exit(setwd(old_wd), add = TRUE)

  result <- system2(
    command = iadhore_bin,
    args    = config_name,
    stdout  = if (verbose) "" else FALSE,
    stderr  = if (verbose) "" else FALSE,
    wait    = TRUE
  )

  if (result == 11L) {
    # Status 11 = SIGSEGV: i-ADHoRe is known to segfault during teardown
    # after successfully writing all output. Warn but do not treat as failure.
    warning(
      "i-ADHoRe exited with status 11 (segmentation fault). ",
      "This is a known issue with the i-ADHoRe binary and typically occurs ",
      "after all output has been written. Check that your output files are ",
      "complete before proceeding.",
      call. = FALSE
    )
  } else if (result != 0L) {
    stop("i-ADHoRe exited with status ", result,
         ". Check the output above for error messages.", call. = FALSE)
  }

  return(invisible(result))
}
