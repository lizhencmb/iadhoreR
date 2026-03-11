#' Read a tab-delimited i-ADHoRe output table
#' @keywords internal
read_iadhore_table <- function(file) {
  if (!file.exists(file)) {
    stop("Output file not found: ", file)
  }
  utils::read.table(file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
}


#' Read the multiplicons table
#'
#' Reads `multiplicons.txt` from an i-ADHoRe output directory.
#'
#' Columns: id, genome_x, list_x, parent, genome_y, list_y, level,
#' number_of_anchorpoints, profile_length, begin_x, end_x, begin_y, end_y,
#' is_redundant
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_multiplicons <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "multiplicons.txt"))
}


#' Read the baseclusters table
#'
#' Reads `baseclusters.txt` from an i-ADHoRe output directory.
#'
#' Columns: id, multiplicon, number_of_anchorpoints, orientation,
#' was_twisted, random_probability
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_baseclusters <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "baseclusters.txt"))
}


#' Read the anchorpoints table
#'
#' Reads `anchorpoints.txt` from an i-ADHoRe output directory.
#'
#' Columns: id, multiplicon, basecluster, gene_x, gene_y, coord_x, coord_y,
#' is_real_anchorpoint
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_anchorpoints <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "anchorpoints.txt"))
}


#' Read the segments table
#'
#' Reads `segments.txt` from an i-ADHoRe output directory.
#'
#' Columns: id, multiplicon, genome, list, first, last, order
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_segments <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "segments.txt"))
}


#' Read the genes table
#'
#' Reads `genes.txt` from an i-ADHoRe output directory.
#'
#' Columns: id, genome, list, coordinate, orientation, remapped_coordinate,
#' is_tandem, is_tandem_representative, tandem_representative, remapped
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_genes <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "genes.txt"))
}


#' Read the list_elements table
#'
#' Reads `list_elements.txt` from an i-ADHoRe output directory.
#'
#' Columns: id, segment, gene, position, orientation
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_list_elements <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "list_elements.txt"))
}


#' Read the clouds table
#'
#' Reads `clouds.txt` from an i-ADHoRe output directory (cloud/hybrid mode
#' only).
#'
#' Columns: id, genome_x, list_x, genome_y, list_y, number_of_anchorpoints,
#' cloud_density, dim_x, dim_y
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_clouds <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "clouds.txt"))
}


#' Read the cloud anchorpoints table
#'
#' Reads `cloudAP.txt` from an i-ADHoRe output directory (cloud/hybrid mode
#' only).
#'
#' Columns: CloudID, gene_x, gene_y, coord_x, coord_y
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A data frame
#' @export
read_cloud_anchorpoints <- function(output_dir) {
  read_iadhore_table(file.path(output_dir, "cloudAP.txt"))
}


#' Read all i-ADHoRe output tables
#'
#' Reads all standard output tables from an i-ADHoRe output directory into a
#' named list. Cloud tables are included only when present.
#'
#' @param output_dir Path to the i-ADHoRe output directory
#'
#' @return A named list with elements: multiplicons, baseclusters,
#'   anchorpoints, segments, genes, list_elements, and optionally clouds and
#'   cloud_anchorpoints
#' @export
#'
#' @examples
#' \dontrun{
#' results <- read_iadhore_output("output/arath/")
#' head(results$multiplicons)
#' head(results$anchorpoints)
#' }
read_iadhore_output <- function(output_dir) {
  result <- list(
    multiplicons   = read_multiplicons(output_dir),
    baseclusters   = read_baseclusters(output_dir),
    anchorpoints   = read_anchorpoints(output_dir),
    segments       = read_segments(output_dir),
    genes          = read_genes(output_dir),
    list_elements  = read_list_elements(output_dir)
  )

  clouds_file <- file.path(output_dir, "clouds.txt")
  cloud_ap_file <- file.path(output_dir, "cloudAP.txt")

  if (file.exists(clouds_file)) {
    result$clouds <- read_clouds(output_dir)
  }
  if (file.exists(cloud_ap_file)) {
    result$cloud_anchorpoints <- read_cloud_anchorpoints(output_dir)
  }

  result
}
