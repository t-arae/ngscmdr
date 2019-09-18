
#' Parameters for general use
#' @param in_dir A path of directory which contains input files
#' @param out_dir A path of directory for saving output files
#' @param head_label sample label
#' @param label sample label
#' @param core_num number of cpu core
param_general <-
  function(
    in_dir,
    out_dir,
    head_label,
    core_num
  ){}

#' Parameters for genome indexing commands
#' @param fasta_path genome fasta file path for indexing
#' @param idx_name index name
#' @param idx_dir directory path for the indexing output
param_genomeidx <-
  function(
    fasta_path,
    idx_name,
    idx_dir
  ){}


