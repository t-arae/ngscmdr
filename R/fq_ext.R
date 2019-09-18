

#' Set file extensions
#' @param fastq_ext fastq file extension
#' @param bam_ext bam file extension
#' @param bigwiggle_ext bigwiggle file extension
#' @export
set_file_ext <-
  function(
    fastq_ext = "fq",
    bam_ext = "bam",
    bigwiggle_ext = "bw"
  ){
    options("ngscmdr_file_ext" =
              list(
                "fastq" = fastq_ext,
                "bam" = bam_ext,
                "bigwig" = bigwiggle_ext
              )
    )
  }

#' Get file extensions
#' @export
get_file_ext <-
  function(){
    getOption("ngscmdr_file_ext")
  }

set_file_ext()
