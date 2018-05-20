
#' featureCounts
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @param core_num .tar file path
#' @param gtf_path .tar file path
#' @export
glue_se_featurecounts <-
  function(
    head_label,
    in_dir,
    out_dir = "./readcount",
    core_num = 2,
    gtf_path = "./TAIR10_GFF3_genes_transposons.gtf"
  ){
    line_end <- "\\"
    pe_flag <- ""
    glue("

mkdir {out_dir}
featureCounts {line_end}{pe_flag}
  -t exon {line_end}
  -g gene_id {line_end}
  -T {core_num} {line_end}
  -a {gtf_path} {line_end}
  -o {out_dir}/{head_label}_gene_counts.txt {line_end}
  {in_dir}/{head_label}.sort.bam {line_end}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_gene_counts.log {line_end}
  3>&2 2>&1 1>&3

featureCounts {line_end}{pe_flag}
  -t exon {line_end}
  -g transcript_id {line_end}
  -O {line_end}
  -T {core_num} {line_end}
  -a {gtf_path} {line_end}
  -o {out_dir}/{head_label}_transcript_counts.txt {line_end}
  {in_dir}/{head_label}.sort.bam {line_end}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_transcript_counts.log {line_end}
  3>&2 2>&1 1>&3

         ")
  }

#' featureCounts for Pair-end reads
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @param core_num .tar file path
#' @param gtf_path .tar file path
#' @export
glue_pe_featurecounts <-
  function(
    head_label,
    in_dir,
    out_dir = "./readcount",
    core_num = 2,
    gtf_path = "./TAIR10_GFF3_genes_transposons.gtf"
  ){
    line_end <- "\\"
    pe_flag <- "\n  -p \\"
    glue("

mkdir {out_dir}
featureCounts {line_end}{pe_flag}
  -t exon {line_end}
  -g gene_id {line_end}
  -T {core_num} {line_end}
  -a {gtf_path} {line_end}
  -o {out_dir}/{head_label}_gene_counts.txt {line_end}
  {in_dir}/{head_label}.sort.bam {line_end}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_gene_counts.log {line_end}
  3>&2 2>&1 1>&3

featureCounts {line_end}{pe_flag}
  -t exon {line_end}
  -g transcript_id {line_end}
  -O {line_end}
  -T {core_num} {line_end}
  -a {gtf_path} {line_end}
  -o {out_dir}/{head_label}_transcript_counts.txt {line_end}
  {in_dir}/{head_label}.sort.bam {line_end}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_transcript_counts.log {line_end}
  3>&2 2>&1 1>&3

         ")
  }


#' Merging read count tsv file from featureCounts
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom readr read_tsv
#' @importFrom dplyr bind_cols
#' @importFrom tibble as.tibble
#' @importFrom stringr str_extract
#' @importFrom stringr str_sub
#' @importFrom magrittr %>%
#' @param in_dir .tar file path
#' @export
merge_featurecount_output <-
  function(
    in_dir = "./readcount"
  ){
    . <- NULL
    read_and_merge_df <-
      function(fpath_li){
        inf <- map(fpath_li, ~ suppressMessages(read_tsv(., skip = 1)))

        df_feature <- inf[[1]][,1:6]
        df_count <- map(inf, ~ .[[7]])
        names(df_count) <-
          map_chr(inf, ~ colnames(.)[7]) %>%
          str_extract("[^/]+?[.]sort[.]bam$") %>%
          str_sub(1, -10)

        bind_cols(df_feature, as.tibble(df_count))
      }

    list(
      "count_by_gene" =
        list.files(path = in_dir, pattern = "gene_counts.txt$", full.names = T) %>%
        read_and_merge_df,

      "count_by_transcript" =
        list.files(path = in_dir, pattern = "transcript_counts.txt$", full.names = T) %>%
        read_and_merge_df
    )
  }
