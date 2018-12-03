
#' featureCounts
#' @importFrom glue glue
#' @inheritParams param_general
#' @param gtf_path .tar file path
#' @param feat_type specify feature type
#' @export
glue_se_featurecounts <-
  function(
    head_label,
    in_dir,
    out_dir = "./readcount",
    core_num = 2,
    gtf_path = "./TAIR10_GFF3_genes_transposons.gtf",
    feat_type = "exon"
  ){
    le <- "\\"
    pe_flag <- ""
    glue("

mkdir {out_dir}
featureCounts {le}{pe_flag}
  -M {le}
  --fraction {le}
  -t {feat_type} {le}
  -g gene_id {le}
  -T {core_num} {le}
  -a {gtf_path} {le}
  -o {out_dir}/{head_label}_gene_counts.txt {le}
  {in_dir}/{head_label}.sort.bam {le}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_gene_counts.log {le}
  3>&2 2>&1 1>&3

featureCounts {le}{pe_flag}
  -M {le}
  --fraction {le}
  -t {feat_type} {le}
  -g transcript_id {le}
  -O {le}
  -T {core_num} {le}
  -a {gtf_path} {le}
  -o {out_dir}/{head_label}_transcript_counts.txt {le}
  {in_dir}/{head_label}.sort.bam {le}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_transcript_counts.log {le}
  3>&2 2>&1 1>&3

         ")
  }

#' featureCounts for Pair-end reads
#' @importFrom glue glue
#' @inheritParams param_general
#' @param gtf_path .tar file path
#' @param feat_type specify feature type
#' @export
glue_pe_featurecounts <-
  function(
    head_label,
    in_dir,
    out_dir = "./readcount",
    core_num = 2,
    gtf_path = "./TAIR10_GFF3_genes_transposons.gtf",
    feat_type = "exon"
  ){
    le <- "\\"
    pe_flag <- "\n  -p \\"
    glue("
mkdir {out_dir}
featureCounts {le}{pe_flag}
  -M {le}
  --fraction {le}
  -t {feat_type} {le}
  -g gene_id {le}
  -T {core_num} {le}
  -a {gtf_path} {le}
  -o {out_dir}/{head_label}_gene_counts.txt {le}
  {in_dir}/{head_label}.sort.bam {le}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_gene_counts.log {le}
  3>&2 2>&1 1>&3

featureCounts {le}{pe_flag}
  -M {le}
  --fraction {le}
  -t {feat_type} {le}
  -g transcript_id {le}
  -O {le}
  -T {core_num} {le}
  -a {gtf_path} {le}
  -o {out_dir}/{head_label}_transcript_counts.txt {le}
  {in_dir}/{head_label}.sort.bam {le}
  3>&2 2>&1 1>&3 | tee {out_dir}/{head_label}_transcript_counts.log {le}
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
#' @inheritParams param_general
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
