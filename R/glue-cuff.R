
#' cufflinks commands
#' @importFrom glue glue
#' @inheritParams param_general
#' @param ref_gtf_path .tar file path
#' @param use_RABT .tar file path
#' @export
glue_cufflinks <-
  function(
    head_label,
    in_dir,
    out_dir = "./cuff/1_links",
    core_num = 2,
    ref_gtf_path = "./TAIR10_GFF3_genes_transposons.gtf",
    use_RABT = T
  ){
    le <- "\\"
    gtf_option <- ifelse(use_RABT, "-g", "-G")
      glue("

mkdir -p {out_dir}
cufflinks {le}
  -o {out_dir} {le}
  -p {core_num} {le}
  {gtf_option} {ref_gtf_path} {le}
  {in_dir}/{head_label}.sort.bam

mv {out_dir}/genes.fpkm_tracking {out_dir}/{head_label}_genes.fpkm_tracking
mv {out_dir}/isoforms.fpkm_tracking {out_dir}/{head_label}_isoforms.fpkm_tracking
mv {out_dir}/skipped.gtf {out_dir}/{head_label}_skipped.gtf
mv {out_dir}/transcripts.gtf {out_dir}/{head_label}_transcripts.gtf

    ")
  }


# glue_hoge <-
#   function(
#     head_label,
#     in_dir,
#     out_dir = "./cuff/1_links",
#     core_num = 2,
#     ref_gtf_path = "./TAIR10_GFF3_genes_transposons.gtf",
#     use_RABT = T
#   ){
#     le <- "\\"
#     gtf_option <- ifelse(use_RABT, "-g", "-G")
#     glue("
# mkdir -p {out_dir}
# files = ({bam_files})
# for temp_file in ${{files}}; do
#
# cufflinks {le}
#   -o {out_dir} {le}
#   -p {core_num} {ref_gtf_path} {le}
#   {in_dir}/${{temp_files}}.sort.bam
#
# done
#     ")
#   }
#
# glue_hoge(paste0("Motomura", 1:8), "hoge")


#' cuffmerge commands
#' @importFrom glue glue
#' @inheritParams param_general
#' @param cuff_gtf_dir .tar file path
#' @param ref_gtf_path .tar file path
#' @export
glue_cuffmerge <-
  function(
    cuff_gtf_dir = "./cuff/1_links",
    out_dir = "./cuff/2_merge",
    core_num = 2,
    ref_gtf_path = "./TAIR10_GFF3_genes_transposons.gtf"
  ){
    le <- "\\"
    glue("

find {cuff_gtf_dir} | grep transcripts.gtf > {cuff_gtf_dir}/assemblies.txt
mkdir -p {out_dir}
cuffmerge {le}
  -o {out_dir} {le}
  -g {ref_gtf_path} {le}
  -s allChr.fas {le}
  -p {core_num} {le}
  {cuff_gtf_dir}/assemblies.txt
rm {cuff_gtf_dir}/assemblies.txt

         ")
  }

#' cuffcompare commands
#' @importFrom glue glue
#' @inheritParams param_general
#' @param merged_gtf_path .tar file path
#' @param ref_gtf_path .tar file path
#' @param use_RABT .tar file path
#' @export
glue_cuffcompare <-
  function(
    merged_gtf_path = "./cuff/2_merge/merged.gtf",
    out_dir = "./cuff/3_merge",
    ref_gtf_path = "./TAIR10_GFF3_genes_transposons.gtf",
    use_RABT = T
  ){
    le <- "\\"
    gtf_option <- ifelse(use_RABT, "-g", "-G")
    glue("

find {cuff_gtf_dir} | grep transcripts.gtf > {cuff_gtf_dir}/assemblies.txt
mkdir -p {out_dir}
cuffcompare {le}
  -s allChr.fas {le}
  -r {ref_gtf_path}
  {merged_gtf_path}

         ")
  }

#' cuffquant commands
#' @importFrom glue glue
#' @inheritParams param_general
#' @param merged_gtf_path .tar file path
#' @export
glue_cuffquant <-
  function(
    head_label,
    in_dir,
    out_dir = "./cuff/4_quant",
    core_num = 2,
    merged_gtf_path = "./cuff/2_merge/merged.gtf"
  ){
    le <- "\\"
    glue("

mkdir {out_dir}
cuffquant{le}
  -o {out_dir} {le}
  -p {core_num} {le}
  {merged_gtf_path} {le}
  {in_dir}/{head_label}.sort.bam
mv {out_dir}/abundances.cxb {out_dir}/{head_label}.cxb

         ")
  }


#' cuffdiff commands
#' @importFrom glue glue
#' @importFrom stringr str_c
#' @inheritParams param_general
#' @param head_labels .tar file path
#' @param merged_gtf_path .tar file path
#' @export
glue_cuffdiff <-
  function(
    head_labels,
    in_dir = "./cuff/4_quant",
    out_dir = "./cuff/5_diff",
    merged_gtf_path = "./cuff/2_merge/merged.gtf",
    core_num = 2
  ){
    . <- NULL
    le <- "\\"
    group <- names(head_labels)

    input_list <- list()
    for(i in unique(group)){
      temp <- head_labels[names(head_labels) == i]
      input_list[[i]] <-
        str_c(temp, ".cxb") %>%
        str_c(in_dir, "/", .) %>%
        str_c(collapse = ",")
    }
    label_string <- names(input_list) %>% str_c(collapse = ",")
    input_string <- input_list %>% str_c(collapse = " \\\n")

    glue("

mkdir {out_dir}
cuffdiff {le}
  -o {out_dir} {le}
  -p {core_num} {le}
  -L {label_string} {le}
  {merged_gtf_path} {le}
  {input_string}

    ")
  }

#' cuffnorm commands
#' @importFrom glue glue
#' @importFrom stringr str_c
#' @inheritParams param_general
#' @param head_labels .tar file path
#' @param merged_gtf_path .tar file path
#' @export
glue_cuffnorm <-
  function(
    head_labels,
    in_dir = "./cuff/4_quant",
    out_dir = "./cuff/6_norm",
    merged_gtf_path = "./cuff/2_merge/merged.gtf",
    core_num = 2
  ){
    . <- NULL
    le <- "\\"
    group <- names(head_labels)

    input_list <- list()
    for(i in unique(group)){
      temp <- head_labels[names(head_labels) == i]
      input_list[[i]] <-
        str_c(temp, ".cxb") %>%
        str_c(in_dir, "/", .) %>%
        str_c(collapse = ",")
    }
    label_string <- names(input_list) %>% str_c(collapse = ",")
    input_string <- input_list %>% str_c(collapse = " \\\n")

    glue("

mkdir {out_dir}
cuffnorm {le}
  -o {out_dir} {le}
  -p {core_num} {le}
  -L {label_string} {le}
  {merged_gtf_path} {le}
  {input_string}

    ")
  }


#' pipeline for cufflinks
#' @importFrom glue glue
#' @importFrom purrr walk
#' @inheritParams param_general
#' @param head_labels .tar file path
#' @param ref_gtf_path .tar file path
#' @export
pipeline_cuffdiff <-
  function(
    head_labels,
    in_dir,
    core_num = 2,
    ref_gtf_path = "./TAIR10_GFF3_genes_transposons.gtf"
  ){
    . <- NULL

    walk(
      head_labels,
      ~ print(glue_cufflinks(
        head_label = .,
        in_dir = in_dir,
        core_num = core_num,
        ref_gtf_path = ref_gtf_path
        ))
    )

    print(
      glue_cuffmerge(
        core_num = core_num,
        ref_gtf_path = ref_gtf_path
      )
    )

    walk(
      head_labels,
      ~ print(glue_cuffquant(head_label = ., in_dir = in_dir, core_num = core_num))
    )

    print(
      glue_cuffdiff(
        head_labels = head_labels,
        core_num = core_num
        )
    )

    print(
      glue_cuffnorm(
        head_labels = head_labels,
        core_num = core_num
        )
    )
  }


#' pipeline for cufflinks
#' @importFrom glue glue
#' @importFrom purrr walk
#' @inheritParams param_general
#' @param head_labels .tar file path
#' @param ref_gtf_path .tar file path
#' @export
pipeline_cuffdiff_noRABT <-
  function(
    head_labels,
    in_dir,
    core_num = 2,
    ref_gtf_path = "./TAIR10_GFF3_genes_transposons.gtf"
  ){
    . <- NULL

    walk(
      head_labels,
      ~ print(glue_cufflinks(
        head_label = .,
        in_dir = in_dir,
        out_dir = "./cuff_noRABT/1_links",
        core_num = core_num,
        ref_gtf_path = ref_gtf_path,
        use_RABT = F
      ))
    )

    walk(
      head_labels,
      ~ print(glue_cuffquant(
        head_label = .,
        in_dir = in_dir,
        out_dir = "./cuff_noRABT/4_quant",
        core_num = core_num,
        merged_gtf_path = ref_gtf_path
      ))
    )

    print(
      glue_cuffdiff(
        head_labels = head_labels,
        in_dir = "./cuff_noRABT/4_quant",
        out_dir = "./cuff_noRABT/5_diff",
        merged_gtf_path = ref_gtf_path,
        core_num = core_num
        )
    )

    print(
      glue_cuffnorm(
        head_labels = head_labels,
        in_dir = "./cuff_noRABT/4_quant",
        out_dir = "./cuff_noRABT/6_norm",
        merged_gtf_path = ref_gtf_path,
        core_num = core_num
        )
    )
  }




