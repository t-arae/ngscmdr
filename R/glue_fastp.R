
#' Filter fastq file by quality
#' @inheritParams param_general
#' @export
cmd_set_fastp_qf <-
  function(
    inf_list,
    out_dir = "fastq_qf"
  ){
    infs <- purrr::map(inf_list, "inf")
    labels <- purrr::map_chr(inf_list, "label")
    new_cmdopt(
      inf = infs,
      label = labels,
      out_dir = out_dir,
      #outf = file.path(out_dir, paste0(label, ".", get_file_ext()$fastq)),
      cmd_id = "fastp_qf"
    )
  }

cmd_get.fastp_qf <-
  function(x){
    cmd_list <- list()
    for(i in x){
      temp_env <- i
      temp_env$fq_ext <- get_file_ext()$fastq
      temp_env$le <- "\\"
      if(temp_env$is_pe){
      cmd_list <- glue::glue("
./fastp -q 20 -u 15 {le}
  -i {inf[1]} {le}
  -I {inf[2]} {le}
  -o {out_dir}/{label}_1.{fq_ext}.gz {le}
  -O {out_dir}/{label}_2.{fq_ext}.gz
    ", .envir = as.environment(x))

      }
      cmd_list <- glue::glue("
fastp -q 20 -u 15 {le}
  -i {inf} {le}
  -o {out_dir}/{label}.{fq_ext}.gz
    ", .envir = as.environment(x)
      )
      invisible(list(inf = x$inf, outf = x$outf))

    }
    x$fq_ext <- get_file_ext()$fastq
    x$le <- "\\"
    print(glue::glue("
mkdir {out_dir}
fastp -q 20 -u 15 {le}
  -i {inf} {le}
  -o {out_dir}/{label}.{fq_ext}.gz
    ", .envir = as.environment(x)
    ))
    invisible(list(inf = x$inf, outf = x$outf))
  }
# cmd_set_fastp_qf("hoge.fq.gz", "hoge") %>% cmd_get()

#' Filter fastq file by quality
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastp_qf <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./fastq_qf"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

mkdir {out_dir}
fastp -q 20 -u 15 {le}
  -i {in_dir}/{head_label}.{fq_ext}.gz {le}
  -o {out_dir}/{head_label}.{fq_ext}.gz

    ")
  }

#' Filter fastq file by quality
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastp_qf_pe <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./fastq_qf"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

mkdir {out_dir}
./fastp -q 20 -u 15 {le}
  -i {in_dir}/{head_label}_1.{fq_ext}.gz {le}
  -I {in_dir}/{head_label}_2.{fq_ext}.gz {le}
  -o {out_dir}/{head_label}_1.{fq_ext}.gz {le}
  -O {out_dir}/{head_label}_2.{fq_ext}.gz

    ")
  }
