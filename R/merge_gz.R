
#' Merge some gzipped files to a gzipped file
#' @importFrom glue glue
#' @param gziped_files .tar file path
#' @param gzip_out .tar file path
#' @export
cmd_set_merge_gz <-
  function(gziped_files, gzip_out){
    new_cmdopt(inf = gziped_files, outf = gzip_out, cmd_id = "merge_gz")
  }

cmd_get.merge_gz <-
  function(x){
    x$inf <- paste(x$inf, collapse = " ")
    print(glue::glue("
cat {inf} > {outf} && rm {inf}
    ", .envir = as.environment(x)
    ))
    invisible(list(inf = x$inf, outf = x$outf))
  }
# cmd_set_merge_gz(c("hoge.gz", "hage.gz"), "hige.gz") %>% cmd_get()
