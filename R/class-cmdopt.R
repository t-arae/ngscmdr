#' Create new cmdopt class object
#' @param ... ...
#' @param cmd_id command type specifier
new_cmdopt <-
  function(..., cmd_id = NULL){
    cmdopt_obj <- list(...)
    class(cmdopt_obj) <- c(cmd_id, "cmdopt", class(cmdopt_obj))
    cmdopt_obj
  }

#' Print cmdopt class object
#' @param x cmdopt class object
#' @param ... ...
#' @export
print.cmdopt <-
  function(x, ...){
    cat("cmdopt class object\n")
    cat("input file directory: ", x$in_dir, "\n")
    cat("input file: ", x$inf, "\n")
    cat("output file directory: ", x$out_dir, "\n")
    cat("output file name: ", x$outf, "\n")
  }

#' cmd_get methods for cmdopt class object
#' @param x cmdopt class object
#' @export
cmd_get <- function(x) UseMethod("cmd_get")


#' Extract files from arachive file (.tar)
#' @param archive_file .tar file path.
#' @export
cmd_set_extract_file_from_tar <-
  function(archive_file){
    extract_to <- dirname(archive_file)
    new_cmdopt(
      inf = archive_file,
      out_dir = extract_to,
      cmd_id = "extract_file_from_tar"
    )
  }

cmd_get.extract_file_from_tar <-
  function(x){
    glue::glue("
tar xvf {inf} -C {out_dir}
      ", .envir = as.environment(x))
  }
# cmd_set_extract_file_from_tar("data/hoge.tar") %>% cmd_get()




#' Merge some fastq files to a fastq file
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
    glue::glue("
cat {inf} > {outf} && rm {inf}
    ", .envir = as.environment(x))
  }
# cmd_set_merge_gz(c("hoge.gz", "hage.gz"), "hige.gz") %>% cmd_get()



