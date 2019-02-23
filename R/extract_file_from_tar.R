
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
