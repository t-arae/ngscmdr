

#' Summarize FastQC summary.txt files
#' @param paths character vector which contains file paths to the FastQC summary.txt
#'
fastqc_summary <-
  function(paths){
    fpath <- status <- NULL
    paths %>%
      purrr::map(
        .f = purrr::partial(readr::read_tsv, col_names = F)
      ) %>%
      dplyr::bind_rows() %>%
      purrr::set_names(nm = c("status", "test", "fpath")) %>%
      tidyr::spread(key = fpath, value = status)
  }

#' Extract information from the FastQC data.txt file
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @param path A file path to the FastQC data.txt
#'
fastqc_data <-
  function(path){
    . <- X1 <- end <- start <- f <- fpath <- NULL

    lines <- path %>% readr::read_lines()

    mod_idx <-
      which(stringr::str_detect(lines, "^>>")) %>%
      matrix(ncol = 2, byrow = T) %>%
      tibble::as_tibble()
    colnames(mod_idx) <- c("start", "end")

    tests <-
      mod_idx$start %>%
      purrr::map(~ readr::read_tsv(path, skip = . - 1L, n_max = 1, col_names = F)) %>%
      dplyr::bind_rows() %>%
      mutate(X1 =
               stringr::str_remove(X1, ">>") %>%
               stringr::str_replace_all(" ", "_"))
    colnames(tests) <- c("test_name", "result")

    mod_idx <- dplyr::mutate(mod_idx, diff = end - start)

    fn_1 <-
      function(chr){
        chr %>%
          stringr::str_remove("[#']") %>%
          stringr::str_replace_all("[ -]", "_")
      }

    fn_2 <-
      function(fpath, start, diff){
        if(diff != 1){
          return(readr::read_tsv(file = fpath, skip = start, n_max = diff - 2))
        } else {
          return(NA)
        }
      }

    li <- list()
    for(i in seq_len(nrow(mod_idx))){
      temp_start <- mod_idx$start[i]
      temp_end <- mod_idx$end[i]
      temp_diff <- mod_idx$diff[i]

      if(i != 9){
        li[[tests$test_name[i]]] <-
          fn_2(path, temp_start, temp_diff)
      }else{
        li[[paste0(tests$test_name[i], "_1")]] <-
          readr::read_tsv(path, skip = temp_start, n_max = 1, col_names = F) %>%
          {tibble::tibble(!!.[[1]] := .[[2]])}
        li[[paste0(tests$test_name[i], "_2")]] <-
          fn_2(path, temp_start + 1, temp_diff - 1)
      }
    }

    purrr::map_if(
      li,
      .p = function(x) is.data.frame(x),
      ~ purrr::set_names(., fn_1(colnames(.)))
    )
  }

keep_by_names <-
  function(list, pattern){
    purrr::keep(list, stringr::str_detect(names(list), pattern))
  }

convert_qc_detail_li2df <-
  function(qc_detail_list){
    df_list <- list()
    for(i in unique(lapply(qc_detail_list, names) %>% unlist)){
      df_list[[i]] <-
        qc_detail_list %>%
        unlist(recursive = F) %>%
        keep_by_names(i) %>%
        purrr::keep(is.data.frame) %>%
        dplyr::bind_rows()
    }
    return(df_list)
  }
