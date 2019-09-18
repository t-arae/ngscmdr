
# print.inf_info <-
#   function(x, ...){
#     cat("label: ", x$label, "; is_gz: ", x$is_gz,"; is_pe: ", x$is_pe, "\n", sep = "")
#     cat("input file: ", paste(x$inf, collapse = ", "), "\n", sep = "")
#   }

#' Set inf class
#' @param x inf class object object
set_inf_class <-
  function(x){
    class(x) <- c("inf", class(x))
    x
  }

#' Convert data.frame to list of inf class object
#' @param df a data.frame
as_input_file <-
  function(df){
    . <- NULL
    if(missing(df)){
      cat("necessary columns 'label' 'inf' 'is_gz' 'is_pe'\n")
    }else{
    if(!all(colnames(df) %in% c("label", "inf", "is_gz", "is_pe"))){
      stop("input file data.frame was wrong.")
    }
    inf_list <-
      df %>%
      split(.$label) %>%
      lapply(as.list) %>%
      lapply(function(x) lapply(x, unique)) %>%
      lapply(set_inf_class)
    names(inf_list) <- NULL
    inf_list
    }
  }

