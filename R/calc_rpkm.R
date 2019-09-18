
#' Calculate Reads Per Million mapped reads (RPM)
#' @param readcount read count. atomic vector
#' @export
calc_rpm <- function(readcount){
  (readcount / sum(readcount, na.rm = T)) * 10^6
}

#' Calculate Reads Per Killobase of exon per Million mappled reads (RPKM)
#' @param readcount read count. atomic vector
#' @param len feature length. atomic vector
#' @export
calc_rpkm <- function(readcount, len){
  readcount * (10^3 / len) * (10^6 / sum(readcount, na.rm = T))
}

#' Calculate TPM
#' @param readcount read count. atomic vector
#' @param len feature length. atomic vector
#' @export
calc_tpm <- function(readcount, len){
  t <- readcount / len * 10^3
  t / sum(t, na.rm = T) * 10^6
}

#' Calculate RPM from featureCounts output
#' @importFrom dplyr mutate
#' @importFrom stringr str_detect
#' @param df data.frame output from merge_featurecount_output()
#' @export
calc_rpm_from_featurecounts <-
  function(df){
    cn <- colnames(df)
    cn <- cn[7:length(cn)]
    cn <- cn[!str_detect(cn, "^(rpm_|rpkm_|tpm_)")]

    for(i in cn){
      df <- mutate(df, hoge = calc_rpm(df[[i]]))
      colnames(df)[colnames(df) == "hoge"] <- paste0("rpm_", i)
    }
    df
  }

#' Calculate RPKM from featureCounts output
#' @importFrom dplyr mutate
#' @importFrom stringr str_detect
#' @param df data.frame output from merge_featurecount_output()
#' @export
calc_rpkm_from_featurecounts <-
  function(df){
    cn <- colnames(df)
    cn <- cn[7:length(cn)]
    cn <- cn[!str_detect(cn, "^(rpm_|rpkm_|tpm_)")]

    for(i in cn){
      df <- mutate(df, hoge = calc_rpkm(df[[i]], df[["Length"]]))
      colnames(df)[colnames(df) == "hoge"] <- paste0("rpkm_", i)
    }
    df
  }

#' Calculate TPM from featureCounts output
#' @importFrom dplyr mutate
#' @importFrom stringr str_detect
#' @param df data.frame output from merge_featurecount_output()
#' @export
calc_tpm_from_featurecounts <-
  function(df){
    cn <- colnames(df)
    cn <- cn[7:length(cn)]
    cn <- cn[!str_detect(cn, "^(rpm_|rpkm_|tpm_)")]

    for(i in cn){
      df <- mutate(df, hoge = calc_tpm(df[[i]], df[["Length"]]))
      colnames(df)[colnames(df) == "hoge"] <- paste0("tpm_", i)
    }
    df
  }
