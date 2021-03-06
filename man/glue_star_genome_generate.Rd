% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glue-star.R
\name{glue_star_genome_generate}
\alias{glue_star_genome_generate}
\title{Generate genome index for STAR}
\usage{
glue_star_genome_generate(idx_dir = "./idx_star",
  fasta_path = "./allChr.fas",
  gff_path = "./TAIR10_GFF3_genes_transposons.gff", core_num = 4)
}
\arguments{
\item{idx_dir}{directory path for the indexing output}

\item{fasta_path}{genome fasta file path for indexing}

\item{gff_path}{.tar file path}

\item{core_num}{number of cpu core}
}
\description{
Generate genome index for STAR
}
