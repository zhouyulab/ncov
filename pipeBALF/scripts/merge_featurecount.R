library(readr)
library(argparse)
library(dplyr)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input.count")
parser$add_argument("--name1", nargs="*", required=TRUE, type="character", dest = "name1", metavar="name1")
parser$add_argument("--name2", nargs="*", required=TRUE, type="character", dest = "name2", metavar="name2")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "out_data", metavar="Count.sample.tsv")
args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

file_num <- length(args$input)
if( (length(args$name1)!=file_num) | (length(args$name2)!=file_num) ) stop()

load_data <- function(f, name1, name2){
  print(c(f, name1, name2))
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  names(df) <- c("Gid" ,sprintf("%s.%s", name1, name2))
  return(df)
}

for(indx in 1:file_num){
  tmp_df <- load_data(args$input[indx], args$name1[indx], args$name2[indx])
  if(indx == 1){
    df <- tmp_df
  }else{
    df <- inner_join(df, tmp_df, by=c("Gid"))
  }
}

write_tsv(df, args$out_data)

