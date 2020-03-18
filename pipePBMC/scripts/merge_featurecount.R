library(readr)
library(argparse)
library(dplyr)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", nargs="*", required=TRUE, type="character", dest = "input", metavar="input.count")
parser$add_argument("--name", nargs="*", required=TRUE, type="character", dest = "name", metavar="name")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "out_data", metavar="Count.sample.tsv")
args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

file_num <- length(args$input)
if( (length(args$name)!=file_num) ) stop()

load_data <- function(f, name){
  print(c(f, name))
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  names(df) <- c("Gid", name)
  return(df)
}

for(indx in 1:file_num){
  tmp_df <- load_data(args$input[indx], args$name[indx])
  if(indx == 1){
    df <- tmp_df
  }else{
    df <- inner_join(df, tmp_df, by=c("Gid"))
  }
}

write_tsv(df, args$out_data)

