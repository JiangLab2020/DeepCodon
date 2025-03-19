#! /usr/lib/R/bin/Rscript --vanilla

require("tAI")		
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the path of the data file")
}
data_path <- args[1]

file_path <- paste(data_path, "tRNAscan.trna", sep = "/")
file_path2 <- paste(data_path, "output_codonM.m", sep = "/")
file_path3 <- paste(data_path, "output.tai", sep = "/")

organism.trna <- scan(file_path)
organism.ws <- get.ws(tRNA=organism.trna, sking=1)
organism.m <- matrix(scan(file_path2), ncol=61, byrow=TRUE)
organism.m <- organism.m[,-33]
organism.tai <- get.tai(organism.m, organism.ws)
write.table (organism.tai, file =file_path3, row.names = FALSE, col.names = FALSE)
