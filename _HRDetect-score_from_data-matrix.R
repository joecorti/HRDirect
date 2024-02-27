library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("ERROR: 1 argument must be supplied: file_with_data-matrix", call.=FALSE)
}

input_file = args[1]
input_matrix = read.table(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
res = HRDetect_pipeline(input_matrix,
                         genome.v = "hg38")

print(res)
