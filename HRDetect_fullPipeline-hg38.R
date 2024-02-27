#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("ERROR: 5 argument must be supplied\n-sample_name\n-SNV_vcf_file\n-SV_bedpe_file\n-Indels_vcf_file\n-CNV_tab_file", call.=FALSE)
}

sample_name = args[1]
SNV_vcf_file = args[2]
SV_bedpe_file = args[3]
Indels_vcf_file = args[4]
CNV_tab_file = args[5]

#options(max.print=50)
library(signature.tools.lib)

## change from "tab" to "vcf"
names(SNV_vcf_file) <- sample_name
names(SV_bedpe_file) <- sample_name
names(Indels_vcf_file) <- sample_name
names(CNV_tab_file) <- sample_name

SNVcat_list <- list()
#convert to SNV catalogue
res <- vcfToSNVcatalogue(SNV_vcf_file ,genome.v = "hg38")
colnames(res$catalogue) <- sample_name
SNVcat_list[[1]] = res$catalogue

#bind the catalogues in one table
SNV_catalogues <- do.call(cbind,SNVcat_list)

#the catalogues can be plotted as follows
#plotSubsSignatures(signature_data_matrix = SNV_catalogues,
#                   plot_sum = TRUE,
#                   output_file = paste("SNV_catalogues_",sample_names[1],".jpg",sep='')
#                  )

##########################################################################################
# signatures SNV
##########################################################################################

#fit the 12 breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)
subs_fit_res <- SignatureFit_withBootstrap_Analysis(outdir = paste("signatureFit", "_", sample_name, sep=''),
                                    cat = SNV_catalogues,
                                    signature_data_matrix = COSMIC30_subs_signatures[,sigsToUse],
                                    type_of_mutations = "subs",
                                    nboot = 100,nparallel = 4)

#The signature exposures can be found here and correspond to the median
#of the boostrapped runs followed by false positive filters. See
#?SignatureFit_withBootstrap_Analysis for details
snv_exp <- subs_fit_res$E_median_filtered

#Initialise feature matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
#input_matrix <- matrix(NA,nrow = length(sample_names),
#                       ncol = length(col_hrdetect),
#                       dimnames = list(sample_names,col_hrdetect))
input_matrix <- matrix(NA,nrow = 1,
                       ncol = length(col_hrdetect),
                       dimnames = list(sample_name,col_hrdetect))

#We have already quantified the amount of SNV signatures in the samples,
#so we can supply these via the input matrix

input_matrix[colnames(snv_exp),"SNV3"] <- snv_exp["Signature.3",]
input_matrix[colnames(snv_exp),"SNV8"] <- snv_exp["Signature.8",]

## SV
sv_bedpe <- read.table(gzfile(SV_bedpe_file), sep = "\t",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE)
SVcat = bedpeToRearrCatalogue(sv_bedpe)

#run the HRDetect pipeline, for more information see ?HRDetect_pipeline
res <- HRDetect_pipeline(input_matrix,
                         genome.v = "hg38",
                         SV_catalogues = SVcat$rearr_catalogue[sample_name],
                         Indels_vcf_files = Indels_vcf_file,
                         CNV_tab_files = CNV_tab_file,
                         nparallel = 2)

#print(res)
print(res$data_matrix)
print(".")
print(input_matrix)
print(res$SV_catalogues)
print(res$exposures_rearr)
print(res$indels_classification_table)
print(res$bootstrap_fit_rearr)
print("")
print(res$hrdetect_output)
