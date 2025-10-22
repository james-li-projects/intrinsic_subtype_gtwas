library(edgeR)
library(data.table)
library(dplyr)

# tutorial link for reference: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#######################################################
# Importing and joining vandy and GTEx count matrices #
#######################################################
# importing vanderbilt count matrix 
setwd("/gpfs/data/huo-lab/Vanderbilt")
vandy_counts <- fread("Komen.RNAseQC.counts.GENCODEv26.hg38.txt") 
vandy_counts <- data.frame(Counts)

# importing GTEx count matrix
gtex_counts <- fread("/gpfs/data/huo-lab/BCAC/james.li/GTEx/GTEx_counts/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct.gz")
gtex_counts <- data.frame(gtex_counts)
gtex_counts <- gtex_counts %>% select(-id,-Description)

# joining the two count matrices and formatting them: making rownames the ENSG IDs and the column names the sample IDs
Counts <- inner_join(vandy_counts,gtex_counts,by=c("Name"))
rownames(Counts) <- Counts$Name
Counts <- Counts %>% select(-Name,-genename)
count_mat <- data.matrix(Counts)

############################################
# Importing Joined Count Matrix into EdgeR #
############################################
# importing the count matrix into an edgeR object called DGEList
dgList <- DGEList(counts=count_mat, genes=rownames(count_mat))

#########################
# Filtering the dataset #
#########################
total_num_samples <- nrow(dgList$samples)
# filtering based on TPM>=0.1 in >=20% samples
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion >= 0.1
keep <- which(rowSums(countCheck) >= 0.2*total_num_samples)
dgList <- dgList[keep,]
# filtering based on Raw read count>=6 in >=20% samples
rawReadCounts <- dgList$counts
readCountCheck <- rawReadCounts >= 6
keep <- which(rowSums(readCountCheck) >= 0.2*total_num_samples)
dgList <- dgList[keep,]

#######################################
# Performing TMM Normalization of TPM #
#######################################
dgList <- calcNormFactors(dgList, method="TMM")
TMM_TPM_Expr <- cpm(dgList)
# outputting TMM normalized TPM 
save(TMM_TPM_Expr,file="/gpfs/data/huo-lab/Vanderbilt/james.li/TMM_TPM_Expr.RData")
