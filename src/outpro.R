## last update: 2020/03/28

args <- commandArgs(T)
job_name = args[1]
tmp = args[2]
outdir = args[3]
timestamp = args[4]
library(data.table)
load(paste(tmp,"/original.RData", sep = ""))
file = paste(tmp,"/scIGANs-", job_name,".csv", sep = "")

d <- fread(file, header = T)
# new add
# gcm 行是细胞，列是基因
# 进行反归一化
# gcm <- t(apply(d, 1, function(x) {
#   original_data <- x * (reads_max_cell - reads_min_cell) + reads_min_cell
#   return(original_data)
# }))

# gcm <- d[,-1]*reads_max_cell
gcm<-d[,-1]
gcm_out <- cbind(Gene_ID = genenames, t(gcm[,1:geneCount]))
colnames(gcm_out) = cellnames

outfile = paste(outdir,"/scIGANs_",job_name,".txt", sep = "")
fwrite(as.data.table(gcm_out), file=outfile, col.names = T, row.names = F, sep = "\t", quote = F)
message(paste("\nCompleted!!!", Sys.time()))
message(paste("Imputed matrix:", outfile))

