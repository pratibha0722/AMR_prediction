#6.1 : Transferase:
directory= setwd("/user/path")
blank_df= data.frame(
  File_name= "No_data_found",
  Predicted.1="No_data_found",
  Predicted.0="No_data_found",
  AMR_gene="No_data_found",
  Resistant_antibiotics="No_data_found",
  Gene_family="No_data_found"
)

tf= read.csv('/user/path/tf_gene_name2.csv')
dir.create("/user/path/output6")
#Transferase
tf_start= Sys.time()
cat("Transferase_start_time",tf_start)
input_dir = ("/user/path/output5/")
getwd()
setwd("/user/path/output5")
pred_file=list.files(pattern = "transferase_amr_genes.csv", full.names=TRUE)
for (i in pred_file){
  tf_result= read.csv(i) 
  if(nrow(tf_result)>0){ 
    gene_tf= strsplit(tf_result$Name.of.Predicted.as.1, ",")
    gene_tf_sep= data.frame(File_name= rep(tf_result$File.Name, sapply(gene_tf, length)),
                            'Predicted 1'=rep(tf_result$Predicted.1, sapply(gene_tf, length)),
                            'Predicted 0'= rep(tf_result$Predicted.0, sapply(gene_tf, length)),
                            'AMR_gene'= unlist(gene_tf)

    )
    res_tf=gene_tf_sep
    res_tf$Resistant_antibiotics= "No resistant antibiotics predicted"
    res_tf$Gene_family="Transferase"
    for (i in 1:length(tf$gene_name)){
      #print(i)
      pos1= grep(tf$gene_name[i], gene_tf_sep$AMR_gene, ignore.case = TRUE)
      #print(pos1)
      pos2=which(length(pos1)>0)
      #print(pos2)
      if(length(pos2)==1){
        res_tf$Resistant_antibiotics[pos1]= tf$antibiotics[i]
      }
    }
  }else {
    print("no_resistant_gene_found")
    res_tf=blank_df
  }
  getwd()
  write.csv(res_tf, paste0(directory,"/output6/","transferase_antibiotics_resistance.csv", sep = ""))
}

rm(list=ls())
#Step 6.2: Hydrolase:
directory1= setwd("/user/path")
#sink("outputfile_sig.txt")
print(directory1)
start_time = Sys.time()
blank_df= data.frame(
  File_name= "No_data_found",
  Predicted.1="No_data_found",
  Predicted.0="No_data_found",
  AMR_gene="No_data_found",
  Resistant_antibiotics="No_data_found",
  Gene_family="No_data_found"
)
tf= read.csv('/user/path/hyd_ab_list.csv')
dir.create("/user/path/output6a/")
tf_start= Sys.time()
cat("Transferase_start_time",tf_start)
input_dir = ("/user/path/output5a/")
setwd("/user/path/output5a/")
getwd()
pred_file=list.files(pattern = "hydrolase_amr_genes.csv", full.names=TRUE)
for (i in pred_file){
  tf_result= read.csv(i)
  #dim(tf_result)
  if(nrow(tf_result)>0){
    #print("yes")
    gene_tf= strsplit(tf_result$Name.of.Predicted.as.1, ",")
    gene_tf_sep= data.frame(File_name= rep(tf_result$File.Name, sapply(gene_tf, length)),
                            'Predicted 1'=rep(tf_result$Predicted.1, sapply(gene_tf, length)),
                            'Predicted 0'= rep(tf_result$Predicted.0, sapply(gene_tf, length)),
                            'AMR_gene'= unlist(gene_tf)

    )
    res_tf=gene_tf_sep
    res_tf$Resistant_antibiotics= "No Resistant antibiotics predicted"
    res_tf$Gene_family="Hydrolase"
    for (i in 1:length(tf$gene_name)){
      #print(i)
      pos1= grep(tf$gene_name[i], gene_tf_sep$AMR_gene, ignore.case = TRUE)
      #print(pos1)
      pos2=which(length(pos1)>0)
      #print(pos2)
      if(length(pos2)==1){
        res_tf$Resistant_antibiotics[pos1]= tf$antibiotics[i]
      }
    }
  }else {
    print("no_resistant_gene_found")
    res_tf=blank_df
  }
  getwd()
  directory1= "/user/path"
  print(paste0(directory1,"/output6a/","hydrolase_antibiotics_resistance.csv"))
  print(directory1)
  write.csv(res_tf, paste0(directory1,"/output6a/","hydrolase_antibiotics_resistance.csv", sep = ""))}

#Step 7: Combining transferase and hydrolase results (req R programming language)
library(stringr)
directory= "/user/path"
dir.create("/user/path/output7/")
#setwd(directory)
getwd()
print(directory)
setwd(paste0(directory,"/output6/"))
pred_file=list.files(pattern = "antibiotics_resistance.csv", full.names=TRUE)
for (i in pred_file){
  tf_result= read.csv(i)
}
setwd(paste0(directory,"/output6a/"))
pred_file=list.files(pattern = "antibiotics_resistance.csv", full.names=TRUE)
for (i in pred_file){
  hyd_result= read.csv(i)
}
df= rbind(tf_result, hyd_result)
final= df[, c("File_name", "AMR_gene", "Resistant_antibiotics", "Gene_family")]
cab1=list.files(path = paste0(directory, "/user/path/"), pattern = ".fasta")
cab2= strsplit(cab1, ".fasta")
cab3= unlist(cab2)
name2=NULL
for (i in cab3){
  name2=str_c(i,"_",  name2)
}
term = c("tf_sigtransferase_predictions.csv", "hyd_sighydrolase_predictions.csv")
pattern <- paste(term, collapse = "|")
final$File_name= str_remove_all(final$File_name, pattern)
setwd(paste0(directory,"/output7/"))
write.csv(final, "final_amr_prediction_antibiotics.csv")

