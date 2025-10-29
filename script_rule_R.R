##set working directory
directory= getwd()
#setwd("path/to/your/directory")
#list fasta files present in a directory, at this step you can also filter your desired files.
mt_15= read.csv("mutation_amr_data.csv") 

dir.create("/user/path/output_mut")
input_dir = ("/user/path/input/")
setwd("/user/path/input")
all_fasta2=list.files(pattern = "*.fasta")

##read the mt_15 file (mutation and wild type 15-mer unique kmer wrt to gene name and species name, this file contain other info also such as )

#load library
library(seqinr)
library(dplyr)
## code run
df_result2=NULL ##compiled table 
for (i in 1:length(all_fasta2)) # will run for all the .fasta file in the folder
{
  ##read .fasta files, 
  fasta_file= read.alignment(file = all_fasta2[i], format = "fasta", whole.header=TRUE)
  for(j in 1:length(mt_15$gene))#check each 15mer from mt_15 in given fasta file
    {
    #match gene_name in mt_15 with header of fasta sequence, this will give only those header which matched with "gene_name"
    gene2= grep(mt_15$gene[j], fasta_file$nam, ignore.case = TRUE)
    ##this loop will run if we find any match of gene_name in header
    if(length(gene2)>0)
    {
    #get the headers of matched gene only (gene2)
    fasta_file_sub_nam = fasta_file$nam[gene2]
    #get the sequence of matched sequence only (gene2)
    fasta_file_sub_seq = fasta_file$seq[gene2]
    ###header and sequence of matched gene(gene2) (subset of sequence which matched with gene_name (gene2)), and create alignment object
    fasta_file_sub = as.alignment(nb=length(fasta_file_sub_nam),nam = fasta_file_sub_nam,seq = fasta_file_sub_seq)
    ## match the 15mer corresponding to gene_name (gene2) with fasta_file_sub (sequence and header of subset )
    ### here in seq2 we will get the no. of matches present in a single gene sequence (for all sequence in a file)
    ####if match is not found it will give -1, else it gives values
    seq2=gregexpr(mt_15$nt_sub[j], fasta_file_sub$seq, ignore.case = TRUE)
    for (k in 1:length(seq2)) ###  this loop will run number of times gene_name present in header of fasta_file_sub  
    { 
      if (seq2[[k]][1] > -1) ### this loop will run only for sequence having given 15-mer
        {
        #count the number of matched 15 mers in a sequence (by counting length of match.length)
        len1= length(attr(seq2[[k]], "match.length"))
        #as we got some matches so write "found" in result variable
        result="found"
        ## bind in a df, where gene_name (which was matched in header), along with result variable, file name, and header name.
        df1 = cbind(mt_15[j, ],result, len1, all_fasta2[i], fasta_file_sub$nam[k])
        # add rows of all genes 
        df_result2 = bind_rows(df_result2,df1)
      }else
      { ## if length is zero(matched 15mer was zero), this else loop will run, write "not_found" in a result variable
        len1=0
        result="not_found"
        #bind in a df, along with gene_name for which match is not found, with result variable, file_name, header name
        df1 = cbind(mt_15[j, ],result, len1, all_fasta2[i], fasta_file_sub$nam[k])
        df_result2 = bind_rows(df_result2,df1)
      }
    }
    }
  }
}
colnames(df_result2)= c("Reported Gene", "Reported Allele", "Reported Species_name", "Product_name", "Antibiotic_class", "Antibiotic_subclass", 
                        "Genbank_protein_id", "Wildtype_AA","AA_position", "Mutation_AA","Nucleotide position", "Mutation_type", "Gene_species_id", "15mer", "Predicted label", "Functional Machinery", "Presence/absence of 15mer", "Frequency of 15mer", "User given fasta file", "User input gene in fasta file")

### here, we extract only those rows which has at least 1 match with 15-mer (using which() we got the position for the same)
df2= which(df_result2$`Frequency of 15mer` >0)
## store it in final df
df_final= df_result2[df2,]
### here df_final contains prediction for both Wild type and mutation type, if you wish to process them separately then run line 71-74, otherwise proceed from 76
mt= df_final[df_final$`Predicted label`=="MT", ]
#write.csv(mt, "mutation_prediction.csv")
user_final= mt[, c("User input gene in fasta file", "User given fasta file", "Predicted label", "Antibiotic_class", "Antibiotic_subclass", "Frequency of 15mer","15mer", "Reported Gene", "Reported Allele", "Functional Machinery")]
colnames(user_final)= c("User input gene in fasta file", "User given fasta file", "Predicted label", "Resistance antibiotic_class", "Resistance antibiotic_subclass", "Frequency of 15mer","15mer", "Reported Gene", "Reported Allele", "Functional Machinery")
write.csv(user_final, paste0(directory, "/user/path/output_mut/", "mutation_prediction.csv", sep = ""), row.names = F)


