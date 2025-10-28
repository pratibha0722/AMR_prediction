##Executable version
1. The developed pipeline comprising of two scripts 1)script_step1_5.py (require python), 2)script_step6_7.R (require R), needs annotated genomes (with functional annotation, example file given here “example.fasta”) in .fasta format as input file.
2. User need to set directories whenever required (here user can replace “/user/path” with the actual directory path). The list of statistically significant kmers, built models are provided here, along with files containing antibiotic and genes information required in step 6(executable2.R). 
3. The whole script is written in two programming languages (python (3.11.2) and R (4.0.4)) and it is mentioned before every step of the required programming language. 
4. User needs to install all the required libraries (provided in the end of the document in the required packages section for both python and R). 
5. Files required: 
“example.fasta” (input file- user’s input file)
“model1.sav”, “model2.sav”, “transferase_sig_fea.csv”, “hydrolase_sig_fea.csv”, “tf_gene_name2.csv”, “hyd_ab_list.csv” (Required file to execute script)
6. Make sure all these mentioned files are in one directory, along with scripts provided here (executable1.py and executable2.R) and run these scripts. 


