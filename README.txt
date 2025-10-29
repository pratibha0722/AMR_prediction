##Executable version
For AMR prediction, we have integrated two approach: 1) Machine leaning based model 2) Rule based alogorithm.
To execute these scripts, user need to set directories whenever required (here user need to set directry where the fasta files are stored, user  can replace “/user/path” with the actual directory path). For ML based prediction, the list of statistically significant kmers, built models are provided here, along with files containing antibiotic and genes information required in step 6(script6_7_R.R). For rule based algorithm, script_ruleR.R and mutation_amr_data.csv is provided. 

For ML based prediction:
    1. The developed pipeline comprising of two scripts 1)script_step1_5.py (require python), 2)script_step6_7.R (require R), needs annotated genomes (with functional annotation, example file given here “example.fasta”) in .fasta format as input file.
    3. The whole script is written in two programming languages (python (3.11.2) and R (4.0.4)) and it is mentioned before every step of the required programming language. 
    4. User needs to install all the required libraries (provided in the end of the document in the required packages section for both python and R). 
    5. Files required: 
    “example.fasta” (input file- user’s input file (.fasta format))
    “model1.sav”, “model2.sav”, “transferase_sig_fea.csv”, “hydrolase_sig_fea.csv”, “tf_gene_name2.csv”, “hyd_ab_list.csv” , "tf_gene_name2.csv", "hyd_ab_list.csv" (Required file to execute script)
    6. Make sure all these mentioned files are in one directory, along with scripts provided here (script_step1_5_python.py and script_step6_7_R.R) and run these scripts. 

For Rule based prediction:
      1. The required files to execute the script (script_ruleR.R) is i) mutation_amr_data.csv (given) ii) fasta file (user input/ "example file.fasta") are provided.
      2. User need to set directory where the fasta files are placed. User needs to replace "/user/path" with the actual directory path. 
      3. The prediction using Rule based algorithm will be saved in directory "output_mut.csv"


### The result using ML based model can be find in "/user/path/output7/ml_amr_prediction_antibiotics.csv" and rule based preediction can be find in "/user/path/output_mut/mutation_prediction.csv" 


