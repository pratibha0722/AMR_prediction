#import required libraries
import os
#writing present directory in variable “directory” where fasta files are given
directory=os.getcwd()
# set directory
os.chdir(directory)
print("directory")
##Function “process_fasta” to read fasta files and identify header (start with > ) replace ‘ ’ and ‘,’ with _ in header
def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = None
        for line in infile:
            line = line.strip()
            if line.startswith('>'): 
                header = line[1:]  
                header = header.replace(' ', '_')
                header = header.replace(',', '_')
                outfile.write(f">{header}\n")
            else:
                outfile.write(f"{line}\n")
# Function to process multiple files at once:
def process_multiple_files(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)
            process_fasta(input_file, output_file)
#set input folder and output folder
if __name__ == "__main__":
    input_folder = "/user/path/"  
    output_folder= "/user/path/output" 
    process_multiple_files(input_folder, output_folder)

#Step 2: Extracting sequence specific to gene families i.e transferase and hydrolase
#2.1 Transferase
#import required libraries
import os
#Function “process_fasta_and_filter” to read fasta sequences and extract those sequences that matches given string in header
def process_fasta_and_filter(input_file, output_file, gene_family):
    sequences = {}
    current_header = None
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_header = line[1:]  # Remove ">" symbol
                sequences[current_header] = ""
            else:
                sequences[current_header] += line
    # function to filter sequences based on header content
    filtered_sequences = {
        header: sequence
        for header, sequence in sequences.items()
        if gene_family.lower() in header.lower()
    }
    with open(output_file, 'w') as outfile:
        for header, sequence in filtered_sequences.items():
            outfile.write(f">{header.replace(' ', '_')}\n")
            outfile.write(f"{sequence}\n")
###Function to process multiple files at once
def process_multiple_files_and_filter(input_folder, output_folder, gene_family):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)    
    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)
            process_fasta_and_filter(input_file, output_file, gene_family)
if __name__ == "__main__":
###set input and output file folder
    input_folder = "/lustre/pratibha/copyright_ex/output"  
    output_folder = "/lustre/pratibha/copyright_ex/output2/" 
    gene_family = "transferase"  
    process_multiple_files_and_filter(input_folder, output_folder, gene_family)

#2.2 Hydrolase:
#import required libraries
import os
#Function “process_fasta_and_filter” to read fasta sequences and extract those sequences that matches given string in 
def process_fasta_and_filter(input_file, output_file, gene_families):
    sequences = {}
    current_header = None
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_header = line[1:]  # Remove ">" symbol
                sequences[current_header] = ""
            else:
                sequences[current_header] += line

    # function to Filter sequences based on header content
    filtered_sequences = {
        header: sequence
        for header, sequence in sequences.items()
        if any(gene_family.lower() in header.lower() for gene_family in gene_families)
    }
    with open(output_file, 'w') as outfile:
        for header, sequence in filtered_sequences.items():
            outfile.write(f">{header.replace(' ', '_')}\n")
            outfile.write(f"{sequence}\n")
###for multiple files at once
def process_multiple_files_and_filter(input_folder, output_folder, gene_families):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)
            process_fasta_and_filter(input_file, output_file, gene_families)
if __name__ == "__main__":
###set input and output file folder
    input_folder = "/user/path/output"  
    output_folder = "/user/path/output2a/" 
    gene_families = ["hydrol", "lactam", "peptidase"]  
    process_multiple_files_and_filter(input_folder, output_folder, gene_families)

#Step 3: Generating statistically significant k-mer features for gene families i.e transferase and hydrolase
 
#3.1 Transferase:
##import required libraries
import pandas as pd
import os
from Bio import SeqIO
# length of kmer
## set directory where list of kmers is given
directory = os.getcwd()
nmer = 9
##read file in which list of kmers is given “transferase_sig_fea.csv”
perms_as_strings_tf = pd.read_csv('/user/path/transferase_sig_fea.csv', header=None)
##function to calculate kmer frequencies using mentioned formula 
def kmer_function(mer, sequence):
    import re
    matches = re.finditer(mer, sequence, re.IGNORECASE)
    x3 = sum(1 for _ in matches)
    N = len(sequence)
    x4 = f"{x3 / (N - nmer + 1):.16f}"
    return x4
##set input output directories
input_dir = "/user/path/output2"
output_dir = "/user/path/output3/"
os.makedirs("/user/path/output3/", exist_ok=True)
os.chdir("/user/path/output2")

text = [f for f in os.listdir() if f.endswith('.fasta')]
for i in text:
    # reading fasta file
    seq_file = SeqIO.parse(i, "fasta")
    sequences = [str(record.seq) for record in seq_file]
    names = [record.id for record in SeqIO.parse(i, "fasta")]
    result = pd.DataFrame({j: [kmer_function(j, seq) for seq in sequences] for j in perms_as_strings_tf[0]})
    result.index = names
    shape = result.shape
    #print("Shape = {}".format(shape))
    x = i.split(".fasta")[0]
    result.to_csv(f"{directory}/output3/{x}tf_sig.csv")

#3.2 Hydrolase:
##import required libraries
import pandas as pd
import os
from Bio import SeqIO
import time
##set directory where list of kmers is given
directory = "/lustre/pratibha/copyright_ex"
nmer = 9
##read file in which list of kmers is given “hydrolase_sig_fea.csv”
perms_as_strings_tf = pd.read_csv('/user/path/hydrolase_sig_fea.csv', header=None)
#function to calculate kmer frequencies using mentioned formula 
def kmer_function(mer, sequence):
    import re
    matches = re.finditer(mer, sequence, re.IGNORECASE)
    x3 = sum(1 for _ in matches)
    N = len(sequence)
    x4 = f"{x3 / (N - nmer + 1):.16f}"
    return x4
##set input output directories
input_dir = "/user/path/output2a"
output_dir = "/user/path/output3a/"
os.makedirs("/user/path/output3a/", exist_ok=True)
os.chdir("/user/path/output2a")
### for all the .fasta files in input directory , read these files, applying kmer_function and write as a separate .csv file for each input fasta output directory 
text = [f for f in os.listdir() if f.endswith('.fasta')]
for i in text:
    # reading fasta file
    seq_file = SeqIO.parse(i, "fasta")
    sequences = [str(record.seq) for record in seq_file]
    names = [record.id for record in SeqIO.parse(i, "fasta")]
    result = pd.DataFrame({j: [kmer_function(j, seq) for seq in sequences] for j in perms_as_strings_tf[0]})
    result.index = names
    shape = result.shape
    #print("Shape = {}".format(shape))
    # writing the result
    x = i.split(".fasta")[0]
    result.to_csv(f"{directory}/output3a/{x}hyd_sig.csv")
	

#Step 4: Testing of WGS using developed models i.e transferase and hydrolase gene families:
#4.1 Transferase:
#Import required libraries
import pandas as pd
import os
from scipy.stats import mode
import pickle
import numpy as np
directory="/lustre/pratibha/copyright_ex"
os.chdir(directory)
## load transferase model (â€œmodel1.savâ€)
model1 =pickle.load(open('/user/path/model1.sav', 'rb'))
##set input directory where significant kmers are generated
input_directory = "/user/path/output3"
# Create a directory to store the output files
output_directory ="/user/path/output4"
os.makedirs(output_directory, exist_ok=True)

# Create a list to store results for each file
#summary_results = []
###for all the .csv files with significant feature can be processed at once
for filename in os.listdir(input_directory):
     if filename.endswith(".csv"):
         file_path = os.path.join(input_directory, filename)
         try:
            # Load the data from the current file
             data = pd.read_csv(file_path, low_memory=False)
             data.columns.values[0] = 'nameseq'
             data.set_index('nameseq', inplace=True)
             #data.set_index('nameseq', inplace=True)
             index_column = data.index
             data=np.array(data)
             predictions3 = model1.predict(data)
             results_df = pd.DataFrame({
                    'svm': predictions3,
             }, index=index_column)
             output_filename = os.path.splitext(filename)[0] + "transferase_predictions.csv"
             output_path = os.path.join(output_directory, output_filename)
             results_df.to_csv(output_path)
             print(f"Processed: {filename}")
         except Exception as e:
             print(f"Error processing {filename}: {str(e)}")

#4.2 Hydrolase:
##import required libraries
import pandas as pd
import os
from scipy.stats import mode
import pickle
import numpy as np
directory="/lustre/pratibha/copyright_ex"
os.chdir(directory)
###load hydrolase model  (â€œmodel2.savâ€)
model2=pickle.load(open('/user/path/model2.sav', 'rb'))
##set input directory where significant kmers are generated
input_directory = "/user/path/output3a"
# Create a directory to store the output files
output_directory ="/user/path/output4a/"
os.makedirs(output_directory, exist_ok=True)
# Create a list to store results for each file
#summary_results = []
##for all the .csv files with significant feature can be processed at once
for filename in os.listdir(input_directory):
     if filename.endswith(".csv"):
         file_path = os.path.join(input_directory, filename)
         try:
             data = pd.read_csv(file_path, low_memory=False)
             data.columns.values[0] = 'nameseq'
             data.set_index('nameseq', inplace=True)
             index_column = data.index
             data=np.array(data)
             predictions3 = model2.predict(data)
             results_df = pd.DataFrame({
                   'svm': predictions3,
                }, index=index_column)
             output_filename = os.path.splitext(filename)[0] + "hydrolase_predictions.csv"
             print(output_filename)
             output_path = os.path.join(output_directory, output_filename)
             print(output_path)
             results_df.to_csv(output_path)
         except Exception as e:
             print(f"Error processing {filename}: {str(e)}")

#Step 5: To extract those genes which are predicted as AMR
#5.1 Transferase:
#import required libraries
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')
###set directory where model results are stored (output4)
directory = "/lustre/pratibha/copyright_ex"
input_directory= f"{directory}/output4/"
print(input_directory)
os.makedirs(input_directory, exist_ok=True)
output_directory=f"{directory}/output5/"
os.makedirs(output_directory, exist_ok=True)
with open("AMR_genes_extraction.txt", "a") as f:
    ###Transferase
    print("Transferase_amr_gene_extraction_started", file=f)
    #tf_start_time = timeit.timeit()
    result_tf = pd.DataFrame(columns=['File Name', 'Predicted 1', 'Predicted 0', 'Name of Predicted as 1'])
### extract rows in which 1 is written in column name â€œsvmâ€ and count how many 1 and 0 are there for each file and write the final results in output5
    # for all csv in directory
    for filename in os.listdir(input_directory):
        if filename.endswith("transferase_predictions.csv"):
            file_path = os.path.join(input_directory, filename)
            # Read the CSV file into a pandas DataFrame
            df_tf = pd.read_csv(file_path)
            df_tf.set_index('nameseq', inplace=True)
           # Count occurrences of '1' and '0' in a specific column
            count_1 = (df_tf['svm'] == 1).sum()
            count_0 = (df_tf['svm'] == 0).sum()
            amr_row_names = df_tf[df_tf['svm'] == 1].index.tolist()
            result_tf= result_tf._append({'File Name': filename, 'Predicted 1': count_1, 'Predicted 0': count_0, 'Name of Predicted as 1': amr_row_names}, ignore_index=True)
            output_filename = os.path.splitext(filename)[0] + "transferase_amr_genes.csv"
            output_path = os.path.join(output_directory, output_filename)
            result_tf.to_csv(output_path)

#5.2: Hydrolase:
###import required libraries
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')
###set directory where model results are stored
directory="/lustre/pratibha/copyright_ex"
input_directory= f"{directory}/output4a/"
os.makedirs(input_directory, exist_ok=True)
output_directory=f"{directory}/output5a/"
os.makedirs(output_directory, exist_ok=True)
with open("AMR_genes_extraction.txt", "a") as f:
    result_tf = pd.DataFrame(columns=['File Name', 'Predicted 1', 'Predicted 0', 'Name of Predicted as 1'])
### extract rows in which 1 is written in column name â€œsvmâ€ and count how many 1 and 0 are there for each file and write the final results in output5
    #for all csv in directory
    for filename in os.listdir(input_directory):
        if filename.endswith("hydrolase_predictions.csv"):
            file_path = os.path.join(input_directory, filename)
            # Read the CSV file into a pandas DataFrame
            df_tf = pd.read_csv(file_path)
            df_tf.set_index('nameseq', inplace=True)
            # Count occurrences of '1' and '0' in a specific column
            count_1 = (df_tf['svm'] == 1).sum()
            count_0 = (df_tf['svm'] == 0).sum()
            # Get the row names (index) for 'a'
            amr_row_names = df_tf[df_tf['svm'] == 1].index.tolist()
            result_tf= result_tf._append({'File Name': filename, 'Predicted 1': count_1, 'Predicted 0': count_0, 'Name of Predicted as 1': amr_row_names}, ignore_index=True)
            output_filename = os.path.splitext(filename)[0] + "hydrolase_amr_genes.csv"
            output_path = os.path.join(output_directory, output_filename)
            result_tf.to_csv(output_path)

