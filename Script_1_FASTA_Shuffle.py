#Script 1:Take an input sequence from SARS-CoV-2 RNA Genome and make 1000 shuffled sequences in FASTA format
#Warning: I wrote this script to run on the terminal, so to run it with your input file go to the terminal and write "python Script_1.py yourfilename.fasta"
#This script takes in one sequence and shuffles it to 1000 sequences, this will be a fairly large file and can take a bit of time to finish

#IMPORTS

import random #for shuffling the list of sequences
from Bio import SeqIO #for reading the FASTA file
import sys

LW = 60

#1 TAKING IN FASTA FORMAT FILE WITH SeqIO PARSE

def read_fasta(file):
    ans = {}
    for fasta in SeqIO.parse(file, "fasta"):
        ans[fasta.id]=str(fasta.seq)
    return ans


#2 OUTPUT FILE FUNCTION

def out_fasta(info,outfyle="shuffled.fasta"):
    with open(outfyle, "w") as g:
        for name,seq in sorted (info.items()):
            g.write(">" + name + "\n")
            for i in range(0,len(seq),LW):
                g.write(seq[i:i+LW] + "\n")

#3 SHUFFLING FUNCTION

def shuffler(seq):
    letters=list(seq) #only lists can be shuffled
    random.shuffle(letters)
    return "".join(letters)

#4 MAIN FUNCTION

def main():
    num_shuffles = 1000
    input_file = sys.argv[1] if len(sys.argv)>1 else "in.fasta"
    #To use this go to the terminal and type "python Script_1.py filename.fasta", or use the direct one below:
    #input_file = filename.fasta
    out_file = "shuffled.fasta"

    sequences = read_fasta(input_file)
    shuffled = {}
    for name, seq in sequences.items():
        for n in range(1, num_shuffles + 1):
            shuffled[f"{name}_shuffle{n}"] = shuffler(seq)

    out_fasta(shuffled, out_file)




if __name__=="__main__":
    main()






