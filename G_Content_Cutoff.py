#This script is for taking in the 1000 shuffled sequences and finding the G-content value which is significant in this sample (p greater than 0.05)
#This will output a file with a column of number of Gs, the number of windows with the number of Gs, and the p value for that G count

#IMPORTS

from Bio import SeqIO #for reading the shuffled FASTA file)
import numpy as np
import sys

window = 50
step = 1
alpha = 0.05


#1 TAKING IN FASTA FILE WITH SeqIO

def read_fasta(file):
    ans = {}
    for fasta in SeqIO.parse(file, "fasta"):
        ans[fasta.id]=str(fasta.seq).upper() # so it doesnt crash
    return ans


#2 COUNTING Gs and SLIDING

def count_gs(sequence, window_size = window, step = step):

    count_g = np.frombuffer(sequence.encode(), dtype=np.uint8) == ord("G") #this counts the G in the given window
    running = np.concatenate(([0], np.cumsum(count_g))) #this makes the count cummulative instead of constant recounting
    counts = running[window_size:] - running[:-window_size]
    return counts[::step]

#this will return an array with number of Gs in the window starting at i+1 so its moving by 1nt

#3 ANALYZING G CONTENT (for finding the cutoff value)

def analyze_g_content(all_counts, alpha=alpha):

    total = len(all_counts)
    #freq counts how many windows have eadch number of G so for example "400,000 windows have 9Gs"
    freq= np.bincount(all_counts, minlength=window +1)
    #p_at_least then calculates the fraction of window with that many or more Gs which is the P value (so if 5% of windows have over or 15Gs, the p value is 0.05)
    p_at_least = freq[::-1].cumsum()[::-1] / total
    significant = [k for k in range(len(p_at_least)) if p_at_least[k] < alpha]
    #it then looks in reverse (bottom up) to find where p value stops being below 0.05 and thats the cutoff
    k = significant[0]
    return k, p_at_least, freq

#it orders it from highest p value to lowest


#4 OUTPUT FILE #output text file

def output_file(freq, p_at_least, threshold, outfile = "g_content_results.tsv"):
    with open(outfile, "w") as g:
        g.write(f"# significant G-content threshold (P<{alpha}): >= {threshold} G per {window}nt\n")
        g.write("G_count\tn_windows\tP(>=G_count)\n")
        for k in range(len(freq)):
            if freq[k] > 0:
                g.write(f"{k}\t{freq[k]}\t{p_at_least[k]:.5f}\n")
#output file is called "g_content_results.tsv"

#5 MAIN FUNCTION

def main():
    infile = sys.argv[1] if len(sys.argv) > 1 else "shuffled.fasta" #input file with any name
    outfile = "g_content_results.tsv"

    sequences = read_fasta(infile)

    g_counts = {}
    for name, seq in sequences.items():
        g_counts[name] = count_gs(seq)

    all_counts = np.concatenate(list((g_counts.values())))

    threshold, p_at_least, freq= analyze_g_content(all_counts)
    percentile_95 = np.percentile(all_counts, 95)

    output_file(freq, p_at_least, threshold, outfile)

if __name__ == "__main__":
    main()





