"""
Primer design SCRARESCROW
Writer: Lea Schröder
Last updated: 17.04.2023
"""
#import argparse to allow more complex arguments. 
import argparse
parser = argparse.ArgumentParser(
                    prog = 'Primer design',
                    description = 'Creates a primer pair for a given gene sequence',)
#necessary argument genefile
parser.add_argument('-g', '--genefile', required=True,                                
                    help='Give a file with a nucleotide sequence in fasta format')
#optional arguments to adapt melting temperature
parser.add_argument('--min', type=int, default=55, help='Give the minimal melting temperature as an integer. Format: --min int')
parser.add_argument('--max', type=int, default=62, help='Give the maximal melting temperature as an integer. Format: --max int')
args = parser.parse_args()


#open sequence file and save content in a variable
file = open(args.genefile, "r")
content = file.read()
file.close()


#create two substrings, one containing only the header, the other one only containing the sequence (all upper case, no new line char)
header = content[0:content.find("\n")]
seq = content[content.find("\n"):].upper()
seq = seq.replace("\n", "")


#find first start codon and first stop codon (in reading frame)
start = -1
start = seq.find("ATG")
if start == -1: 
    print("No start codon could be found.")                        #error message if no start codon is found
    exit()
i = start
stop_pos = -1
while i <= len(seq):
    if seq[i:i+3] == "TAG" or seq[i:i+3] == "TGA" or seq[i:i+3] == "TAA":
        stop_pos = i
        stop_codon = seq[i:i+3]
        break
    i += 3
    
#error message if no stop codon is found    
if stop_pos == -1:
    print("No stop codon could be found.") 
    exit()
    
#error message if coding region is too short
#shortest described protein has 11 amino acids (Su et al.:Small proteins: untapped area of potential biological importance. Front Genet. 2013 Dec 16;4:286. doi: 10.3389/fgene.2013.00286.)
if stop_pos-start <= 11:
    print("The coding sequence is too short.")  
    exit()
    

#20-23 nucleotide length, 55°C < Tm < 62°C (default settings), forward and reverse primer should not diverge more than 4°C 
#formula for primer melting temp: Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC), found on https://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic

#find forward primer including the start codon
for_primers = []
for lf in range (start-20, start):                                 #primer can start max up to 20 nt before the start codon and the last possible start of the primer is the start codon itself
    for vf in range (0, 4):                                        #to get 20-23 nt long primers
          
        A = seq[lf:lf+20+vf].count("A")
        T = seq[lf:lf+20+vf].count("T")
        G = seq[lf:lf+20+vf].count("G")
        C = seq[lf:lf+20+vf].count("C")
        Tm = (64.9+41*(G+C-16.4)/(A+T+G+C))
        
        if Tm > args.min and Tm < args.max and "ATG" in seq[lf:lf+20+vf]:      #check Tm and make sure ATG is in the sequence (important for the first primers, where ATG is in the end and not included in the first ones)
            for_primers.append([seq[lf:lf+20+vf], round(Tm, 2)])   #write all possible forward primers and their melting temp in a list
        
        
#find reverse primer including the stop codon
rev_primers = []
for lr in range (stop_pos-20, stop_pos):                           #primer can start max up to 20 nt before the stop codon and the last possible start of the primer is the stop codon itself
    for vr in range (0, 4):                                        #to get 20-23 nt long primers
        
        A = seq[lr:lr+20+vr].count("A")
        T = seq[lr:lr+20+vr].count("T")
        G = seq[lr:lr+20+vr].count("G")
        C = seq[lr:lr+20+vr].count("C")
        Tm = (64.9+41*(G+C-16.4)/(A+T+G+C))
        
        if Tm > args.min and Tm < args.max and stop_codon in seq[lr:lr+20+vr]: #check Tm and make sure the above found stop codon is in the sequence (read above)
            rev_primers.append([seq[lr:lr+20+vr], round(Tm, 2)])   #write all possible reverse primers and their melting temp in a list

#error messages if no suitable primers can be found, possibility to adapt temperature settings
if len(for_primers) == 0:
    print("No forward primers could be found, you may alter the temperature settings.\nRegion around the start codon: ")
    print(seq[start-20:start]+'\033[1m'+ seq[start:start+3] + '\033[0m' + seq[start+3:start+23])
    exit()       
    
if len(rev_primers) == 0:
    print("No reverse primers could be found, you may alter the temperature settings.\nRegion around the stop codon: ")
    print(seq[stop_pos-20:stop_pos] +'\033[1m'+ seq[stop_pos:stop_pos+3] + '\033[0m' + seq[stop_pos+3:stop_pos+23])
    exit()


#interrate over both primer lists and find suitable ones, print the first pair that comes up (less or equal to 4°C difference in melting temp)
found_primer_pair = False

for r in range (len(rev_primers)):
    Tm_ref = rev_primers[r][1]
    for f in range (len(for_primers)):
        if for_primers[f][1] >= Tm_ref-4 and for_primers[f][1] <= Tm_ref+4:
            print("A suitable forward primer is " + for_primers[f][0] + ". It has a melting temperature of " + str(for_primers[f][1]) + "°C.")
            #reverse primer should be complemented reverse: switch order of string and make it opposite strand
            to_adapt = rev_primers[r][0]
            adapted_rev_p = ""
            a = len(to_adapt)-1
            while a >= 0:
                if to_adapt[a] == "A":
                    adapted_rev_p = adapted_rev_p + "T"
                elif to_adapt[a] == "T":
                    adapted_rev_p = adapted_rev_p + "A"
                elif to_adapt[a] == "C":
                    adapted_rev_p = adapted_rev_p + "G"
                elif to_adapt[a] == "G":
                    adapted_rev_p = adapted_rev_p + "C"
                a -= 1
            print("A suitable reverse primer is " + adapted_rev_p + ". It has a melting temperature of " + str(rev_primers[r][1]) + "°C.")
            found_primer_pair = True               #if no suitable primer pair is found, this remains false and an error message can be printed
            break
    break
#without breaks all suitable pairs could be printed, but it might be too much      
            
if not found_primer_pair:
    print("No suitable primer pair could be found")
