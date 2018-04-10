from Bio.Seq import Seq
import csv

# my sequence (MlaC +/- 100bp)
gene_seq = "TGACAATAAGAATAGTGGCGATGCGCCAGCTGCTGCGCCAGGTAATAATGAAACCACTGAACCTGTGGGTACAACGAAATAATTTCAGGAGAACCGACGCATGT" \
           "TTAAACGTTTAATGATGGTCGCTTTGCTGGTGATTGCACCTCTGAGTGCGGCAACCGCGGCAGACCAGACCAATCCGTATAAGCTGATGGACGAGGCGGCGCAG" \
           "AAAACGTTCGATCGCCTGAAGAATGAGCAACCGCAAATTCGGGCCAACCCGGATTATCTGCGTACCATTGTTGATCAGGAACTGCTGCCATACGTACAGGTGAA" \
           "ATACGCCGGTGCGCTGGTGCTGGGCCAGTATTACAAGAGTGCGACCCCTGCTCAACGTGAAGCCTACTTTGCCGCTTTCCGTGAGTACCTGAAGCAGGCTTACGG" \
           "TCAGGCGCTGGCGATGTATCACGGTCAAACCTATCAGATTGCGCCAGAACAGCCGCTGGGCGATAAAACCATTGTGCCTATTCGCGTTACCATTATTGACCCGAA" \
           "TGGCCGTCCGCCGGTGCGTCTGGACTTCCAGTGGCGTAAAAACTCCCAGACGGGCAATTGGCAGGCTTACGACATGATTGCTGAAGGCGTCAGTATGATCACCAC" \
           "CAAACAAAACGAGTGGGGAACGCTGCTGCGTACCAAAGGTATCGACGGCCTGACTGCGCAACTGAAATCGATTTCTCAACAGAAAATCACTCTGGAAGAGAAAAA" \
           "ATAATGAGCGAGTCACTGAGCTGGATGCAGACGGGTGACACGCTGGCGTTATCCGGAGAGCTGGATCAGGACGTTTTGCTACCGCTTTGGGAAATGCGTGAGGA"

# insert the positions of the start and stop, make sure it is -1 the actual, since python indexes from 0 onward.
start_codon_position = 100
stop_codon_position = 733

# codon number is the length of the gene (from start to stop), divided by three, excluding the stop.
codon_number = (733 - 100)/3
primers_ = {}

# function for going codon by codon, getting sequence -/+ 15bp upstream codon, mergeing with NNS and printing the seq
#for loop to go codon by codon (start to stop)
    #defines the position of current codon, each loop increases by three
    #get 15bp upstream of codon
    #get 15bp downstream of codon
    #merge upstream, NNS, and downstream sequences
    #convert string to sequence (biopython)
    #count number ATGC's
    #calculate the Tm (http://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic)
    #variable used for counting if Tm needs to be adjusted
    #if Tm is  less than 60 enter the while loop. which adds 1 base at each end to increse, until >60


def Codon_by_codon(start, end, codon_number):
    for x_ in range(0, int(codon_number)):
        start_ = start + 3*x_
        up_ = gene_seq[(start_ - 15):start_]
        down_ = gene_seq[(start_ + 3):start_ + 18]
        P = (up_ + 'NNS' + down_)
        P = Seq(P)
        y, z, w, x,  = P.count('G'), P.count('C'), P.count('A'), P.count('T')
        Tm = 64.9 + 41 * (1 + y + z - 16.4) / (3 + w + x + y + z)
        add_temp = 1
        if Tm <=60:
            while Tm <=60:
                up_ = gene_seq[(start_ - (15 + add_temp)):start_]
                down_ = gene_seq[(start_ + 3):start_ + (18 +add_temp)]
                P = (up_ + 'NNS' + down_)
                y, z, w, x, = P.count('G'), P.count('C'), P.count('A'), P.count('T')
                Tm = 64.9 + 41 * (1+ y + z - 16.4) / (3+ w + x + y + z)
                P = Seq(P)
                add_temp += 1
        if Tm >65:
            while Tm >65:
                up_ = gene_seq[(start_ - (15 - add_temp)):start_]
                down_ = gene_seq[(start_ + 3):start_ + (18 - add_temp)]
                P = (up_ + 'NNS' + down_)
                y, z, w, x, = P.count('G'), P.count('C'), P.count('A'), P.count('T')
                Tm = 64.9 + 41 * (1+ y + z - 16.4) / (3+ w + x + y + z)
                P = Seq(P)
                add_temp += 1
        P_rc = P.reverse_complement()
        y, z, w, x, = P.count('G'), P.count('C'), P.count('A'), P.count('T')
        Tm_rc = 64.9 + 41 * (1 + y + z - 16.4) / (3 + w + x + y + z)
        forward = 'F' + str(x_+1)
        reverse = 'R' + str(x_+1)
        primers_[forward] = str(P), round(Tm), len(P)
        primers_[reverse] = str(P_rc), round(Tm_rc), len(P_rc)
        # print(P, "\t",  round(Tm), "\t", P_rc, "\t", round(Tm_rc))


Codon_by_codon(start_codon_position,stop_codon_position, codon_number)
w = csv.writer(open("Primers.csv", "w"))
for key, val in primers_.items():
    w.writerow([key, val[0], val[1], val[2]])
