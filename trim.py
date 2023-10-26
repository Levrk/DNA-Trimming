from Bio import SeqIO
from Bio.Seq import Seq

def main (seq, count, width, threshold):
    output = leading(seq, threshold)
    output = trailing(output,threshold)
    output = count_n(output,count)
    output = trim_n(output)
    output = slidingWindow(output, width, threshold)
    printInfo(output)
    return output

def printInfo(seq):
    print("\n")
    print("ID: " + seq.id)
    print("Sequence: " + seq.seq)
    print(seq.letter_annotations['phred_quality'])
    print("\n")
    return seq

def slidingWindow (seq, width, threshold):
    #modify seq
    return seq

def leading (seq, threshold):
    print(threshold)
    return seq

def trailing (seq, threshold):
    #modify seq
    print(threshold)
    return seq

def count_n_helper (seq, count, currCount, index):
    if currCount > count:
        seq.seq = Seq("")
        return seq
    if index == -1:
        return seq
    if seq.seq[index] == "N":
        return count_n_helper(seq, count, currCount + 1, index - 1)
    else:
        return count_n_helper(seq, count, currCount, index - 1)


def count_n (seq, count):
     return count_n_helper (seq, count, 0, len(seq.seq) - 1)   

def trim_n (seq):
    #modify seq
    return seq

def quality(char):
    output = ord(char) - 33
    return output

