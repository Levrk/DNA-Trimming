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

def count_n (seq, count):
    #modify seq
    print(count)
    return seq

def trim_n (seq):
    printInfo(seq)
    qualityScores = seq.letter_annotations['phred_quality']
    sequence = seq.seq
    
    # Create a copy of qualityScores list 
    qualityScores_copy = qualityScores.copy()
    
    for i in range(len(sequence)):
        if sequence[i] == "N": 
            # clear letter_annotations first
            seq.letter_annotations = {}
            seq.seq = sequence[1::]
            seq.letter_annotations['phred_quality'] = qualityScores_copy[1::]
        else:
            break
    
    for j in range(len(sequence)):
        if sequence[-j-1] == "N":
            seq.seq = sequence[0:-1]
            seq.letter_annotations['phred_quality'] = qualityScores_copy[0:-1]
        else:
            break
    # Return modified seq 
    return seq

def quality(char):
    output = ord(char) - 33
    return output

