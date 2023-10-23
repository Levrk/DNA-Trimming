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

def leading(seq, threshold):
    printInfo(seq)
    qualityScores = seq.letter_annotations['phred_quality']
    sequence = seq.seq

    # Create a copy of the qualityScores list
    qualityScores_copy = qualityScores.copy()

    count = 0
    for i in qualityScores:
        if i >= threshold:
            # Clear the letter_annotations before updating the sequence
            seq.letter_annotations = {}
            seq.seq = sequence[count:]
            seq.letter_annotations['phred_quality'] = qualityScores_copy[count:]
            break
        count += 1

    # Return the modified seq
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
    #modify seq
    return seq

def quality(char):
    output = ord(char) - 33
    return output

