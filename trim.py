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

def trailing(seq, threshold):
    qualityScores = seq.letter_annotations['phred_quality']
    sequence = seq.seq

    # Create a copy of the qualityScores list
    qualityScores_copy = qualityScores.copy()

    count = len(qualityScores) - 1  # Start from the end of the list
    while count >= 0:
        if qualityScores[count] >= threshold:
            # Clear the letter_annotations before updating the sequence
            seq.letter_annotations = {}
            seq.seq = sequence[:count + 1]  # Include the found quality score
            seq.letter_annotations['phred_quality'] = qualityScores_copy[:count + 1]
            break
        count -= 1

    # Return the modified seq
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

