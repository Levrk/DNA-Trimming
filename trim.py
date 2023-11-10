from Bio import SeqIO
from Bio.Seq import Seq

def main (seq, count, width, threshold):
    output = leading(seq, threshold)
    output = trailing(output,threshold)
    output = count_n(output,count)
    output = trim_n(output)
    #printInfo(output)
    output = slidingWindow(seq, width, threshold)
    #printInfo(output)
    return output

def printInfo(seq):

    if (seq.seq == ""):
        return seq
    
    print("\n")
    print("ID: " + seq.id)
    print("Sequence: " + seq.seq)
    print(seq.letter_annotations['phred_quality'])
    print("\n")
    return seq

def slidingWindow (seq, width, threshold):
    # slidingWindow iterates over the sequence looking at an index 
    # of [width] bases and calculates their average quality score. 
    # Once the average quality score drops below a certain threshold indicated 
    # by [threshold] the function cuts off the remaining bases in the sequence. 
    
    if (seq.seq == ""):
        return seq

    qualityScores = seq.letter_annotations['phred_quality']
    sequence = seq.seq

    count = 0

    # Create a copy of the qualityScores list
    qualityScores_copy = qualityScores.copy()
    for i in range(len(qualityScores)):
        #determine window & window avg
        window = qualityScores[i:i+width]
        windowAVG = sum(window)/len(window)
        if (windowAVG < threshold):
            #rewrite sequence based on cutoff
            seq.letter_annotations = {}
            seq.seq = sequence[:count]
            seq.letter_annotations['phred_quality'] = qualityScores_copy[:count]
            break
        count+=1
    return seq

def leading(seq, threshold):
    # leading iterates over the sequence until it reaches a base which 
    # matches the quality score indicated by [threshold]. It then cuts 
    # off every base prior to the one which surpassed the threshold.


    if (seq.seq == ""):
        return seq
     
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

    # trailing iterates over the sequence backwards until it reaches a base which 
    # matches the quality score indicated by [threshold]. It then cuts off 
    # every base after the one which surpassed the threshold.
    
    if (seq.seq == ""):
        return seq
    
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
    if (seq.seq == ""):
            return seq

    if currCount > count: 
        seq.letter_annotations = {}
        seq.seq = ""
        return seq
    if index == -1:
        return seq
    if seq.seq[index] == "N":
        return count_n_helper(seq, count, currCount + 1, index - 1)
    else:
        return count_n_helper(seq, count, currCount, index - 1)


def count_n (seq, count):
    #count_n cuts the entire sequence if it contains more than [count] n’s 
    return count_n_helper (seq, count, 0, len(seq.seq) - 1)   

def trim_n (seq):
    if (seq.seq == ""):
        return seq
    #trim_n trims all N’s off of the beginning and end of the sequence until a non N base is hit on both sides
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

