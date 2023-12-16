from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main (seq, count, width, thresholdLeading, thresholdTrailing, thresholdSliding, order):
    # order will be a string which contains the functions to be applied in order seperated by '+'
    orderArray = order.split('+')
    output = seq
    for i in orderArray :
        if (i == "leading"):
            output = leading(output, thresholdLeading)
        elif (i == "trailing"):
            output = trailing(output,thresholdTrailing)
        elif (i == "count"):
            output = count_n(output,count)
        elif (i == "trim"):
            output = trim_n(output)
        elif (i == "sliding"):
            output = slidingWindow(seq, width, thresholdSliding)
        elif (i == "combineTL"):
            output = combineTL(seq,thresholdLeading)
    return output

def printInfo(seq):
    #Stop right here if the sequence is empty
    if (seq.seq == ""):
        return seq
    
    print("\n")
    print("ID: " + seq.id)
    print("Sequence: " + seq.seq)
    print(seq.letter_annotations['phred_quality'])
    print("\n")
    return seq

def leading(seq, threshold):
    # leading iterates over the sequence until it reaches a base which 
    # matches the quality score indicated by [threshold]. It then cuts 
    # off every base prior to the one which surpassed the threshold.

    #Stop right here if the sequence is empty
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
    
    #Stop right here if the sequence is empty
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

def combineTL(seq, threshold):
    #Stop right here if the sequence is empty
    if (seq.seq == ""):
        return seq
    
    qualityScores = seq.letter_annotations['phred_quality']
    sequence = seq.seq

    # Create a copy of the qualityScores list
    qualityScores_copy = qualityScores.copy()

    trailingThres = False
    leadingThres = False
    
    thresholdTIndex = 0
    thresholdLIndex = 0
    
    for i in range(len(qualityScores)):
        if trailingThres and leadingThres:
            seq.letter_annotations = {}
            seq.seq = sequence[thresholdLIndex:thresholdTIndex + 1]  # Include the found quality score
            seq.letter_annotations['phred_quality'] = qualityScores_copy[thresholdLIndex:thresholdTIndex + 1]
            break
        else: 
            if qualityScores[i] >= threshold:
                leadingThres = True
                thresholdLIndex = i
            if qualityScores[-1 - i] >= threshold:
                trailingThres = True
                thresholdTIndex = -1 - i
        
    # Return the modified seq
    return seq

def count_n (seq, count):
    #count_n cuts the entire sequence if it contains more than [count] n’s 
   
    currCount = 0
    
    #Stop right here if the sequence is empty
    if (seq.seq == ""):
        return seq
    
    for i in range(len(seq.seq)):
        #print(i + currCount)
        if currCount > count:
            #print("Threshold Reached: " + i) 
            empty = Seq("")
            empty_record = SeqRecord(empty)
            empty_record.letter_annotations["phred_quality"] = []
            return empty_record
        if seq.seq[i] == "N":
            #print("N Found")
            currCount += 1
    #print("Sequence Returned: ")        
    return seq


def slidingWindow (seq, width, threshold):
    # slidingWindow iterates over the sequence looking at an index 
    # of [width] bases and calculates their average quality score. 
    # Once the average quality score drops below a certain threshold indicated 
    # by [threshold] the function cuts off the remaining bases in the sequence. 

    #Stop right here if the sequence is empty
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