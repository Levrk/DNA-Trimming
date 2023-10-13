from Bio import SeqIO

def main (seq):
    output = slidingWindow(trim_n(count_n(trailing(leading(seq)))))
    #printInfo(output)
    return output

def printInfo(seq):
    print("\n")
    print("ID: " + seq.id)
    print("Sequence: " + seq.seq)
    print("Description: " + seq.description)
    print("\n")
    return seq


def slidingWindow (seq):
    #modify seq
    return seq

def leading (seq):
    #modify seq
    return seq

def trailing (seq):
    #modify seq
    return seq

def count_n (seq):
    #modify seq
    return seq

def trim_n (seq):
    #modify seq
    return seq
