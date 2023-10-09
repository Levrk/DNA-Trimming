from Bio import SeqIO


def main (seq):
    output = f3(f2(f1(seq)))
    printInfo(output)
    return output


def printInfo(seq):
    print("\n")
    print("ID: " + seq.id)
    print("Sequence: " + seq.seq)
    print("Description: " + seq.description)
    print("\n")
    return seq


def f1 (seq):
    #modify seq
    return seq


def f2 (seq):
    #modify seq
    return seq


def f3 (seq):
    #modify seq
    return seq
