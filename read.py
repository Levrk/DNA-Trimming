from Bio import SeqIO
import trim

for i in range (10):
    read = next(SeqIO.parse("samplefiles/sample.fastq", "fastq"))
    trim.f1(read)
