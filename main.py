from Bio import SeqIO
import trim
import argparse

#set up parser
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file path')
parser.add_argument('--output', type=str, help='Output file path')
args = parser.parse_args()

# Check if 'input' and 'output' were provided
if not (args.input and args.output):
    parser.error("Both --input and --output arguments are required.")

#set up filenames as variables and initialize count
fastq_file_in = "input/" + args.input + ".fastq"
fastq_file_out = "output/" + args.output + ".fastq"
count = 0


for read in SeqIO.parse(fastq_file_in, "fastq"):
    trim.main(read)
    count+=1
    if count > 10 :
        break
