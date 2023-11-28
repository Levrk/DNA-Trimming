from Bio import SeqIO
import sys
import trim
import time
import argparse

#set up parser
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file path')
parser.add_argument('--output', type=str, help='Output file path')
parser.add_argument('--count', type=int, help='count for trim functions')
parser.add_argument('--width', type=int, help='width for trim functions')
parser.add_argument('--thresholdL', type=int, help='threshold for leading function')
parser.add_argument('--thresholdT', type=int, help='threshold for trailing function')
parser.add_argument('--thresholdS', type=int, help='threshold for sliding function')
parser.add_argument('--order', type=str, help='order for function application')
parser.add_argument('--cutoff', type=int, help='how many reads would you like to process')
args = parser.parse_args()

# Check if 'input' and 'output' were provided
if not (args.input and args.output and args.count and args.width and args.thresholdL and args.thresholdT and args.thresholdS and args.order):
    parser.error("--input, --output, --count, --width, --thresholdL, --thresholdT, --thresholdS, and --order arguments are required.")

#set up filenames as variables and initialize count
fastq_file_in = "input/" + args.input + ".fastq"
fastq_file_out = "output/" + args.output + ".fastq"
countN = args.count
width = args.width
thresholdLeading = args.thresholdL
thresholdTrailing = args.thresholdT
thresholdSliding = args.thresholdS
order = args.order

if (args.cutoff):
    cutoff = args.cutoff
else:
    cutoff = sys.maxsize

reads = []

#count in place while we work on efficiency
count = 0

#set up time
startTime = time.time()
print("\nInitiating Sequence Trimming on ", fastq_file_in, "\n")

#iterate through each line of the file calling trim.main on each read 
for read in SeqIO.parse(fastq_file_in, "fastq"):
    
    #append the modified read to the list of reads
    newRead = trim.main(read, countN, width, thresholdLeading, thresholdTrailing, thresholdSliding, order) #need to update these defaults to make them args
    if (newRead.seq != ""):
        reads.append(newRead)
    count+=1
    if (count%500 == 0):
        print(count, "/???")
    if (count >= cutoff):
        break


#write to output .fastq file
SeqIO.write(reads, fastq_file_out, "fastq")

#calculate total time
endTime = time.time()
totalTime = endTime - startTime
print("\nTotal execution time for " + str(count) + " reads: " + str(totalTime) + " seconds")