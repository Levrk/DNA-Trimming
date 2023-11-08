from Bio import SeqIO
import trim
import time
import argparse

#set up parser
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file path')
parser.add_argument('--output', type=str, help='Output file path')
parser.add_argument('--count', type=int, help='count for trim functions')
parser.add_argument('--width', type=int, help='width for trim functions')
parser.add_argument('--threshold', type=int, help='threshold for trim functions')
parser.add_argument('--cutoff', type=int, help='how many reads would you like to process')
args = parser.parse_args()

# Check if 'input' and 'output' were provided
if not (args.input and args.output and args.count and args.width and args.threshold):
    parser.error("--input, --output, --count, --width, and --threshold arguments are required.")

#set up filenames as variables and initialize count
fastq_file_in = "input/" + args.input + ".fastq"
fastq_file_out = "output/" + args.output + ".fastq"
count = args.count
width = args.width
threshold = args.threshold

if (args.cutoff):
    cutoff = args.cutoff
else:
    cutoff = 10

reads = []

#count in place while we work on efficiency
count = 0

#set up time
startTime = time.time()

#iterate through each line of the file calling trim.main on each read 
for read in SeqIO.parse(fastq_file_in, "fastq"):
    #append the modified read to the list of reads
    reads.append(trim.main(read, count, width, threshold)) #need to update these defaults to make them args
    count+=1
    if (count >= cutoff):
        break


#write to output .fastq file
SeqIO.write(reads, fastq_file_out, "fastq")

#calculate total time
endTime = time.time()
totalTime = endTime - startTime
print("Total execution time for " + str(count) + " reads: " + str(totalTime) + " seconds")