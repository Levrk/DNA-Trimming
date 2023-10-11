from Bio import SeqIO
import trim
import time
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
reads = []
#count in place while we work on efficiency
count = 0

#set up time
startTime = time.time()

#iterate through each line of the file calling trim.main on each read 
for read in SeqIO.parse(fastq_file_in, "fastq"):
    #append the modified read to the list of reads
    reads.append(trim.main(read))
    count+=1
    #if count > 10 :
        #break


#write to output .fastq file
SeqIO.write(reads, fastq_file_out, "fastq")

#calculate total time
endTime = time.time()
totalTime = endTime - startTime
print("Total execution time for " + str(count) + " reads: " + str(totalTime) + " seconds")