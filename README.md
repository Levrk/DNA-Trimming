# DNA-Trimming
 
## main.py

Main.py iterates through each read within a .fastq file calling functions from trim.py on each and then appending the modified read to a list. After iteration is completed, the program will then write the entire list of DNA reads to a new .fastq file

This program is executed from the command line with 9 arguments specifying the names of the .fastq files that will be processed from the input folder and written to the output folder, the parameters for the trimming functions, the order for trimming functions to be applied, and the number of reads to be processed. The final argument is optional and the program will process the entire input file if no cutoff is provided. Here is an example of how you might call main.py from the command line:

> main.py --input sample --output test1 --count 3 --width 5 --thresholdL 8  --thresholdT 8  --thresholdS 4 --order leading+trailing+count+trim+sliding --cutoff 100


## trim.py

trim.py contains a number of functions used for modifying DNA reads. These functions are all called within the main function on a single argument of type SeqRecord from the biopython library. 


## The SeqRecord Object
The following descriptions outline all attributes of the SeqRecord object.

.seq – The sequence itself, typically a Seq object.

.id – The primary ID used to identify the sequence – a string. In most cases this is something like an accession number.

.name – A “common” name/id for the sequence – a string. In some cases this will be the same as the accession number, but it could also be a clone name. I think of this as being analogous to the LOCUS id in a GenBank record.

.description – A human readable description or expressive name for the sequence – a string.

.letter_annotations – Holds per-letter-annotations using a (restricted) dictionary of additional information about the letters in the sequence. The keys are the name of the information, and the information is contained in the value as a Python sequence (i.e. a list, tuple or string) with the same length as the sequence itself. This is often used for quality scores (e.g. Section 20.1.6) or secondary structure information (e.g. from Stockholm/PFAM alignment files).
32

.annotations – A dictionary of additional information about the sequence. The keys are the name of the information, and the information is contained in the value. This allows the addition of more “unstructured” information to the sequence.

.features – A list of SeqFeature objects with more structured information about the features on a sequence (e.g. position of genes on a genome, or domains on a protein sequence). The structure of sequence features is described below in Section 4.3.

.dbxrefs - A list of database cross-references as strings.

[source](http://biopython.org/DIST/docs/tutorial/Tutorial-1.81.pdf)
