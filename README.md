# MatePEreads

This script will go over two fastq files simultaneously and extract true pairs, while saving unpaired reads into a separate file.

    python MatePEreads.py <base_file_name>

## Running the code

Lets say you have your `R1` and `R2` file:

    $ ls
    seqs_R1.fastq.gz
    seqs_R2.fastq.gz

Run the code like this:

    python MatePEreads.py seqs

The program will print some messages, and you will have at the end these files:

    $ ls
    seqs_R1.fastq.gz
    seqs_R2.fastq.gz
    seqs_paired_R1.fastq.gz
    seqs_paired_R2.fastq.gz
    seqs_unpaired.fastq.gz

The suffix `_R1/R2.fastq.gz` is necessary unless you want to tweak the code to match your needs.

Currently, it will only take gzipped or `gz` fastq files as input, and will only save gzipped files as output.

The code also uses `gzcat` or `zcat` internally, so change that line (62) to match your system.