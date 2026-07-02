# Greedy-Algorithm Superstring Assembly

## Background:
This repository contains a Python script for assembling a superstring fasta sequences from a fasta file with multiple sequences. Sequences that overlap are stringed together to build the assembled superstring.

## Script:
1. Use `parser()` to parse through the input fasta file to clean, strip, and join the header and nucleotide sequence into the `sequences` variable. 
2. `remove_contained` is used to discard sequences that are contained in larger sequences, preventing redundancy.
3. `read_overlap` is implemented to build pairwise suffix-prefix overlaps from each of the reads from the cleaned input file.
4. The `greedy_algorithm` function is used to merge the pair with the largest overlap amongst the reads.
5. This process is repeated until a uniform superstring is produced, where no further overlaps occur.

## Usage:
Run the following command in the terminal via Bash to run the script (this utilizes the current file in `sample_sequences`):

```bash
python shortestSS.py sample_sequences/dataset1_shortest_superstring_1000.fasta
```

The conda environment for reproducibility is made available for installation via the conda_environment.yaml file.

## Input:
A sample fasta file is located in the `sample_sequences` directory. This can be changed to a different input fasta fie by altering the input when running the command function from **Usage**.

## Output:
The output of the script will be placed in an initially empty directory labeled `outputs`. This can be downloaded from the repository.


## Author:

Jack Cleary