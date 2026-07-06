# Task: Use a greedy shortest-superstring algorithm that creates a long contig by merging individual DNA reads with the strongest overlap

# 1. Import packages and libraries as necessary

import sys

path = "/Users/jackcleary/Desktop/Bioinformatics/BIOI621/dataset1_shortest_superstring_1000.fasta"

sequences = []

def parser(x): # Input the path for the fasta file as x
    sequences = []
    seq_parts = []
    header = None

    with open(x,"r") as fasta_file: # "r" = read (default)
        for line in fasta_file:
            line = line.strip().upper()
            if not line: # if line is empty after cleaning, it is skipped
                continue
            if line.startswith(">"):
                if header != None: # allows program to wait until the real sequence record is reached to save information to the sequences list
                    sequences.append((header, "".join(seq_parts)))
                header = line[1:] # slices off the < of the fasta header (obtains only the ID)
                seq_parts = []
            else:
                seq_parts.append(line) # stores the actual nucleotide sequence (will not start with <)
    if header != None:
        sequences.append((header,"".join(seq_parts)))
    return sequences


def write_fasta(filepath, header, sequence):
    line_width = 80
    with open(filepath,"w") as greedy_fasta:
        greedy_fasta.write(f">{header}\n") # adds the header and then moves to next line with \n
        for i in range(0,len(sequence), line_width):
            greedy_fasta.write(sequence[i:i+line_width]+"\n") # skips down to a new line after 80 nucleotides are written from the final consensus sequence


def read_overlap(a,b, min_len=1): # requires the overlap between read A and read B to be at least one nucleotide in length
    max_overlap = min(len(a),len(b)) # The shortest read is the threshold for the largest overlap
    for length in range(max_overlap,min_len -1, -1): # count backwards to check for longest overlap first
        if a.endswith(b[:length]): # checks if the suffix of A matches the prefix of B
            return length
    return 0


def remove_contained(sequences): # Goal is to remove reads that are inside larger reads
    seqs = list(sequences)
    seqs.sort(key=len, reverse=True) # searches for longest sequences first
    filtered_seqs = []
    for i,seq in enumerate(seqs):
        contained = False
        for j, t in enumerate(seqs):
            if i != j and len(t) >= len(seq) and seq in t:
                contained = True
                break
        if not contained:
            filtered_seqs.append(seq)
    return filtered_seqs

def overlap_matrix(filtered_seqs):
    seqs = filtered_seqs
    n = len(seqs)
    scores = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                scores[(i,j)] = read_overlap(seqs[i],seqs[j])
    return scores



def greedy_algorithm(sequences): # This will utilize the functions we defined prior in a pipeline
    seqs = list(sequences)
    filt_seqs = remove_contained(seqs)
    print(f"Number of Filtered Sequences: {len(filt_seqs)}", flush=True)

    n = len(filt_seqs)
    scores = overlap_matrix(filt_seqs)

    while len(filt_seqs) > 1:
        best_score = -1
        best_i, best_j = -1, -1
        for (i, j), s in scores.items():
            if s > best_score:
                best_score = s
                best_i, best_j = i, j

        # Merge the sequences as a string
        a = filt_seqs[best_i]
        b = filt_seqs[best_j]
        merged_seq = a + b[best_score:]

        # Remove the scores to clear the matrix
        new_idx = len(filt_seqs)
        sequences_to_remove = [k for k in scores if best_i in k or best_j in k]
        for k in sequences_to_remove:
            del scores[k]

        remaining_idx = [i for i in range(len(filt_seqs)) if i != best_i and i != best_j]

        for idx in remaining_idx:
            scores[(new_idx, idx)] = read_overlap(merged_seq, filt_seqs[idx]) # Runs the overlap function on the appended sequence prior to the rest of the sequences
            scores[(idx, new_idx)] = read_overlap(filt_seqs[idx], merged_seq)

        # Redefine the new set of sequences after merging once
        new_seqs = [filt_seqs[i] for i in remaining_idx] + [merged_seq]
        old_to_new = {old: new for new, old in enumerate(remaining_idx)}
        old_to_new[new_idx] = len(remaining_idx)

        new_scores = {}
        for (i,j), s in scores.items():
            new_i = old_to_new.get(i)
            new_j = old_to_new.get(j)
            if new_i != None and new_j != None:
                new_scores[(new_i, new_j)] = s

        filt_seqs = new_seqs
        scores = new_scores

    return filt_seqs[0]


def main(): # This will be the command used to run the script in terminal
    if len(sys.argv) < 2:
        print("Usage: python shortestSS.py <input.fasta> [output.fasta]")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2] if len(sys.argv) > 2 else "outputs/shortestSS.fasta"

    print(f"Reading sequences from: {input_fasta}")
    records = parser(input_fasta)
    print(f"Loaded {len(records)} sequences.")

    sequences = [seq for _, seq in records]

    print("Running greedy shortest superstring assembly...")
    superstring = greedy_algorithm(sequences)

    print(f"\nAssembly complete!")
    print(f"Superstring length: {len(superstring)} bp")
    print(f"Writing output to: {output_fasta}")

    write_fasta(output_fasta, "shortestSS", superstring)
    print(superstring)  # also print to stdout as required

if __name__ == "__main__":
    main()






















