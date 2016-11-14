import sys


def map_dna(input_file):
    # Parse input
    sequence = None
    k = 0

    # Build Suffix Tree
    suffix_tree = None

    # Map DNA sequence
    find_dna_mapping(suffix_tree, sequence, k)

    # Output results


def find_dna_mapping(suffix_tree, sequence, k):
    # Prepare pattern matching algorithm

    def find_pattern_matches(kmer):
        return []

    # Prepare list of kmers. Associate each kmer with its relative position within the sequence.
    kmer_indices = {}
    
    for i in xrange(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        
        if kmer not in kmer_indices:
            kmer_indices[kmer] = []
            
        kmer_indices[kmer].append(i)
        
    # Perform pattern matching using suffix_tree with each kmer storing results for each matching kmer index relative
    # to its location within the sequence
    candidate_mapping_indices = {}
    for kmer, relative_indices in kmer_indices.iteritems():
        matching_indices = find_pattern_matches(kmer)

        for matching_index in matching_indices:
            for relative_index in relative_indices:
                potential_sequence_index = matching_index - relative_index + 1

                if potential_sequence_index not in candidate_mapping_indices:
                    candidate_mapping_indices[potential_sequence_index] = 0

                candidate_mapping_indices[potential_sequence_index] += 1

    return candidate_mapping_indices


if __name__ == "__main__":
    map_dna(sys.argv)
