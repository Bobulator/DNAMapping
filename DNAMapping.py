import sys
from SuffixTree import SuffixTree


def map_dna(input_file):
    # Parse input
    sequence = "mississippi$"
    k = 0

    # Build Suffix Tree
    suffix_tree = SuffixTree(sequence)
    #print suffix_tree.nodes

    # Map DNA sequence
    find_dna_mapping(suffix_tree, sequence, k)

    # Output results


def find_dna_mapping(suffix_tree, sequence, k):
    # Prepare pattern matching algorithm
    suffix_array = suffix_array_from_suffix_tree(suffix_tree)
    bwt = []
    first_occurences = {}
    counts = {}

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


def suffix_array_from_suffix_tree(suffix_tree):
    suffix_array = []

    def traverse_suffix_tree(node, edge_length):
        # sort edges
        edges = sorted(node.get_edges())

        # is node a leaf node?
        if len(edges) == 0:
            suffix_array.append(node.start - edge_length)
        else:
            for edge in edges:
                child_node = suffix_tree.nodes[node.get_edge(edge)]
                traverse_suffix_tree(child_node, edge_length + (node.end - node.start))

    traverse_suffix_tree(suffix_tree.get_root_node(), 0)

    return suffix_array


if __name__ == "__main__":
    map_dna(sys.argv)
