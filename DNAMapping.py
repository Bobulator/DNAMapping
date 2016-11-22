import sys
from SuffixTree import SuffixTree
from Parser import Parser
from SAM import SAM


def map_dna(args):
    # Parse input
    parser = Parser()
    sequence_name, sequence = parser.parse_fasta_sequence(args[1])
    read_file = args[2]
    k = int(args[3])

    # Build Suffix Tree
    print("Building Suffix Tree...",)
    suffix_tree = SuffixTree(sequence)
    print("DONE")

    print("Building Suffix Array...",)
    suffix_array = suffix_array_from_suffix_tree(suffix_tree)
    print("DONE")

    del suffix_tree

    print("Building BWT...",)
    bwt = bwt_from_suffix_array(suffix_array, sequence)
    print("DONE")

    print("Building First Occurrence Table...")
    first_occurrences = build_first_occurrence(bwt)
    print("DONE")

    print("Building Counts Table...")
    counts = build_counts(bwt)
    print("DONE")

    read_file_name = read_file.split(".")[0]
    sam = SAM(filename="%s.SAM" % sequence_name, sequence_length=len(sequence))

    reads = parser.parse_fasta_reads(read_file)

    # Map DNA sequence
    print("Beginning Mapping Process")
    counter = 0
    for read_name, read in reads.iteritems():
        counter += 1

        print("Mapping Read " + str(counter) + " ")

        candidate_index = find_dna_mapping(suffix_array, bwt, first_occurrences, counts, read, k)

        print(sam.append_sam_output(read_name=read_name, sequence_name=sequence_name,
                                    position=candidate_index+1, read=read))


def find_dna_mapping(suffix_array, bwt, first_occurrences, counts, read, k):

    def find_pattern_matches(kmer):
        indices = []

        top = 0
        bottom = len(bwt) - 1
        kmer_index = len(kmer) - 1

        while top <= bottom:
            # print top, bottom
            if kmer_index >= 0:
                symbol = kmer[kmer_index]
                kmer_index -= 1

                if symbol in bwt[top:bottom + 1]:
                    top = first_occurrences[symbol] + counts[symbol][top]
                    bottom = first_occurrences[symbol] + counts[symbol][bottom + 1] - 1
                else:
                    break
            else:
                for i in xrange(top, bottom + 1):
                    indices.append(suffix_array[i])
                break

        return indices

    # Prepare list of kmers. Map each kmer with its position(s) within the read.
    kmer_indices = {}
    
    for i in xrange(len(read) - k + 1):
        kmer = read[i:i + k]
        
        if kmer not in kmer_indices:
            kmer_indices[kmer] = []
            
        kmer_indices[kmer].append(i)
        
    # Perform pattern matching using suffix_tree with each kmer storing results for each matching kmer index relative
    # to its location within the read, i.e. matching_index - relative_index = potential_read_mapping_index
    candidate_mapping_indices = {}
    max_index = -1
    max_score = -1
    for kmer, relative_indices in kmer_indices.iteritems():
        matching_indices = find_pattern_matches(kmer)

        for matching_index in matching_indices:
            for relative_index in relative_indices:
                potential_read_index = matching_index - relative_index

                if potential_read_index not in candidate_mapping_indices:
                    candidate_mapping_indices[potential_read_index] = 0

                candidate_mapping_indices[potential_read_index] += 1

                if candidate_mapping_indices[potential_read_index] > max_score:
                    max_index = potential_read_index
                    max_score = candidate_mapping_indices[potential_read_index]

    return max_index


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


def bwt_from_suffix_array(suffix_array, sequence):
    bwt = []

    for index in suffix_array:
        if index == 0:
            bwt.append(sequence[len(suffix_array) - 1])
        else:
            bwt.append(sequence[index - 1])

    return bwt


def build_first_occurrence(bwt):
    first_occurrence = {}
    sorted_bwt = sorted(bwt)

    for i in xrange(len(sorted_bwt)):
        if sorted_bwt[i] not in first_occurrence:
            first_occurrence[sorted_bwt[i]] = i

    return first_occurrence


def build_counts(bwt):
    counts = {}

    for i in xrange(len(bwt)):
        if bwt[i] not in counts:
            counts[bwt[i]] = [0] * (i + 1)

        for key in counts.keys():
            if key == bwt[i]:
                counts[key].append(counts[key][-1] + 1)
            else:
                counts[key].append(counts[key][-1])

    return counts


if __name__ == "__main__":
    map_dna(sys.argv)
