class Parser:

    def __init__(self):
        pass

    def parse_fasta_reads(self, filename):
        results = {}

        with open(filename, 'r') as f:
            data = f.read().split()
            results = {data[i][1:]: data[i+1] for i in range(0, len(data) - 1, 2)}

        return results

    def parse_fasta_sequence(self, filename):
        results = ()

        with open(filename, 'r') as f:
            data = f.read().split()
            name = data[0][1:]
            sequence = "".join(data[1:])
            results = (name, sequence)

        return results
