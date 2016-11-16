class Parser:

    def __init__(self):
        pass

    def parse_fasta_reads(self, filename):
        results = {}

        with open(filename, 'r') as f:
            data = f.read().split()
            labels = [x[1:] for x in data[::2]]
            sequences = data[1::2]

            pairs = zip(labels, sequences)
            results = {x[0]: x[1] for x in pairs}

        return results

    def parse_fasta_sequence(self, filename):
        results = ()

        with open(filename, 'r') as f:
            data = f.read().split()
            name = data[0][1:]
            sequence = data[1]
            results = (name, sequence)

        return results