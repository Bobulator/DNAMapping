class Parser:

    def __init__(self):
        pass

    def parse_fasta_reads(self, filename):
        results = {}

        with open(filename, 'r') as f:
            for line in f:
                if line[0] == '>':
                    read_name = line[1:].strip()
                    results[read_name] = ''

                else:
                    results[read_name] += line.strip().upper()

        return results

    def parse_fasta_sequence(self, filename):
        results = ()

        with open(filename, 'r') as f:
            data = f.read().split()
            name = data[0][1:]
            sequence = "".join(data[1:]).upper()
            results = (name, sequence)

        return results