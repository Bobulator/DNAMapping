class SAM:

    def __init__(self, filename, sequence_name, sequence_length):
        self.filename = filename

        header = "@HD\tVN:1.0\tSO:unknown\n@SQ\tSN:%s\tLN:%s" % (sequence_name, sequence_length)
        with open(filename, 'w') as f:
            f.write(header)


    def append_sam_output(self, read_name, cigar, sequence_name, position):
        flag = 0
        mapping_quality = 255
        rnext = "*"
        pnext = 0
        tlen = 0
        quality = "*"
        read = "*"

        result = [read_name, flag, sequence_name, position, mapping_quality, cigar, rnext, pnext, tlen, read, quality]
        result_string = "\t".join([str(x) for x in result])

        with open(self.filename, 'a') as f:
            f.write("\n%s" % result_string)

        return result_string

