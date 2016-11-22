class SAM:

    def __init__(self, filename, sequence_length):
        self.filename = filename

        header = "@HD\tVN:1.0\tSO:unknown\n@SQ\tSN:%s\tLN:%s" % (filename.split(".")[0], sequence_length)
        with open(filename, 'w') as f:
            f.write(header)


    def append_sam_output(self, read_name, sequence_name, position, read):
        flag = 0
        mapping_quality = 255
        cigar = "%dM" % len(read)
        rnext = "*"
        pnext = 0
        tlen = 0
        quality = "*"

        result = [read_name, flag, sequence_name, position, mapping_quality, cigar, rnext, pnext, tlen, read, quality]
        result_string = "\t".join([str(x) for x in result])

        with open(self.filename, 'a') as f:
            f.write("\n%s" % result_string)

        return result_string