class SAM:

    def generate_sam_output(self, read_name, sequence_name, position, read):
        flag = 0
        mapping_quality = 255
        cigar = "9M"
        rnext = "*"
        pnext = 0
        tlen = 0
        quality = "*"

        result = [read_name, flag, sequence_name, position, mapping_quality, cigar, rnext, pnext, tlen, read, quality]

        return "\t".join([str(x) for x in result])