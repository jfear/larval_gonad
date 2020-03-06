from collections import namedtuple

import pysam

fname = snakemake.input[0]
oname = snakemake.output[0]

repeat = 'AAGAG' * 3
Record = namedtuple('Record', 'read_id cell_id umi seq')


def main():
    bam = pysam.AlignmentFile(fname, 'rb')
    fo = open(oname, 'w')
    fo.write('\t'.join(Record._fields) + '\n')
    for seq in bam:
        if repeat in seq.seq:
            read_id = seq.qname
            if seq.has_tag('CB'):
                cell_id = seq.get_tag('CB').rstrip('-1')
            else:
                cell_id = seq.get_tag('CR')

            if seq.has_tag('UB'):
                umi = seq.get_tag('UB')
            else:
                umi = seq.get_tag('UR')

            _seq = seq.seq
            fo.write('\t'.join(Record(read_id, cell_id, umi, _seq)) + '\n')


if __name__ == '__main__':
    main()
