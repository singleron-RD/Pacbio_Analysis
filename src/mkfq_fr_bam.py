import pysam
import sys

bamfile = sys.argv[1]
outfile = sys.argv[2]

bf = pysam.AlignmentFile(bamfile,"rb",check_sq=False)
obf = pysam.AlignmentFile(outfile, "wb", template = bf)

c = 0

for line in bf:
    correct_name = line.get_tag('XC')+'_'+line.get_tag('XM')+'_'+line.qname
    line.qname = correct_name
    obf.write(line)
    c += 1

print("Corrected items: "+str(c))
