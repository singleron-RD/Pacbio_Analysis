import pysam
import sys

bc_corre = sys.argv[1]

correct_bc = {}
with open(bc_corre) as infile:
    for line in infile:
        nline = line.strip().split('\t')
        bcl = nline[0]
        if bcl not in correct_bc.keys():
            correct_bc[bcl] = nline[1]

#bamfile = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/4.Remove_PolyA_and_Artificial_Concatemers/L210721027-L/fltnc.bam"
#outbam = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/L210721027-L/fltnc.sgr.bam"

bamfile = sys.argv[2]
outbam = sys.argv[3]

bf = pysam.AlignmentFile(bamfile,"rb",check_sq=False)
obf = pysam.AlignmentFile(outbam, "wb", template = bf)

i = 0
c = 0
for line in bf:
    ident = line.get_tag('XC')
    i+=1
    if ident in correct_bc.keys():
        line.set_tag('XC',correct_bc[ident])
        obf.write(line)
        c+=1

with open(outbam +".stat.txt","a") as out:
    out.write("BAM has items:\t"+str(i)+"\n")
    stat_percent = round(c/i,4)*100
    out.write("XC corrected:\t"+str(c)+"\t("+str(stat_percent)+"%)")

    

    