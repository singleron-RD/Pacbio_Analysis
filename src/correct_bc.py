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
        sgr_bc = tuple(['XC',correct_bc[ident]])
        line.tags = [line.tags[0],line.tags[1],sgr_bc,line.tags[2],line.tags[4],line.tags[5],line.tags[6],
        line.tags[7],line.tags[8],line.tags[9],line.tags[10],line.tags[11],line.tags[12],line.tags[13],
        line.tags[3],line.tags[15],line.tags[16]]
        obf.write(line)
        c+=1

with open("/Personal/fuxin/dfuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/stat.txt","a") as out:
    out.write("BAM has "+str(i)+" items.")
    out.write(str(c)+" item's XC are corrected.")
    

    