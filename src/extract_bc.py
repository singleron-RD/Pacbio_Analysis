import pysam
import sys
#import pandas as pd

#bamfile = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/4.Remove_PolyA_and_Artificial_Concatemers/L210721027-L/fltnc.bam"
#out_bc8 = open("/Personal/fuxin/dfuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/align_bc_8/align_bc_8.fa","a")
#out_bc57 = open("/Personal/fuxin/dfuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/align_bc_57/align_bc_57.fa","a")

bamfile = sys.argv[1]
outfile = sys.argv[2]
out_bc8 = open(outfile,"a")

bf = pysam.AlignmentFile(bamfile,"rb",check_sq=False)

i=1
base_reverse = {"A":"T","T":"A","G":"C","C":"G"}

for line in bf:
    ident = line.get_tag('XC')
    ident_fa = '>'+str(i)+'_'+ident
    bc0 = ''.join(list(map(lambda x: base_reverse[x],[y for y in ident]))[::-1])
    bc1 = bc0[0:8]
    bc2 = bc0[24:32]
    bc3 = bc0[48:56]

    #out_bc57.write(ident_fa+'\n')
    #out_bc57.write(bc0+'\n')

    out_bc8.write(ident_fa+'_1\n')
    out_bc8.write(bc1+'\n')
    out_bc8.write(ident_fa+'_2\n')
    out_bc8.write(bc2+'\n')
    out_bc8.write(ident_fa+'_3\n')
    out_bc8.write(bc3+'\n')

    i+=1



    
