import sys

in_file = "bclist"
out_file = "bclist.fa"

with open(out_file,"a") as out:
    with open(in_file) as infile:
        i=1
        for line in infile:
            out.write('>'+str(i)+'\n')
            out.write(line)
            i+=1

