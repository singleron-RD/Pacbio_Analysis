import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

i = 1
bclist = {}
with open(in_file) as infile:
    for line in infile:
        bclist[str(i)] = line.strip()
        i+=1

blast = {}
with open(in_file) as infile:
    for line in infile:
        nline = line.strip().split('\t')
        if int(nline[3])>6:
            nfa = nline[0].split('_')
            seqid = nfa[0]
            ident = nfa[1]
            bcid = nfa[2]
            if seqid in blast.keys():
                if bcid in blast[seqid].keys():
                    pass
                else:
                    blast[seqid][bcid] = bclist[nline[1]]
            else:
                blast[seqid] = {"ident":ident,bcid:bclist[nline[1]]}

with open(out_file,"a") as out:
    for item in blast:
        if len(blast[item].keys()) == 4:
            newline = blast[item]['ident']+'\t'+blast[item]['1']+blast[item]['2']+blast[item]['3']
            out.write(newline+'\n')