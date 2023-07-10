import sys

bclist_file = sys.argv[1]
in_file = sys.argv[2]
out_file = sys.argv[3]

i = 1
bclist = {}
with open(bclist_file) as infile:
    line = infile.readline()
    while line:
        if line.startswith('>'):
            id = line.strip()[1:]
            seq = infile.readline().strip()
            bclist[id] = seq
        line = infile.readline()

blast = {}
with open(in_file) as infile:
    for line in infile:
        nline = line.strip().split('\t')
        if int(nline[3])>7:
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

with open(out_file,"w") as out:
    for item in blast:
        if len(blast[item].keys()) == 4:
            newline = blast[item]['ident']+'\t'+blast[item]['1']+blast[item]['2']+blast[item]['3']
            out.write(newline+'\n')