from asyncio import subprocess
import os,sys
import json
import subprocess

out_dir = sys.argv[1]
sample = sys.argv[2]

def out_stat(c):
    with open(summary_stat,"a") as out_sta:
        out_sta.write(c+'\n')

summary_stat = out_dir+"/"+ sample + "Pacbio_qc.stat.txt"
if os.path.exists(summary_stat):
    cmd = f'rm {summary_stat}'
    subprocess.check_call(cmd, shell = True)

summary_01 = out_dir+"/01.ccs/ccs.stat.txt"
if os.path.exists(summary_01):
    with open(summary_01,"r") as s01:
        ccs_read = s01.readline().strip()
    ccs_stat  = "Read count in CCS bam file:\t"+str(ccs_read)
    out_stat("\n# 01.ccs")
    out_stat(ccs_stat)

summary_02 = out_dir+"/02.remove.primers/remove.primers.lima.summary"
if os.path.exists(summary_02):
    out_stat("\n# 02.remove.primers")
    with open(summary_02,"r") as s02:
        for line in s02:
            out_stat(line.strip())

summary_03 = out_dir+"/03.pattern.bc.umi/stat.txt"
if os.path.exists(summary_03):
    out_stat("\n# 03.pattern.bc.umi")
    with open(summary_03,"r") as s03:
        pattern_bu = s03.readline().strip()
    pbu_stat  = "Read count match BC/UMI pattern:\t"+str(pattern_bu)
    out_stat(pbu_stat)

summary_04 = out_dir+"/04.flc/fltnc.filter_summary.json"
if os.path.exists(summary_04):
    out_stat("\n# 04.flc")
    with open(summary_04,"r") as s04:
        s04_dict = json.load(s04)
        if 'attributes' in s04_dict:
            num_reads_fl = s04_dict['attributes'][1]["value"]
            num_reads_flnc = s04_dict['attributes'][2]["value"]
            s04_dict['num_reads_flnc'] = num_reads_flnc
            s04_dict['num_reads_fl'] = num_reads_fl
        stat_1_percent = round(s04_dict['num_reads_flnc']/s04_dict['num_reads_fl'],4)*100
        stat_1 = "num_reads_flnc:\t"+str(s04_dict['num_reads_flnc'])+"/"+str(s04_dict['num_reads_fl'])+"\t("+str(stat_1_percent)+"%)"
    out_stat(stat_1)

summary_05 = out_dir+"/05.split.linker/fltnc.sgr.bam.stat.txt"
if os.path.exists(summary_05):
    out_stat("\n# 05.split.linker")
    with open(summary_05,"r") as s05:
        for line in s05:
            out_stat(line.strip())

summary_06 = out_dir+"/06.dedup/dedup.bam.stat"
if os.path.exists(summary_06):
    out_stat("\n# 06.dedup")
    with open(summary_06,"r") as s06:
        read_dedup = s06.readline().strip()
    dedup_stat  = "Read count after dedup:\t"+str(read_dedup)
    out_stat(dedup_stat)