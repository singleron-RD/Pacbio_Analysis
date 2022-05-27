#6
sample=$1
ccs_bam=$2

out_dir=$3

if [ ! -d  ${out_dir} ];then
    mkdir ${out_dir}
fi

if [ ! -d ${out_dir}/res ];then
    mkdir ${out_dir}/res
else
    rm -r ${out_dir}/res
    mkdir ${out_dir}/res
fi

res_path=${out_dir}/res
step1_script=${out_dir}/${sample}.s1.sh
primer_fasta=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/2.Detect_and_remove_primers/primers.fasta
  


#ccs
if [ ! -d ${res_path}/01.ccs ];then
    mkdir ${res_path}/01.ccs
    ln -s ${ccs_bam} ${res_path}/01.ccs/ccs.bam
fi

#Detect and remove primers
if [ ! -d ${res_path}/02.remove.primers ];then
    mkdir ${res_path}/02.remove.primers
fi
echo "lima --isoseq --dump-clips ${ccs_bam} ${primer_fasta} ${res_path}/02.remove.primers/remove.primers.bam">${step1_script}

#Detect UMI BC
if [ ! -d ${res_path}/03.pattern.bc.umi ];then
    mkdir ${res_path}/03.pattern.bc.umi
fi
echo "isoseq3 tag ${res_path}/02.remove.primers/remove.primers.5p--3p.bam ${res_path}/03.pattern.bc.umi/det_umi_bc.5p--3p.fl.bam --design T-12U-57B">>${step1_script}


#Remove PolyA and Artificial Concatemers
if [ ! -d ${res_path}/04.flc ];then
    mkdir ${res_path}/04.flc
fi
echo "isoseq3 refine ${res_path}/03.pattern.bc.umi/det_umi_bc.5p--3p.fl.bam ${primer_fasta} ${res_path}/04.flc/fltnc.bam --require-polya">>${step1_script} 


#fix barcode
extract_bc=/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/extract_bc.py
blastn=/SGRNJ/Database/script/soft/scISA-Tools/tools/blastn
parse_blastn=/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/parse_blastn.py
correct_bc=/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/correct_bc.py
bclist=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/bclist.fa

if [ ! -d ${res_path}/05.split.linker ];then
    mkdir ${res_path}/05.split.linker
fi
echo "python ${extract_bc} ${res_path}/04.flc/fltnc.bam ${res_path}/05.split.linker/align_bc_8.fa">>${step1_script}
echo "${blastn} -query ${res_path}/05.split.linker/align_bc_8.fa -db ${bclist} -outfmt 6 -word_size 6 -num_threads 8 >${res_path}/05.split.linker/whitelist_8.m6">>${step1_script}
echo "python ${parse_blastn} ${res_path}/05.split.linker/whitelist_8.m6 ${res_path}/05.split.linker/bc_correct.txt">>${step1_script}
echo "python ${correct_bc} ${res_path}/05.split.linker/bc_correct.txt ${res_path}/04.flc/fltnc.bam ${res_path}/05.split.linker/fltnc.sgr.bam">>${step1_script}

#Cluster_reads_by_UFM
if [ ! -d ${res_path}/06.dedup ];then
    mkdir ${res_path}/06.dedup
fi
echo "isoseq3 dedup ${res_path}/05.split.linker/fltnc.sgr.bam ${res_path}/06.dedup/dedup.bam --max-tag-mismatches 1 --max-tag-shift 0">>${step1_script}
