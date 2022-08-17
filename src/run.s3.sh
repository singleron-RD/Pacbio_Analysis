sample=$1
out_dir=$2

homo_fa=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa
homo_gtf=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf

mkfq_fr_bam=/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/mkfq_fr_bam.py

rsp=${out_dir}/res
script_file=${out_dir}/${sample}.s3.sh

#echo "export PYTHONPATH=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake">${script_file}

if [ ! -d  ${rsp}/07.featurecount ];then
    mkdir ${rsp}/07.featurecount
fi
echo "python ${mkfq_fr_bam} ${rsp}/06.dedup/dedup.bam  ${rsp}/07.featurecount/dedup.fix_id.bam" >>${script_file}
echo "samtools sort -o ${rsp}/07.featurecount/dedup.fix_id.sort.bam ${rsp}/07.featurecount/dedup.fix_id.bam" >>${script_file}

echo "samtools index -b ${rsp}/07.featurecount/dedup.fix_id.sort.bam" >> ${script_file}

#featurecount
echo "samtools view ${rsp}/07.featurecount/dedup.fix_id.bam -O SAM |  awk -F\"\t\" '{print \">\"\$1\"\n\"\$10}' > ${rsp}/07.featurecount/dedup.fix_id.bam.fa" >>${script_file}
echo "minimap2 -t 30 -ax splice -uf --secondary=no -C5 ${homo_fa} ${rsp}/07.featurecount/dedup.fix_id.bam.fa > ${rsp}/07.featurecount/dedup.fix_id.bam.fa.sam 2> ${rsp}/07.featurecount/dedup.fix_id.bam.fa.sam.log" >>${script_file}

echo "celescope rna featureCounts --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --input ${rsp}/07.featurecount/dedup.fix_id.bam.fa.sam --gtf_type gene --outdir ${rsp}/07.featurecount --sample ${sample}" >>${script_file}

