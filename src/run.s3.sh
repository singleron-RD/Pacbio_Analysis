sample=$1

res_path=/Personal/fuxin/dfuxin/PROJECTS/pacbio_pipe/res
homo_fa=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa
homo_gtf=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf

mkfq_fr_bam=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/0.script/mkfq_fr_bam.py

rsp=${res_path}/${sample}
outdir=${rsp}/celescope
script_file=${sample}.s3.sh

echo "export PYTHONPATH=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake">${script_file}

echo "python ${mkfq_fr_bam} ${rsp}/dedup.bam  ${outdir}/dedup.fix_id.bam" >>${script_file}
echo "samtools sort -o ${outdir}/dedup.fix_id.sort.bam ${outdir}/dedup.fix_id.bam" >>${script_file}

echo "samtools index -b ${outdir}/dedup.fix_id.sort.bam" >> ${script_file}

#featurecount
echo "samtools view ${outdir}/dedup.fix_id.bam -O SAM |  awk -F\"\t\" '{print \">\"\$1\"\n\"\$10}' > ${outdir}/dedup.fix_id.bam.fa" >>${script_file}
echo "minimap2 -t 30 -ax splice -uf --secondary=no -C5 ${homo_fa} ${outdir}/dedup.fix_id.bam.fa > ${outdir}/dedup.fix_id.am.fa.sam 2> ${outdir}/dedup.fix_id.bam.fa.sam.log" >>${script_file}
#echo "celescope rna featureCounts --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --input ${outdir}/dedup.fix_id.bam.fa.sam --gtf_type gene --outdir ${outdir} --assay rna --sample ${sample}" >>${script_file}

