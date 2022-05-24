sample=$1

res_path=/Personal/fuxin/dfuxin/PROJECTS/pacbio_pipe/res
homo_fa=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa
homo_gtf=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf
genomeDir=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92

rsp=${res_path}/${sample}
outdir=${rsp}/celescope
script_file=${sample}.s4.sh

echo "celescope rna featureCounts --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --input ${outdir}/dedup.fix_id.bam.fa.sam --gtf_type gene --outdir ${outdir}/featurecount --assay rna --sample ${sample}" >${script_file}

featurecount_path=${rsp}/celescope/featurecount
count_path=${rsp}/celescope/count
seurat_path=${rsp}/celescope/seurat

if [ ! -d ${count_path} ];then
    mkdir ${count_path}
fi
if [ ! -d ${seurat_path} ];then
    mkdir ${seurat_path}
fi


echo "samtools sort -o ${featurecount_path}/dedup.fix_id.bam.fa.sam.featureCounts.sort.bam ${featurecount_path}/dedup.fix_id.bam.fa.sam.featureCounts.bam" >> ${script_file}
echo "celescope rna count --genomeDir ${genomeDir} --outdir ${count_path} --assay rna --sample ${sample} --bam ${featurecount_path}/dedup.fix_id.bam.fa.sam.featureCounts.sort.bam" >> ${script_file}
echo "Rscript /SGRNJ/Public/Software/conda_env/celescope1.6.2/lib/python3.6/site-packages/celescope/tools/run_analysis.R --sample ${sample} --outdir ${seurat_path} --matrix_file ${count_path}/${sample}_matrix_10X --mt_gene_list None --save_rds True">>${script_file}
