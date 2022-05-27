sample=$1
out_dir=$2

homo_fa=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa
homo_gtf=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf
genomeDir=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92

rsp=${out_dir}/res
script_file=${out_dir}/${sample}.s4.sh

#echo "celescope rna featureCounts --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --input ${outdir}/dedup.fix_id.bam.fa.sam --gtf_type gene --outdir ${outdir}/featurecount --assay rna --sample ${sample}" >${script_file}

#featurecount_path=${rsp}/celescope/featurecount
#count_path=${rsp}/celescope/count
#seurat_path=${rsp}/celescope/seurat

echo "samtools sort -o ${rsp}/07.featurecount/dedup.fix_id.bam.fa.sam.featureCounts.sort.bam ${rsp}/07.featurecount/dedup.fix_id.bam.fa.sam.featureCounts.bam" >> ${script_file}

if [ ! -d ${rsp}/08.count ];then
    mkdir ${rsp}/08.count
fi
echo "celescope rna count --genomeDir ${genomeDir} --outdir ${rsp}/08.count --sample ${sample} --bam ${rsp}/07.featurecount/dedup.fix_id.bam.fa.sam.featureCounts.sort.bam" >> ${script_file}

if [ ! -d ${rsp}/09.seurat ];then
    mkdir ${rsp}/09.seurat
fi
echo "Rscript /SGRNJ/Public/Software/conda_env/celescope1.6.2/lib/python3.6/site-packages/celescope/tools/run_analysis.R --sample ${sample} --outdir ${rsp}/09.seurat --matrix_file ${rsp}/08.count/${sample}_matrix_10X --mt_gene_list None --save_rds True">>${script_file}
