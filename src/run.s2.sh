sample=$1
out_dir=$2

homo_fa=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa
homo_gtf=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf
sqanti3_qc=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_qc.py
sqanti3_rulerfilter=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_RulesFilter.py
collate=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py
make_seurat_input=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/make_seurat_input.py

rsp=${out_dir}/res
script_file=${out_dir}/${sample}.s2.sh

echo "export PYTHONPATH=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/sequence/" > ${script_file}

#minimap
if [ ! -d  ${rsp}/10.isoform ];then
    mkdir ${rsp}/10.isoform
fi

echo "minimap2 -t 30 -ax splice -uf --secondary=no -C5 ${homo_fa} ${rsp}/06.dedup/dedup.fasta > ${rsp}/10.isoform/dedup.fasta.sam 2> ${rsp}/10.isoform/dedup.fasta.sam.log" >> ${script_file}

#collapse
echo "sort -k 3,3 -k 4,4n ${rsp}/10.isoform/dedup.fasta.sam > ${rsp}/10.isoform/dedup.fasta.sorted.sam" >> ${script_file}
echo "collapse_isoforms_by_sam.py --input ${rsp}/06.dedup/dedup.fasta -s ${rsp}/10.isoform/dedup.fasta.sorted.sam -c 0.99 -i 0.95 --gen_mol_count -o ${rsp}/10.isoform/dedup.5merge" >> ${script_file}

#annotation
echo "python ${sqanti3_qc} ${rsp}/10.isoform/dedup.5merge.collapsed.gff ${homo_gtf} ${homo_fa} --fl_count ${rsp}/10.isoform/dedup.5merge.collapsed.abundance.txt -d ${rsp}/10.isoform" >> ${script_file}

#filter
echo "python ${sqanti3_rulerfilter}  ${rsp}/10.isoform/dedup.5merge.collapsed_classification.txt ${rsp}/10.isoform/dedup.5merge.collapsed_corrected.fasta ${rsp}/10.isoform/dedup.5merge.collapsed.gff" >> ${script_file}

if [ ! -d  ${rsp}/11.annotation ];then
    mkdir ${rsp}/11.annotation
fi
#error correction
echo "python /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/make_csv_for_dedup.py ${rsp}/06.dedup"
echo "gzip ${rsp}/06.dedup/dedup.info.csv" >> ${script_file}
echo "python ${collate} ${rsp}/10.isoform/dedup.5merge.collapsed.group.txt ${rsp}/06.dedup/dedup.info.csv.gz ${rsp}/10.isoform/dedup.5merge.collapsed_classification.filtered_lite_classification.txt ${rsp}/11.annotation/dedup.annotated.csv" >> ${script_file}

#produce seurat
echo "python ${make_seurat_input} -i ${rsp}/11.annotation/dedup.annotated.csv -a ${rsp}/10.isoform/dedup.5merge.collapsed_classification.filtered_lite.gtf -o ${rsp}/11.annotation" >> ${script_file}