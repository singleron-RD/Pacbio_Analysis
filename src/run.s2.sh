sample=$1

res_path=/Personal/fuxin/dfuxin/PROJECTS/pacbio_pipe/res
homo_fa=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa
homo_gtf=/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf
sqanti3_qc=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_qc.py
sqanti3_rulerfilter=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_RulesFilter.py
collate=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py
make_seurat_input=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/make_seurat_input.py

rsp=${res_path}/${sample}
script_file=${sample}.s2.sh

echo "export PYTHONPATH=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/sequence/" > ${script_file}

#minimap
echo "minimap2 -t 30 -ax splice -uf --secondary=no -C5 ${homo_fa} ${rsp}/dedup.fasta > ${rsp}/dedup.fasta.sam 2> ${rsp}/dedup.fasta.sam.log" >> ${script_file}

#collapse
echo "sort -k 3,3 -k 4,4n ${rsp}/dedup.fasta.sam > ${rsp}/dedup.fasta.sorted.sam" >> ${script_file}
echo "collapse_isoforms_by_sam.py --input ${rsp}/dedup.fasta -s ${rsp}/dedup.fasta.sorted.sam -c 0.99 -i 0.95 --gen_mol_count -o ${rsp}/dedup.5merge" >> ${script_file}

#annotation
echo "python ${sqanti3_qc} ${rsp}/dedup.5merge.collapsed.gff ${homo_gtf} ${homo_fa} --fl_count ${rsp}/dedup.5merge.collapsed.abundance.txt -d ${rsp}/" >> ${script_file}

#filter
echo "python ${sqanti3_rulerfilter}  ${rsp}/dedup.5merge.collapsed_classification.txt ${rsp}/dedup.5merge.collapsed_corrected.fasta ${rsp}/dedup.5merge.collapsed.gff" >> ${script_file}

#error correction
echo "gzip ${rsp}/dedup.info.csv" >> ${script_file}
echo "python ${collate} ${rsp}/dedup.5merge.collapsed.group.txt ${rsp}/dedup.info.csv.gz ${rsp}/dedup.5merge.collapsed_classification.filtered_lite_classification.txt ${rsp}/dedup.annotated.csv" >> ${script_file}

#produce seurat
echo "python ${make_seurat_input} -i ${rsp}/dedup.annotated.csv -a ${rsp}/dedup.5merge.collapsed_classification.filtered_lite.gtf -o ${rsp}" >> ${script_file}