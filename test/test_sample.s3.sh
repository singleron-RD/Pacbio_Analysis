export PYTHONPATH=/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake
python /Personal/fuxin/dfuxin/PROJECTS/Pacbio/0.script/mkfq_fr_bam.py /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/06.dedup/dedup.bam  /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam
samtools sort -o /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.sort.bam /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam
samtools index -b /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.sort.bam
samtools view /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam -O SAM |  awk -F"\t" '{print ">"$1"\n"$10}' > /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa
minimap2 -t 30 -ax splice -uf --secondary=no -C5 /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa > /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.am.fa.sam 2> /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa.sam.log
celescope rna featureCounts --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --input /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa.sam --gtf_type gene --outdir /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount --assay rna --sample test_sample
python /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/mkfq_fr_bam.py /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/06.dedup/dedup.bam  /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam
samtools sort -o /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.sort.bam /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam
samtools index -b /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.sort.bam
samtools view /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam -O SAM |  awk -F"\t" '{print ">"$1"\n"$10}' > /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa
minimap2 -t 30 -ax splice -uf --secondary=no -C5 /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa > /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa.sam 2> /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa.sam.log
celescope rna featureCounts --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --input /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount/dedup.fix_id.bam.fa.sam --gtf_type gene --outdir /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/res/07.featurecount --sample test_sample
