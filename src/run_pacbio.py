import os
import subprocess
import argparse


class pacbio_analysis():
    def __init__(self):
        self.sample = sample
        self.outdir = outdir
        self.ccs_bam = ccs_bam
        self.primer_fasta = primer_fasta
        self.bu_pattern = bu_pattern
        self.blastn = blastn
        self.bclist = bclist
        self.extract_bc = extract_bc
        self.parse_blastn = parse_blastn        
        self.src_correct_bc = src_correct_bc
        self.mkfq_fr_bam = mkfq_fr_bam
        self.homo_fa = homo_fa
        self.homo_gtf = homo_gtf
        self.ensembl_92 = ensembl_92
        self.genomeDir = genomeDir
        self.sqanti3_qc = sqanti3_qc
        self.sqanti3_rulerfilter = sqanti3_rulerfilter
        self.collate = collate
        self.make_seurat_input = make_seurat_input
        self.run_analysis = run_analysis
        self.make_csv_for_dedup_py = make_csv_for_dedup_py
        self.collate_FLNC_gene_info_py = collate_FLNC_gene_info_py
        self.cDNA_Cupcake_sequence_path = cDNA_Cupcake_sequence_path
        self.collapse_isoforms_by_sam_py = collapse_isoforms_by_sam_py
        self.src_summary_stat = src_summary_stat

    def ccs(self):
        out_dir = f'{self.outdir}/01.ccs'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        if not self.ccs_bam:
            exit("Error: CCS bam file is required.")
    
        self.summary_01 = f'{out_dir}/ccs.stat.txt'

        cmd = f'ln -s {self.ccs_bam} {out_dir}'
        cmd1 = f'samtools view {self.ccs_bam}|wc -l > {self.summary_01}'
        print("run...")
        print(cmd)  
        subprocess.check_call(cmd, shell =True)
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell =True)

    
    def remove_primers(self):
        out_dir = f'{self.outdir}/02.remove.primers'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        self.remove_primers_5p3p_bam = f'{out_dir}/remove.primers.5p--3p.bam'
        self.summary_02 = f'{out_dir}/remove.primers.lima.summary'

        cmd = (
            f'lima --isoseq --dump-clips {ccs_bam} {primer_fasta} {out_dir}/remove.primers.bam'
        )
        print("run...")
        print(cmd)
        subprocess.check_call(cmd, shell = True)
        
    
    def pattern_detection(self):
        out_dir = f'{self.outdir}/03.pattern.bc.umi'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        
        self.det_umi_bc_5p3p_fl_bam = f'{out_dir}/det_umi_bc.5p--3p.fl.bam'
        self.summary_03 = f'{out_dir}/stat.txt'

        cmd = (f'isoseq3 tag {self.remove_primers_5p3p_bam} {self.det_umi_bc_5p3p_fl_bam} --design {self.bu_pattern}')
        cmd1 = (f'samtools view {self.det_umi_bc_5p3p_fl_bam}|wc -l > {self.summary_03}')
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd, shell = True)
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)

    def flc(self):
        out_dir = f'{self.outdir}/04.flc'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        
        self.fltnc_bam = f'{out_dir}/fltnc.bam'
        self.summary_04 = f'{out_dir}/fltnc.filter_summary.json'

        cmd = (f'isoseq3 refine {self.det_umi_bc_5p3p_fl_bam} {primer_fasta} {self.fltnc_bam}')
        print("run...")
        print(cmd)
        subprocess.check_call(cmd, shell = True)
    
    def split_linker(self):
        out_dir = f'{self.outdir}/05.split.linker'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        self.align_bc_8_fa = f'{out_dir}/align_bc_8.fa'
        self.white_list = f'{out_dir}/whitelist_8.m6'
        self.f_bc_correct = f'{out_dir}/bc_correct.txt'
        self.fltnc_sgr_bam = f'{out_dir}/fltnc.sgr.bam'
        self.summary_05 = f'{out_dir}/fltnc.sgr.bam.stat.txt'

        cmd1 = (f'python {self.extract_bc} {self.fltnc_bam} {self.align_bc_8_fa}')
        cmd2 = (
            f'{self.blastn} -query {self.align_bc_8_fa} -db {self.bclist} -outfmt 6 -word_size 6 -num_threads 8 '
            f'>{self.white_list}'
        )
        cmd3 = (f'python {self.parse_blastn} {self.white_list} {self.f_bc_correct}')
        cmd4 = (f'python {self.src_correct_bc} {self.f_bc_correct} {self.fltnc_bam} {self.fltnc_sgr_bam}')
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("run...")
        print(cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("run...")
        print(cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("run...")
        print(cmd4)
        subprocess.check_call(cmd4, shell = True)

    def dedup(self):
        out_dir = f'{self.outdir}/06.dedup'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        
        self.dedup_bam = f'{out_dir}/dedup.bam'
        self.dedup_fa = f'{out_dir}/dedup.fasta'
        self.summary_06 = f'{out_dir}/dedup.bam.stat.txt'

        cmd = (f'isoseq3 dedup {self.fltnc_sgr_bam} {self.dedup_bam} --max-tag-mismatches 1 --max-tag-shift 0')
        cmd1 = (f'samtools view {self.dedup_bam}|wc -l > {out_dir}/dedup.bam.stat')
        print("run...")
        print(cmd)
        subprocess.check_call(cmd, shell = True)
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)
    
    def freaturecount(self):
        out_dir = f'{self.outdir}/07.featurecount'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        
        self.dedup_fixid_bam = f'{out_dir}/dedup.fix_id.bam'
        self.dedup_fixid_sort_bam = f'{out_dir}/dedup.fix_id.sort.bam'
        self.dedup_fixid_sort_bam_fa = f'{out_dir}/dedup.fix_id.bam.fa'
        self.minimap_sam = f'{out_dir}/dedup.fix_id.bam.fa.sam'
        self.fc_bam = f'{out_dir}/dedup.fix_id.bam.fa.sam.featureCounts.bam'
        self.fc_sort_bam = f'{out_dir}/dedup.fix_id.bam.fa.sam.featureCounts.sort.bam'
        self.summary_07 = f'{out_dir}/'

        cmd1 = (f'python {self.mkfq_fr_bam} {self.dedup_bam} {self.dedup_fixid_bam}')
        cmd2 = (f'samtools sort -o {self.dedup_fixid_sort_bam} {self.dedup_fixid_bam}')
        cmd3 = (f'samtools index -b {self.dedup_fixid_sort_bam}')
        awk_cmd = "awk -F\"\t\" '{print \">\"\$1\"\n\"\$10}'"
        cmd4 = (
            f'samtools view {self.dedup_fixid_sort_bam} -O -SAM|'
            f'{awk_cmd} > {self.dedup_fixid_sort_bam_fa}'
            )
        cmd5 = (f'minimap2 -t 30 -ax splice -uf --secondary=no -C5 {self.homo_fa} {self.dedup_fixid_sort_bam_fa} > {self.minimap_sam}')
        cmd6 = (
            f'celescope rna featureCounts '
            f'--genomeDir {self.ensembl_92} ' 
            f'--input {self.minimap_sam} '
            f'--gtf_type gene '
            f'--outdir {out_dir} '
            f'--sample {self.sample} '
        )
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("run...")
        print(cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("run...")
        print(cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("run...")
        print(cmd4)
        subprocess.check_call(cmd4, shell = True)
        print("run...")
        print(cmd5)
        subprocess.check_call(cmd5, shell = True)
        print("run...")
        print(cmd6)
        subprocess.check_call(cmd6, shell = True)

    def count(self):
        out_dir = f'{self.outdir}/08.count'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
       
        self.matrix_cele = f'{out_dir}/{self.sample}_matrix_10X'
        self.summary_08 = f'{out_dir}/stat.txt'

        cmd1 = (f'samtools sort -o {self.fc_bam} {self.fc_sort_bam}')
        cmd2 = (
            f'celecope rna count '
            f'--genomeDir {self.genomeDir} '
            f'--outdir {out_dir} '
            f'--sample {self.sample} '
            f'--bam {self.fc_sort_bam}'
        )
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("run...")
        print(cmd2)
        subprocess.check_call(cmd2, shell = True)
    
    def run_seurat(self):
        out_dir = f'{self.outdir}/09.seurat'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        self.summary_09 = f'{out_dir}/stat.txt'

        cmd = (
            f'Rscript {self.run_analysis} '
            f'--sample {self.sample} '
            f'--outdir {out_dir} '
            f'--matrix_file {self.matrix_cele} '
            f'--mt_gene_list None '
            f'--save_rds True'
        )
        print("run...")
        print(cmd)
        subprocess.check_call(cmd, shell = True)

    def run_isoform(self):
        out_dir = f'{self.outdir}/10.isoform'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        
        self.dedup_fa_sam = f'{out_dir}/dedup.fasta.sam'
        self.dedup_fa_sort_sam = f'{out_dir}/dedup.fasta.sorted.sam'
        self.summary_10 = ""

        cmd0 = (f'export PYTHONPATH={self.cDNA_Cupcake_sequence_path}')
        cmd1 = (
            f'minimap2 -t 30 -ax splice -uf --secondary=no -C5 {self.homo_fa} {self.dedup_fa} > {self.dedup_fa_sam}'
        )
        cmd2 = (f'sort -k 3,3 -k 4,4n {self.dedup_fa_sam} > {self.dedup_fa_sort_sam}')
        cmd3 = (
            f'{self.collapse_isoforms_by_sam_py} '
            f'--input {self.dedup_fa} '
            f'-s {self.dedup_fa_sort_sam} -c 0.99 -i 0.95 --gen_mol_count '
            f'-o {out_dir}/{self.sample}'
            )
        cmd4 = (
            f'python {self.sqanti3_qc} ' 
            f'{out_dir}/{self.sample}.collapsed.gff '
            f'{homo_gtf} {homo_fa} '
            f'--fl_count {out_dir}/{self.sample}.collapsed.abundance.txt '
            f'-d {out_dir}'
            )
        cmd5 = (
            f'python {self.sqanti3_rulerfilter} '
            f'{self.sample}.collapsed_classification.txt '
            f'{out_dir}/{self.sample}.collapsed_corrected.fasta '
            f'{out_dir}/{self.sample}.collapsed.gff'
            )
        print("run...")
        print(cmd0)
        subprocess.check_call(cmd0, shell = True)
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("run...")
        print(cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("run...")
        print(cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("run...")
        print(cmd4)
        subprocess.check_call(cmd4, shell = True)
        print("run...")
        print(cmd5)
        subprocess.check_call(cmd5, shell = True)

    def isoform_annotation(self):
        out_dir = f'{self.outdir}/11.annotation'
        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        self.dir_dedup = f'{outdir}/06.dedup'
        self.dir_isoform = f'{outdir}/10.isoform'
        self.dedup_info = f'{self.dir_dedup}/dedup.info.csv'


        cmd1 = (f'python {self.make_csv_for_dedup_py} {self.dir_dedup} ')
        cmd2 = (f'gzip {self.dedup_info}')
        cmd3 = (
            f'python {self.collate_FLNC_gene_info_py} '
            f'{self.dir_isoform}/{self.sample}.collapsed.group.txt '
            f'{out_dir}/{self.dedup_info}.gz '
            f'{self.dir_isoform}/{self.sample}.collapsed_classification.filtered_lite_classification.txt '
            f'{out_dir}/dedup.annotated.csv'
            )
        cmd4 = (
            f'python {make_seurat_input} '
            f'-i {out_dir}/dedup.annotated.csv '
            f'-a {self.dir_isoform}/{self.sample}.collapsed_classification.filtered_lite.gtf '
            f'-o {out_dir}'
            )
        print("run...")
        print(cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("run...")
        print(cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("run...")
        print(cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("run...")
        print(cmd4)
        subprocess.check_call(cmd4, shell = True)

    def summary(self):
        #QC_stat
        cmd = f'python {self.src_summary_stat} {self.outdir}'
        subprocess.check_call(cmd, shell = True)
        #celescope_stat
        #AS_report

    def run_pacbio(self):
        #self.ccs()
        self.remove_primers()
        self.pattern_detection()
        self.flc()
        self.split_linker()
        self.dedup()
        #self.freaturecount()
        #self.count()
        #self.run_seurat()
        self.run_isoform()

if __name__ == "__main__":

    sample = "test"
    outdir = "/Personal/fuxin/dfuxin/Github_repo/Pacbio_Analysis/test/test_res_1"
    ccs_bam = "/SGRNJ03/DATA03/2111/20211124_4/L210924011-L_CCS/m64236_211121_102448.hifi_reads.bam"
    primer_fasta = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/2.Detect_and_remove_primers/primers.fasta"
    bu_pattern = "T-12U-57B"
    blastn = "/SGRNJ/Database/script/soft/scISA-Tools/tools/blastn"
    bclist = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/bclist.fa"
    extract_bc = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/extract_bc.py"
    parse_blastn = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/parse_blastn.py"
    src_correct_bc = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/correct_bc.py"
    mkfq_fr_bam = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/mkfq_fr_bam.py"
    homo_fa = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa"
    homo_gtf = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf"
    sqanti3_qc = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_qc.py"
    sqanti3_rulerfilter = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_RulesFilter.py"
    collate = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py"
    make_seurat_input = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/make_seurat_input.py"
    make_csv_for_dedup_py = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/make_csv_for_dedup.py"
    ensembl_92 = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92"
    genomeDir = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92"
    run_analysis = "/SGRNJ/Public/Software/conda_env/celescope1.6.2/lib/python3.6/site-packages/celescope/tools/run_analysis.R"
    collate_FLNC_gene_info_py = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py"
    cDNA_Cupcake_sequence_path = "/Personal/fuxin/dfuxin/PROJECTS/Pacbio/cDNA_Cupcake/sequence/"
    collapse_isoforms_by_sam_py = "collapse_isoforms_by_sam.py"
    src_summary_stat = "/Personal/fuxin/dfuxin/Github_repo/Pacbio_Analysis/src/summary_stat.py"

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample')
    parser.add_argument('--outdir')
    parser.add_argument('--ccs_bam')
    parser.add_argument('--primer_fasta')
    parser.add_argument('--bu_pattern')
    parser.add_argument('--bclist')
    parser.add_argument('--genomeDir')
    args = parser.parse_args()

    sample = args.sample
    outdir = args.outdir
    ccs_bam = args.ccs_bam
    if args.primer_fasta:
        primer_fasta = args.primer_fasta
    if args.bu_pattern:
        print("Set pattern: "+args.bu_pattern)
        bu_pattern = args.bu_pattern

    PA = pacbio_analysis()
    PA.run_pacbio()

    