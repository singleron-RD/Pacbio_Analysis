import os
import subprocess
import argparse


class pacbio_analysis():
    def __init__(self):
        #args
        self.sample = sample
        self.outdir = outdir
        self.ccs_bam = ccs_bam
        self.primer_fasta = primer_fasta
        self.bu_pattern = bu_pattern
        self.blastn = blastn
        self.bclist = bclist
        self.genomeDir = genomeDir
        self.barcode_match = barcode_match
        self.steps= step_list
        self.report = report

        #config
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
        self.featurecount_bam_py = featurecount_bam_py
        self.make_gene_seurat = make_gene_seurat
        self.hgnc_gene_set = hgnc_gene_set

    def check_outdir(self):
        if not os.path.exists(self.outdir):
            cmd_dir = (f'mkdir ${self.outdir}')

    def ccs(self):
        out_dir = f'{self.outdir}/01.ccs'

        self.summary_01 = f'{out_dir}/ccs.stat.txt'
        self.ccs_bam_ls = out_dir+"/"+self.ccs_bam.split('/')[-1]

        if os.path.exists(self.ccs_bam_ls):
            return('CCS bam file already exist!')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)
        if not self.ccs_bam:
            exit("Error: CCS bam file is required.")

        cmd = f'ln -s {self.ccs_bam} {out_dir}'
        cmd1 = f'samtools view {self.ccs_bam}|wc -l > {self.summary_01}'
        print("Running command: ",cmd) 
        subprocess.check_call(cmd, shell =True)
        print("Running command: ",cmd1) 
        subprocess.check_call(cmd1, shell =True)
        return('Get ccs bam file!')
    
    def remove_primers(self):
        out_dir = f'{self.outdir}/02.remove.primers'

        self.remove_primers_5p3p_bam = f'{out_dir}/remove.primers.5p--3p.bam'
        self.summary_02 = f'{out_dir}/remove.primers.lima.summary'

        self.ccs_bam_ls = self.outdir+"/01.ccs/"+self.ccs_bam.split('/')[-1]
        print(self.ccs_bam_ls)
        if not os.path.exists(self.ccs_bam_ls):
            return('No ccs file found! Consider add --steps ccs')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd = (
            f'lima --isoseq --dump-clips {ccs_bam} {primer_fasta} {out_dir}/remove.primers.bam'
        )
        print("Running command: ",cmd)
        subprocess.check_call(cmd, shell = True)
        return('Remove primers succeed!')
        
    
    def pattern_detection(self):
        out_dir = f'{self.outdir}/03.pattern.bc.umi'

        self.remove_primers_5p3p_bam = f'{self.outdir}/02.remove.primers/remove.primers.5p--3p.bam'
        self.det_umi_bc_5p3p_fl_bam = f'{out_dir}/det_umi_bc.5p--3p.fl.bam'
        self.summary_03 = f'{out_dir}/stat.txt'

        if not os.path.exists(self.remove_primers_5p3p_bam):
            return('Remove primer step should be run before this, Consider add --steps remove_primers')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)    

        cmd = (f'isoseq3 tag {self.remove_primers_5p3p_bam} {self.det_umi_bc_5p3p_fl_bam} --design {self.bu_pattern}')
        cmd1 = (f'samtools view {self.det_umi_bc_5p3p_fl_bam}|wc -l > {self.summary_03}')
        print("Running command: ",cmd)
        subprocess.check_call(cmd, shell = True)
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        return('Pattern detection succeed!')

    def flc(self):
        out_dir = f'{self.outdir}/04.flc'

        self.det_umi_bc_5p3p_fl_bam = f'{self.outdir}/03.pattern.bc.umi/det_umi_bc.5p--3p.fl.bam'
        self.fltnc_bam = f'{out_dir}/fltnc.bam'
        self.summary_04 = f'{out_dir}/fltnc.filter_summary.json'

        if not os.path.exists(self.det_umi_bc_5p3p_fl_bam):
            return('Pattern detection step should be run before this, Consider add --steps pattern_detection')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd = (f'isoseq3 refine {self.det_umi_bc_5p3p_fl_bam} {primer_fasta} {self.fltnc_bam}')
        print("Running command: ",cmd)
        subprocess.check_call(cmd, shell = True)
        return('Flc step succeed!')
    
    def split_linker(self):
        out_dir = f'{self.outdir}/05.split.linker'

        self.fltnc_bam = f'{self.outdir}/04.flc/fltnc.bam'
        self.align_bc_8_fa = f'{out_dir}/align_bc_8.fa'
        self.white_list = f'{out_dir}/whitelist_8.m6'
        self.f_bc_correct = f'{out_dir}/bc_correct.txt'
        self.fltnc_sgr_bam = f'{out_dir}/fltnc.sgr.bam'
        self.summary_05 = f'{out_dir}/fltnc.sgr.bam.stat.txt'

        if not os.path.exists(self.fltnc_bam):
            return('Flc step should be run before this, Consider add --steps flc')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd1 = (f'python {self.extract_bc} {self.fltnc_bam} {self.align_bc_8_fa}')
        cmd2 = (
            f'{self.blastn} -query {self.align_bc_8_fa} -db {self.bclist} -outfmt 6 -word_size 6 -num_threads 8 '
            f'>{self.white_list}'
        )
        cmd3 = (f'python {self.parse_blastn} {self.bclist} {self.white_list} {self.f_bc_correct}')
        cmd4 = (f'python {self.src_correct_bc} {self.f_bc_correct} {self.fltnc_bam} {self.fltnc_sgr_bam}')
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("Running command: ",cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("Running command: ",cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("Running command: ",cmd4)
        subprocess.check_call(cmd4, shell = True)
        return('Split linker succeed!')

    def dedup(self):
        out_dir = f'{self.outdir}/06.dedup'

        self.fltnc_sgr_bam = f'{self.outdir}/05.split.linker/fltnc.sgr.bam'
        self.dedup_bam = f'{out_dir}/dedup.bam'
        self.dedup_fa = f'{out_dir}/dedup.fasta'
        self.summary_06 = f'{out_dir}/dedup.bam.stat.txt'

        if not os.path.exists(self.fltnc_sgr_bam):
            return('Split linker step should be run before this, Consider add --steps split_linker')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd = (f'isoseq3 dedup {self.fltnc_sgr_bam} {self.dedup_bam} --max-tag-mismatches 1 --max-tag-shift 0')
        cmd1 = (f'samtools view {self.dedup_bam}|wc -l > {out_dir}/dedup.bam.stat')
        print("Running command: ",cmd)
        subprocess.check_call(cmd, shell = True)
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        return('Dedup succeed!')
    
    def featurecount(self):
        out_dir = f'{self.outdir}/07.featurecount'

        self.dedup_bam = f'{self.outdir}/06.dedup/dedup.bam'
        self.dedup_fixid_bam = f'{out_dir}/dedup.fix_id.bam'
        self.dedup_fixid_sort_bam = f'{out_dir}/dedup.fix_id.sort.bam'
        self.dedup_fixid_sort_bam_fa = f'{out_dir}/dedup.fix_id.bam.fa'
        self.minimap_sam = f'{out_dir}/dedup.fix_id.bam.fa.sam'
        self.fc_bam = f'{out_dir}/dedup.fix_id.bam.fa.sam.featureCounts.bam'
        self.summary_07 = f'{out_dir}/'

        if not os.path.exists(self.dedup_bam):
            return('Dedup step should be run before this, Consider add --steps dedup')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd1 = (f'python {self.mkfq_fr_bam} {self.dedup_bam} {self.dedup_fixid_bam}')
        cmd2 = (f'samtools sort -o {self.dedup_fixid_sort_bam} {self.dedup_fixid_bam}')
        cmd3 = (f'samtools index -b {self.dedup_fixid_sort_bam}')
        awk_cmd = "awk -F\"\\t\" '{print \">\"$1\"\\n\"$10}'"
        cmd4 = (
            f'samtools view {self.dedup_fixid_sort_bam} -O SAM|'
            f'{awk_cmd} > {self.dedup_fixid_sort_bam_fa}'
            )
        cmd5 = (f'minimap2 -t 30 -ax splice -uf --secondary=no -C5 {self.homo_fa} {self.dedup_fixid_sort_bam_fa} > {self.minimap_sam}')
        #cmd6 = (
        #    f'celescope rna featureCounts '
        #    f'--genomeDir {self.ensembl_92} ' 
        #    f'--input {self.minimap_sam} '
        #    f'--gtf_type gene '
        #    f'--outdir {out_dir} '
        #    f'--sample {self.sample} '
        #)
        cmd6 = (
            'featureCounts '
            '-s 1 '
            f'-a {self.homo_gtf} '
            f'-o {out_dir}/{self.sample} '  # not bam
            '-R CORE '
            f'-T 1 '
            f'-t gene '
            f'{self.minimap_sam} '
            f'-L'
        )
        cmd7 = (
            f'python {self.featurecount_bam_py} ' 
            f'{self.homo_gtf} '
            f'{self.minimap_sam}.featureCounts '
            f'{self.minimap_sam}' 
        )
        cmd8 = (
            f'samtools view -bS '
            f'{self.minimap_sam}.featureCounts.temp_sam > '
            f'{self.minimap_sam}.featureCounts.bam'
        )
        cmd9 = (
            f'samtools sort '
            f'{self.minimap_sam}.featureCounts.bam > '
            f'{self.minimap_sam}.featureCounts.sort.bam'
        )
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("Running command: ",cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("Running command: ",cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("Running command: ",cmd4)
        subprocess.check_call(cmd4, shell = True)
        print("Running command: ",cmd5)
        subprocess.check_call(cmd5, shell = True)
        print("Running command: ",cmd6)
        subprocess.check_call(cmd6, shell = True)
        print("Running command: ",cmd7)
        subprocess.check_call(cmd7, shell = True)
        print("Running command: ",cmd8)
        subprocess.check_call(cmd8, shell = True)
        print("Running command: ",cmd9)
        subprocess.check_call(cmd9, shell = True)

        return('Featurecount succeed!')

    def count(self):
        out_dir = f'{self.outdir}/08.count'

        self.fc_bam = f'{self.outdir}/07.featurecount/dedup.fix_id.bam.fa.sam.featureCounts.bam'
        self.fc_sort_bam = f'{self.outdir}/07.featurecount/dedup.fix_id.bam.fa.sam.featureCounts.sort.bam'
        self.matrix_cele = f'{out_dir}/{self.sample}_matrix_10X'
        self.summary_08 = f'{out_dir}/stat.txt'

        if not os.path.exists(self.fc_bam):
            return('Featurecount step should be run before this, Consider add --steps featurecount')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd1 = (f'samtools sort -o {self.fc_bam} {self.fc_sort_bam}')
        cmd2 = (
            f'celescope rna count '
            f'--genomeDir {self.genomeDir} '
            f'--outdir {out_dir} '
            f'--sample {self.sample} '
            f'--bam {self.fc_sort_bam}'
        )
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("Running command: ",cmd2)
        subprocess.check_call(cmd2, shell = True)
        return('Count succeed!')
    
    def run_seurat(self):
        out_dir = f'{self.outdir}/09.seurat'

        self.matrix_cele =f'{self.outdir}/08.count/{self.sample}_filtered_feature_bc_matrix'
        self.summary_09 = f'{out_dir}/stat.txt'

        if not os.path.exists(self.matrix_cele):
            return('Count step should be run before this, Consider add --steps count')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        cmd = (
            f'Rscript {self.run_analysis} '
            f'--sample {self.sample} '
            f'--outdir {out_dir} '
            f'--matrix_file {self.matrix_cele} '
            f'--mt_gene_list None '
            f'--save_rds True'
        )
        print("Running command: ",cmd)
        subprocess.check_call(cmd, shell = True)
        return('Run seurat succeed!')

    def run_isoform(self):
        out_dir = f'{self.outdir}/10.isoform'

        self.dedup_fa = f'{self.outdir}/06.dedup/dedup.fasta'
        self.dedup_fa_sam = f'{out_dir}/dedup.fasta.sam'
        self.dedup_fa_sort_sam = f'{out_dir}/dedup.fasta.sorted.sam'
        self.summary_10 = ""

        if not os.path.exists(self.dedup_fa):
            return('Dedup step should be run before this, Consider add --steps dedup')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

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
            f'{out_dir}/{self.sample}.collapsed_classification.txt '
            f'{out_dir}/{self.sample}.collapsed_corrected.fasta '
            f'{out_dir}/{self.sample}.collapsed.gff'
            )
        print("Set PYTHONPATH...")
        os.environ['PYTHONPATH'] = self.cDNA_Cupcake_sequence_path
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("Running command: ",cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("Running command: ",cmd3)
        subprocess.check_call(cmd3, shell = True)
        print("Running command: ",cmd4)
        subprocess.check_call(cmd4, shell = True)
        print("Running command: ",cmd5)
        subprocess.check_call(cmd5, shell = True)
        return('Run isoform succeed!')

    def isoform_annotation(self):
        out_dir = f'{self.outdir}/11.annotation'

        self.dedup_fa = f'{self.outdir}/06.dedup/dedup.fasta'
        self.dedup_bam = f'{self.outdir}/06.dedup/dedup.bam'
        self.dir_dedup = f'{self.outdir}/06.dedup'
        self.dir_isoform = f'{self.outdir}/10.isoform'
        self.dedup_info = f'{self.dir_dedup}/dedup.info.csv'

        if not os.path.exists(self.dedup_fa):
            return('Dedup step should be run before this, Consider add --steps dedup')

        if not os.path.exists(self.dedup_bam):
            return('Dedup step should be run before this, Consider add --steps dedup')

        if not os.path.exists(out_dir):
            cmd_dir = (f'mkdir {out_dir}')
            subprocess.check_call(cmd_dir, shell = True)

        if not os.path.exists(f'{out_dir}/known_isoform_seurat'):
            cmd_dir1 = (f'rm -r {out_dir}/known_isoform_seurat')
            subprocess.check_call(cmd_dir1, shell = True)
        if not os.path.exists(f'{out_dir}/with_novel_isoform_seurat'):
            cmd_dir2 = (f'rm -r {out_dir}/with_novel_isoform_seurat')
            subprocess.check_call(cmd_dir2, shell = True)
        if not os.path.exists(f'{out_dir}/gene_seurat'):
            cmd_dir3 = (f'rm -r {out_dir}/gene_seurat')
            subprocess.check_call(cmd_dir3, shell = True)

        cmd1 = (f'python {self.make_csv_for_dedup_py} {self.dir_dedup}')
        cmd2 = (f'gzip {self.dedup_info}')
        cmd3 = (
            f'python {self.collate_FLNC_gene_info_py} '
            f'{self.dir_isoform}/{self.sample}.collapsed.group.txt '
            f'{self.dedup_info}.gz '
            f'{self.dir_isoform}/{self.sample}.collapsed_classification.filtered_lite_classification.txt '
            f'{out_dir}/dedup.annotated.csv'
            )
        cmd4 = (
            f'python {make_seurat_input} '
            f'-i {out_dir}/dedup.annotated.csv '
            f'-a {self.dir_isoform}/{self.sample}.collapsed_classification.filtered_lite.gtf '
            f'-o {out_dir} '
            )
        cmd5 = (
            f'mv {out_dir}/isoforms_seurat {out_dir}/known_isoform_seurat'
        )
        cmd6 = (
            f'python {make_seurat_input} '
            f'-i {out_dir}/dedup.annotated.csv '
            f'-a {self.dir_isoform}/{self.sample}.collapsed_classification.filtered_lite.gtf '
            f'-o {out_dir} '
            f'--keep_novel'
        )
        cmd7 = (
            f'mv {out_dir}/isoforms_seurat {out_dir}/with_novel_isoform_seurat'
        )
        cmd8 = (
            f'mkdir {out_dir}/gene_seurat'
        )
        cmd9 = (
            f'Rscript {self.make_gene_seurat} '
            f'--known_isoform_seurat {out_dir}/known_isoform_seurat '
            f'--outdir {out_dir}/gene_seurat '
            f'--gene_name_set {self.hgnc_gene_set}'
        )
        print("Running command: ",cmd1)
        subprocess.check_call(cmd1, shell = True)
        print("Running command: ",cmd2)
        subprocess.check_call(cmd2, shell = True)
        print("Running command: ",cmd3)
        if os.path.exists(f'{out_dir}/dedup.annotated.csv'):
            print('Annotation file already exists')
        else:
            subprocess.check_call(cmd3, shell = True)
        print("Running command: ",cmd4)
        print("Running command: ",cmd5)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)
        print("Running command: ",cmd6)
        print("Running command: ",cmd7)
        subprocess.check_call(cmd6, shell = True)
        subprocess.check_call(cmd7, shell = True)
        print("Running command: ", cmd8)
        print("Running command: ", cmd9)
        subprocess.check_call(cmd8, shell = True)
        subprocess.check_call(cmd9, shell = True)
        return('Isoform annotation succeed!')

    def summary(self):
        #QC_stat
        cmd = f'python {self.src_summary_stat} {self.outdir} {self.sample}'
        print("Running command: ",cmd)
        subprocess.check_call(cmd, shell = True)
        #celescope_stat
        #AS_report
        return("QC summary done!")

    def run_pacbio(self):
        if 'ccs' in self.steps:
            rec = self.ccs()
            print(rec)

        if 'remove_primer' in self.steps:
            rec = self.remove_primers()
            print(rec)
        if 'pattern_detection' in self.steps:
            rec = self.pattern_detection()
            print(rec)
        if 'flc' in self.steps:
            rec = self.flc()
            print(rec)
        if 'split_linker' in self.steps:
            rec = self.split_linker()
            print(rec)
        if 'dedup' in self.steps:
            rec = self.dedup()
            print(rec)
        if 'run_isoform' in self.steps:
            rec = self.run_isoform()
            print(rec)
            rec = self.isoform_annotation()
            print(rec)
            cmd_report = f'cp {self.outdir}/10.isoform/{self.sample}.collapsed_classification.filtered_lite_SQANTI3_report.html {self.outdir}/{self.sample}.SQANTI3_report.html'
            subprocess.check_call(cmd_report, shell = True)

        if 'featurecount' in self.steps:
            rec = self.featurecount()
            print(rec)
        if 'count' in self.steps:
            rec = self.count()
            print(rec)
        if 'run_seurat' in self.steps:
            rec = self.run_seurat()
            print(rec)
        if self.report == "True":
            rec = self.summary()    
            print(rec)

if __name__ == "__main__":
    cw_dir = os.getcwd()
    print("Work directory:"+cw_dir)
    sample = "sample"
    #outdir = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/test_res_1"
    #ccs_bam = "/SGRNJ03/DATA03/2111/20211124_4/L210924011-L_CCS/m64236_211121_102448.hifi_reads.bam"
    primer_fasta = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/2.Detect_and_remove_primers/primers.fasta"
    bu_pattern = "T-12U-57B"
    blastn = "/SGRNJ/Database/script/soft/scISA-Tools/tools/blastn"
    bclist = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/5.Split_linker_and_barcode/bclist.fa"
    extract_bc = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/extract_bc.py"
    parse_blastn = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/parse_blastn.py"
    src_correct_bc = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/correct_bc.py"
    mkfq_fr_bam = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/mkfq_fr_bam.py"
    homo_fa = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.fa"
    homo_gtf = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf"
    sqanti3_qc = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_qc.py"
    sqanti3_rulerfilter = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/9.Compare_Against_Annotation/SQANTI3-4.2/sqanti3_RulesFilter.py"
    collate = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py"
    make_seurat_input = "/Personal/fuxin/dfuxin/Github_repo/dev/Pacbio_Analysis/src/make_seurat_input_update_mtx.py"
    make_csv_for_dedup_py = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/make_csv_for_dedup.py"
    ensembl_92 = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92"
    genomeDir = "/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92"
    run_analysis = "/SGRNJ/Public/Software/conda_env/celescope1.6.2/lib/python3.6/site-packages/celescope/tools/run_analysis.R"
    collate_FLNC_gene_info_py = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py"
    cDNA_Cupcake_sequence_path = "/SGRNJ03/randd/user/fuxin/PROJECTS/Pacbio/cDNA_Cupcake/sequence/"
    collapse_isoforms_by_sam_py = "collapse_isoforms_by_sam.py"
    src_summary_stat = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/summary_stat.py"
    featurecount_bam_py = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/featurecount_bam.py"
    make_gene_seurat = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/src/make_gene_seurat.R"
    hgnc_gene_set = "/SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/config/hgnc_complete_set.txt"

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample')
    parser.add_argument('--outdir')
    parser.add_argument('--ccs_bam')
    parser.add_argument('--primer_fasta')
    parser.add_argument('--bu_pattern')
    parser.add_argument('--bclist')
    parser.add_argument('--genomeDir')
    parser.add_argument('--barcode_match',default = "None")
    parser.add_argument('--steps')
    parser.add_argument('--report',default = "True")
    parser.add_argument('--mapfile',default = "None")
    parser.add_argument('--pacbio_source_path',default = "None")
    args = parser.parse_args()

    sample = args.sample
    outdir = args.outdir
    ccs_bam = args.ccs_bam
    report = args.report
    
    if args.primer_fasta:
        primer_fasta = args.primer_fasta
    if args.pacbio_source_path:
        pacbio_source_path = args.pacbio_source_path
        extract_bc = pacbio_source_path + "/src/extract_bc.py"
        parse_blastn = pacbio_source_path + "/src/parse_blastn.py"
        src_correct_bc = pacbio_source_path + "/src/correct_bc.py"
        mkfq_fr_bam = pacbio_source_path + "/src/mkfq_fr_bam.py"
        sqanti3_qc = pacbio_source_path + "/tools/SQANTI3-4.2/sqanti3_qc.py"
        sqanti3_rulerfilter = pacbio_source_path + "/tools/SQANTI3-4.2/sqanti3_RulesFilter.py"
        collate = pacbio_source_path + "/tools/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py"
        make_seurat_input = pacbio_source_path + "/src/make_seurat_input_update_mtx.py"
        make_csv_for_dedup_py = pacbio_source_path + "/src/make_csv_for_dedup.py"
        collate_FLNC_gene_info_py = pacbio_source_path + "/tools/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py"
        cDNA_Cupcake_sequence_path = pacbio_source_path + "/tools/cDNA_Cupcake/sequence/"
        src_summary_stat = pacbio_source_path + "/src/summary_stat.py"
        featurecount_bam_py = pacbio_source_path + "/src/featurecount_bam.py"
        make_gene_seurat = pacbio_source_path + "/src/make_gene_seurat.R"
        hgnc_gene_set = pacbio_source_path + "/config/hgnc_complete_set.txt"


    if args.bu_pattern:
        print("Set pattern: "+args.bu_pattern)
        bu_pattern = args.bu_pattern
    if args.bclist:
        print("Set bclist file: "+args.bclist)
        bclist = args.bclist
    if args.genomeDir:
        print("Set geonome directory: ",args.genomeDir)
        genomeDir = args.genomeDir
    if args.barcode_match:
        print("Set match barcode in NGS: ",args.barcode_match)
        barcode_match = args.barcode_match
    if args.steps:
        if args.steps == 'all':
            step_list = ["ccs","remove_primer","pattern_detection","flc","split_linker","dedup","featurecount","count","run_seurat","run_isoform"]
        else:
            step_list = args.steps.strip().split(',')
        print("Set steps: ",step_list)

    PA = pacbio_analysis()
    PA.run_pacbio()

    