python ../src/run_pacbio.py --sample test_sample --outdir test_res_all/ --ccs_bam test.bam --bu_pattern T-12U-57B --steps all

python ../src/run_pacbio.py --sample test_sample --outdir test_res_isoform/ --ccs_bam test.bam --bu_pattern T-12U-57B --steps ccs,remove.primers,pattern.bc.umi,flc,split.linker,dedup,run_isoform
