python ../src/run_pacbio.py --sample test_sample --outdir test_res_all/ --ccs_bam /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/test.bam --bu_pattern T-12U-57B --steps all

python ../src/run_pacbio.py --sample test_sample --outdir test_res_isoform/ --ccs_bam /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/test.bam --bu_pattern T-12U-57B --steps ccs,remove_primer,pattern_detection,flc,split_linker,dedup,run_isoform
