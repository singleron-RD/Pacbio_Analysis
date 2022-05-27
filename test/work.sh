#!usr/bin/bash

#make test bam file
samtools view -bh --subsample 0.01 -o test.bam m64236_211121_102448.hifi_reads.bam

#generate script for step1
bash ../src/run.s1.sh test_sample /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test/test.bam /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test

#generate script for step2
bash ../src/run.s2.sh test_sample /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test

#generate script for step3
bash ../src/run.s3.sh test_sample /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test

#generate script for step4
bash ../src/run.s4.sh test_sample /SGRNJ03/randd/user/fuxin/Github_repo/Pacbio_Analysis/test
