##### Download [ASprofile](https://ccb.jhu.edu/software/ASprofile/)

`wget https://ccb.jhu.edu/software/ASprofile/ASprofile.tar.gz`

`gunzip < ASprofile.tar.gz | tar -xvf -`

##### Generate referen genome for AS profile

`python <Pacbio_analysis_path>/src/ASprofile/get_ASprofile_ref_hdrs.py path/species.genome.fa species`

##### Extract AS with ASprofile

`ASprofile.b-1.0.4/extract-as <pacbio_analysis_outpath>/10.isoform/sample.collapsed_classification.filtered_lite.gtf sapiens.fa.hdrs > sample_as_events.txt`

`perl ASprofile.b-1.0.4/summarize_as.pl <pacbio_analysis_outpath>/10.isoform/sample.collapsed_classification.filtered_lite.gtf sample_as_events.txt -p sample`

##### Downstream analysis

Refer to [Statistics_of_structure_features_byASprofile](https://github.com/singleron-RD/Pacbio_Analysis/blob/main/notebooks/Statistics_of_structure_features_byASprofile.ipynb)




