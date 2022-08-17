# Pacbio_Analysis
Third-generation single cell sequencing (Pacbio) analysis pipeline.

## Introduction

This pipeline is developed to analysis third-generation sequencing data produced by Pacbio platform. It takes CCS bam file as input file and generate output files which can be used for downstream data analysis as well as a summary QC statistic file and two html report of transcriptome and alternative splicing atlas.

The two main functions of this pipeline are to analyze transcript characteristics and alternative splicing characteristics of Pacbio dataset. The pipeline consists of several steps and the first 6 steps were identical pre-processing steps, refer to [ISOseq guidline](https://github.com/Magdoll/cDNA_Cupcake/wiki/Iso-Seq-Single-Cell-Analysis:-Recommended-Analysis-Guidelines):
`ccs` HiFi reads which produced using circular consensus sequencing (CCS) mode on PacBio long-read systems are took as input file and this step generate the directory structure for CCS bam file.

`remove_primer` [lima](https://lima.how/) is used to remove the 5' and 3' cDNA primers.

`pattern_detection` UMIs and cell barcodes are clipped from the reads and associated with the reads for later deduplication, refer to [isoseq3 tag](https://isoseq.how/umi/cli-workflow.html#step-3---tag).

`flc` [isoseq3 refine](https://isoseq.how/umi/cli-workflow.html) is used to remove polyA tail and artificial concatemers.

`split_linker` Remove linker sequence in barcode.

`dedup` This step performs PCR deduplicatation via clustering by UMI and cell barcodes, refer to [isoseq3 dedup](https://isoseq.how/umi/cli-workflow.html#step-5---deduplication).

To analysis transcript characteristics, refer to [celescope](https://github.com/singleron-RD/CeleScope):
`featurecount` Assigning uniquely mapped reads to genomic features with FeatureCounts.

`count` Distinguish cell barcodes from background barcodes and Generate expression matrix.

`run_seurat` Cell clustering with Seurat and Calculate the marker gene of each cluster.

To analysis alternative splicing characteristics:
`run_isoform`   Applying [Cupcake](https://github.com/Magdoll/cDNA_Cupcake/) to collapse the mapped dedup-ed reads into unique transcripts and [SQANTI3](https://github.com/ConesaLab/SQANTI3) to annotate each unique transcript against an annotation.

## Installation

###### Clone git
`git clone git@github.com:singleron-RD/Pacbio_Analysis.git`

###### Create conda env
`cd Pacbio_Analysis`

`conda create -n pacbio_v1 -y --file conda_pkgs.txt`

###### Install dependency software
[Celescope](https://github.com/singleron-RD/CeleScope/blob/master/docs/installation.md)

[cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake)

[SQANTI3](https://github.com/ConesaLab/SQANTI3)

## Usage

##### Command

`python src/run_pacbio.py --ccs_bam input_bam_file --sample sample_name <optional parameters>`

##### Parameters

`--ccs_bam` CCS bam file, required

`--sample`  Sample name/Output file prefix, required

`--outdir`  Output file directory

`--primer_fasta`    primer sequence file, fasta fomat, optional   

`--bu_pattern`  Barcord and UMI pattern expression, default T-12U-57B

`--blastn`  Blastn tool path

`--bclist`  Barcode sequence file

`--genomeDir`   Reference genome directory, default human genome ensembl 92

`--steps`   set to run specific steps, optional steps include: ccs, remove_primer, pattern_detection, flc, split_linker, dedup, featurecount, count, run_seurat, run_isoform 



##### Exsample
```
conda activate pacbio_v1
python src/run_pacbio.py 
    --ccs_bam  test.bam
    --sample test_sample
    --outdir test_res_all
    --steps all
```

To analysis isoform only:

```
conda activate pacbio_v1
python src/run_pacbio.py 
    --ccs_bam  test.bam
    --sample test_sample
    --outdir test_res_isoform
    --steps ccs,remove_primer,pattern_detection,flc,split_linker,dedup,run_isoform
```

## Output

If all the steps finished, we'll get:

```
outdir
|--01.ccs
|--02.remove.primers
|--03.pattern.bc.umi
|--04.flc
|--05.split.linker
|--06.dedup
|--07.featurecount
|--08.count
|--09.seurat
|--10.isoform
|--11.annotation
|--sample.Pacbio_qc.stat.txt
|--sample_report.html
|--sample.SQANTI3_report.html
```





