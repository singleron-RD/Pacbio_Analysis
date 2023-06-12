# Pacbio_Analysis
Third-generation single cell sequencing (Pacbio) analysis pipeline.

## Introduction

This pipeline is developed to analysis third-generation sequencing data produced by Pacbio platform. It takes CCS bam file as input file and generate output files which can be used for downstream data analysis as well as a summary QC statistic file and two html report of transcriptome and alternative splicing atlas.

The two main functions of this pipeline are to analyze transcript characteristics and alternative splicing characteristics of Pacbio dataset. The pipeline consists of several steps and the first 6 steps were identical pre-processing steps, refer to [ISOseq guidline](https://github.com/Magdoll/cDNA_Cupcake/wiki/Iso-Seq-Single-Cell-Analysis:-Recommended-Analysis-Guidelines):

* `ccs` HiFi reads which produced using circular consensus sequencing (CCS) mode on PacBio long-read systems are took as input file and this step generate the directory structure for CCS bam file.

* `remove_primer`   [lima](https://lima.how/) is used to remove the 5' and 3' cDNA primers.

* `pattern_detection` UMIs and cell barcodes are clipped from the reads and associated with the reads for later deduplication, refer to [isoseq3 tag](https://isoseq.how/umi/cli-workflow.html#step-3---tag).

* `flc` [isoseq3 refine](https://isoseq.how/umi/cli-workflow.html) is used to remove polyA tail and artificial concatemers.

* `split_linker`    Remove linker sequence in barcode.

* `dedup`   This step performs PCR deduplicatation via clustering by UMI and cell barcodes, refer to [isoseq3 dedup](https://isoseq.how/umi/cli-workflow.html#step-5---deduplication).

To analysis transcript characteristics, refer to [celescope](https://github.com/singleron-RD/CeleScope):

* `featurecount`    Assigning uniquely mapped reads to genomic features with FeatureCounts.

* `count` Distinguish cell barcodes from background barcodes and Generate expression matrix.

* `run_seurat`  Cell clustering with Seurat and Calculate the marker gene of each cluster.

To analysis alternative splicing characteristics:

* `run_isoform` Applying [Cupcake](https://github.com/Magdoll/cDNA_Cupcake/) to collapse the mapped dedup-ed reads into unique transcripts and [SQANTI3](https://github.com/ConesaLab/SQANTI3) to annotate each unique transcript against an annotation.

## Installation

###### Version

v1.1


###### Create conda env

It is recommended to use mamba, [download and install mamba](https://mamba.readthedocs.io/en/latest/install).

```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

###### Clone git and install dependency software

`git clone https://github.com/singleron-RD/Pacbio_Analysis.git`

`cd Pacbio_Analysis`

`mamba create -n pacbio_isoseq -y --file conda_pkgs.txt`

`mamba activate pacbio_isoseq`

###### Install [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake)

`cd Pacbio_Analysis/tools/cDNA_Cupcake`

`python setup.py build`

`python setup.py install`

`export PYTHONPATH=$PYTHONPATH:<path_to_Cupcake>/sequence/`

## Usage

##### Command

`python src/run_pacbio.py --ccs_bam input_bam_file -- genomeDir genome_directory --sample sample_name <optional parameters>`

##### Parameters

`--ccs_bam` CCS bam file, required

`--genomeDir`   Reference genome directory, required

`--sample`  Sample name/Output file prefix, required

`--pacbio_source_path` Absolute path to Pacbio_analysis

`--steps`   set specific steps to run, "all" to run all the steps, optional steps include: ccs, remove_primer, pattern_detection, flc, split_linker, dedup, featurecount, count, run_seurat, run_isoform, required

`--outdir`  Output file directory, default current

`--primer_fasta`    Primer sequence file, fasta fomat, optional   

`--bu_pattern`  Barcord and UMI pattern expression, default T-12U-57B

`--blastn`  Blastn tool path

`--bclist`  Barcode sequence file, used in linker_split step

`--report`   Whether to generate summary statistical report, default True



##### Exsample

To analysis isoform only (Recommend):

```
cd test
mkdir test_res_isoform
conda activate pacbio_v1
python src/run_pacbio.py 
    --ccs_bam  <path>/test.bam
    --genomeDir genome_directory
    --pacbio_source_path <path>/Pacbio_Analysis
    --sample test_sample
    --outdir test_res_isoform
    --steps ccs,remove_primer,pattern_detection,flc,split_linker,dedup,run_isoform
```

To run all the steps (Calculate expression matrix with Celescope):
```
cd test
mkdir test_res_all
conda activate pacbio_v1
python src/run_pacbio.py 
    --ccs_bam  <path>/test.bam
    --genomeDir genome_directory
    --pacbio_source_path <path>/Pacbio_Analysis
    --sample test_sample
    --outdir test_res_isoform
    --steps all
```

## Output

If all the steps completed, we'll get:

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

## Notice

Provide the absolute path of ccs bam file to avoid unexpected errors.





