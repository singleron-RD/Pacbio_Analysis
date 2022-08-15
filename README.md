# Pacbio_Analysis
Third-generation single cell sequencing (Pacbio) analysis pipeline.

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

`--steps`   set to run specific steps, optional steps include: ccs, remove.primers, pattern.bc.umi, flc, split.linker, dedup, featurecount, count, run_seurat, run_isoform 



##### Exsample
```
conda activate pacbio_v1
python src/run_pacbio.py 
    --ccs_bam  test.bam
    --sample test_sample
    --outdir test_res_1
    --steps all
```

To analysis isoform only:

```
conda activate pacbio_v1
python src/run_pacbio.py 
    --ccs_bam  test.bam
    --sample test_sample
    --outdir test_res_1
    --steps ccs,remove.primers,pattern.bc.umi,flc,split.linker,dedup,run_isoform
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





