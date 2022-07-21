# Pacbio_Analysis
Third-generation single cell sequencing (Pacbio) analysis pipeline.

## Install

###### 克隆仓库
`git clone git@github.com:singleron-RD/Pacbio_Analysis.git`

###### 创建 conda 环境，安装依赖包
`cd Pacbio_Analysis`
`conda create -n pacbio_v1 -y --file conda_pkgs.txt`
###### 安装依赖软件
[Celescope](https://github.com/singleron-RD/CeleScope/blob/master/docs/installation.md)
[cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake)
[SQANTI3](https://github.com/ConesaLab/SQANTI3)

## Usage

###### 命令示例
`python src/run_pacbio.py --ccs_bam input_bam_file --sample sample_name <optional parameters> `

###### 参数说明

```parameters:
    --ccs_bam   输入的bam文件(required)
    --sample    样本名称/输出文件前缀(required)
    --outdir    输出文件目录
    --primer_fasta	自定义引物序列文件，fasta格式
    --bu_pattern	自定义Barcode与UMI文件，默认T-12U-57B
    --blastn	自定义blastn工具路径
    --bclist	自定义barcode序列文件
    --genomeDir	自定义参考基因组，默认人类基因组ensembl 92
```





