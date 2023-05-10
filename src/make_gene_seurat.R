library(tidyverse)
library(Seurat)
library(Matrix)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--known_isoform_seurat", help="known_isoform_seurat path")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--gene_name_set", help="Correspondence file between gene symbol and Ensembl id, separated by tab")
argv <- parse_args(argv)

outdir = argv$outdir
indir = argv$known_isoform_seurat
hgnc_gene_set = argv$gene_name_set

print(indir)
gene_df = readr::read_tsv(hgnc_gene_set)%>%select(symbol,ensembl_gene_id)

genes = readr::read_tsv(paste0(indir,"/genes.tsv"),
                        col_names = FALSE)%>%
    left_join(gene_df,by = c("X2" = "ensembl_gene_id"))


mtx = Read10X(indir,gene.column = 1)

data.frame(mtx)%>%mutate(pb_id = rownames(mtx))%>%
    left_join(genes,by = c("pb_id" = "X1"))%>%
    filter(!(is.na(symbol)))->mtx_df_1

aggregate(x = mtx_df_1%>%select(-X2,-pb_id,-symbol), 
    by=list(symbol = mtx_df_1$symbol), 
    FUN=sum, 
    simplify = TRUE)->mtx_df_2

rownames(mtx_df_2) <- mtx_df_2$symbol
data.frame(symbol = mtx_df_2$symbol)%>%left_join(gene_df)%>%select(ensembl_gene_id,symbol)->gene_df_2
mtx_df_2%>%select(-symbol) ->mtx_df_2

sparse_mm = Matrix(as.matrix(mtx_df_2),sparse = TRUE)
write.table(gene_df_2,paste0(outdir,"/genes.tsv"),sep = "\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
write(x = sparse_mm@Dimnames[[2]],file = paste0(outdir,"/barcodes.tsv"))
writeMM(obj = sparse_mm,file =  paste0(outdir,"/matrix.mtx"))