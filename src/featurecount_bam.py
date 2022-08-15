#generate bam file for aligned sequence, and add gene id 
#input: 1.featurecount CORE file; 2.original sam file
#output: _tmp.sam

import re
from collections import Counter, defaultdict
import pysam
import sys

import gzip

import sys
from collections import Counter, defaultdict
from datetime import timedelta
from functools import wraps

from Bio.Seq import Seq
import pandas as pd
import pysam

def generic_open(file_name, *args, **kwargs):
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj

class Gtf_dict(dict):
    '''
    key: gene_id
    value: gene_name
    If the key does not exist, return key. This is to avoid the error:
        The gtf file contains one exon lines with a gene_id, but do not contain a gene line with the same gene_id. FeatureCounts 
        work correctly under this condition, but the gene_id will not appear in the Gtf_dict.
    '''

    def __init__(self, gtf_file):
        super().__init__()
        self.gtf_file = gtf_file
        self.load_gtf()

    def load_gtf(self):
        """
        get gene_id:gene_name from gtf file
            - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
            - one gene_id with multiple gene_name: error.
            - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
            - no gene_name: gene_id will be used as gene_name.

        Returns:
            {gene_id: gene_name} dict
        """

        gene_id_pattern = re.compile(r'gene_id "(\S+)";')
        gene_name_pattern = re.compile(r'gene_name "(\S+)"')
        id_name = {}
        c = Counter()
        with generic_open(self.gtf_file) as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    continue
                tabs = line.split('\t')
                gtf_type, attributes = tabs[2], tabs[-1]
                if gtf_type == 'gene':  
                    gene_id = gene_id_pattern.findall(attributes)[-1]
                    gene_names = gene_name_pattern.findall(attributes)
                    if not gene_names:
                        gene_name = gene_id
                    else:
                        gene_name = gene_names[-1]
                    c[gene_name] += 1
                    if c[gene_name] > 1:
                        if gene_id in id_name:
                            assert id_name[gene_id] == gene_name, (
                                'one gene_id with multiple gene_name '
                                f'gene_id: {gene_id}, '
                                f'gene_name this line: {gene_name}'
                                f'gene_name previous line: {id_name[gene_id]}'
                            )
                            self.load_gtf.logger.warning(
                                'duplicated (gene_id, gene_name)'
                                f'gene_id: {gene_id}, '
                                f'gene_name {gene_name}'
                            )
                            c[gene_name] -= 1
                        else:
                            gene_name = f'{gene_name}_{c[gene_name]}'
                    id_name[gene_id] = gene_name
        self.update(id_name)

    def __getitem__(self, key):
        '''if key not exist, return key'''
        return dict.get(self, key, key)

def get_assign_dict(featurecount_file):
    featurecount_assign_dict = {}
    with open(featurecount_file) as infile:
        for line in infile:
            nline = line.strip().split('\t')
            featurecount_assign_dict[nline[0]] = {'status' : nline[1],'flag':nline[2],'gene_id':nline[3]}
    return(featurecount_assign_dict)

#def get_count_dict(featurecount_count_file):
#    featurecount_count_dict = {}  



def add_tag(gtf_file,featurecount_file,minimap_sam):
    """
    - CB cell barcode
    - UB UMI
- GN gene name
    - GX gene id
    """

    temp_sam_file = featurecount_file+".temp_sam"

    gtf_dict = Gtf_dict(gtf_file)
    featurecount_assign_dict = get_assign_dict(featurecount_file)

    with pysam.AlignmentFile(minimap_sam, "rb") as original_sam:
        header = original_sam.header
        with pysam.AlignmentFile(temp_sam_file, "w", header=header) as temp_sam:
            for read in original_sam:
                if featurecount_assign_dict[read.query_name]['status'] == 'Assigned':
                    attr = read.query_name.split('_')
                    barcode = attr[0]
                    umi = attr[1]
                    read.set_tag(tag='CB', value=barcode, value_type='Z')
                    read.set_tag(tag='UB', value=umi, value_type='Z')
                    read.set_tag(tag='XT', value=featurecount_assign_dict[read.query_name]['gene_id'], value_type='Z')
                    read.set_tag(tag='XS', value=featurecount_assign_dict[read.query_name]['status'], value_type='Z')
                    #read.set_tag(tag='XN', value=featurecount_count_dict[gene_id]['count'], value_type='i')                    
                    # assign to some gene
                    if read.has_tag('XT'):
                        gene_id = read.get_tag('XT')
                        # if multi-mapping reads are included in original bam,
                        # there are multiple gene_ids
                        if ',' in gene_id:
                            gene_name = [gtf_dict[i] for i in gene_id.split(',')]
                            gene_name = ','.join(gene_name)
                        else:
                            gene_name = gtf_dict[gene_id]
                        read.set_tag(tag='GN', value=gene_name, value_type='Z')
                        read.set_tag(tag='GX', value=gene_id, value_type='Z')
                    temp_sam.write(read)

add_tag(gtf_file = sys.argv[1], featurecount_file = sys.argv[2], minimap_sam=sys.argv[3])