# coding=utf-8
# pzw
# 20241220
# 创建autocnv的数据库，又不是不能用😁

from multiprocessing import Pool
import pandas as pd
import numpy as np
import allel
import os
pd.set_option('display.max_columns', None)
from tqdm import tqdm
from gtfparse import read_gtf

## 用到的软件
bgzip = 'bgzip'
tabix = 'tabix'

## 下载这些数据库
exon_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz'
cyto_band_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz'
gene_info_url = 'https://ftp.ncbi.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz'
ref_gene_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
clingen_gene_curation_url = 'https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv'
clingene_region_curation_url = 'https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv'
clinvar_url = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
hi_pred_url = 'https://www.deciphergenomics.org/files/downloads/HI_Predictions_Version3.bed.gz'
gnomad_lof_url = 'https://datasetgnomad.blob.core.windows.net/dataset/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'
gnomad_control_only_url = 'https://datasetgnomad.blob.core.windows.net/dataset/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.bed.gz'
hgnc_gene_fam_url = 'https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/family.csv'
DGV_GS_Gain_Loss_url = 'https://dgv.tcag.ca/dgv/docs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3'
omim_gene_list_url = 'https://www.omim.org/static/omim/data/mim2gene.txt'

# 将下载下来的数据库统一放在raw_data文件夹下
raw_data_dir = 'raw_data'
result_data_dir = 'data'
exon_ori_file = os.path.join(raw_data_dir, 'hg19.refGene.gtf.gz')
cyto_band_ori_file = os.path.join(raw_data_dir, 'cytoBand.txt.gz')
gene_info_ori_file = os.path.join(raw_data_dir, 'Homo_sapiens.gene_info.gz')
ref_gene_ori_file = os.path.join(raw_data_dir, 'refGene.txt.gz')
clingen_gene_ori_file = os.path.join(raw_data_dir, 'ClinGen_gene_curation_list_GRCh37.tsv')
clingen_region_ori_file = os.path.join(raw_data_dir, 'ClinGen_region_curation_list_GRCh37.tsv')
clinvar_ori_vcf_file = os.path.join(raw_data_dir, 'clinvar.vcf.gz')
hi_pred_ori_file = os.path.join(raw_data_dir, 'HI_Predictions_Version3.bed.gz')
gnomad_lof_ori_file = os.path.join(raw_data_dir, 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
gnomad_control_ori_file = os.path.join(raw_data_dir, 'gnomad_v2.1_sv.controls_only.sites.bed.gz')

# 准备外部omim基因列表
omim_gene_list_file = os.path.join(raw_data_dir, 'mim2gene.txt')
dgv_ori_file = os.path.join(raw_data_dir, 'DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3')

# 输出文件
gene_file = os.path.join(result_data_dir, 'gene.sorted.bed')
omim_gene_file = os.path.join(result_data_dir, 'omim-gene.sorted.bed')
clinvar_file = os.path.join(result_data_dir, 'clinvar-pathogenic.sorted.vcf')
decipher_gene_file = os.path.join(result_data_dir, 'decipher-gene.sorted.bed')
dgv_gain_file = os.path.join(result_data_dir, 'dgv-gain.sorted.bed')
dgv_loss_file = os.path.join(result_data_dir, 'dgv-loss.sorted.bed')
func_region_file = os.path.join(result_data_dir, 'func-region.sorted')
gnomad_del_file = os.path.join(result_data_dir, 'gnomad-del.sorted.bed')
gnomad_dup_file = os.path.join(result_data_dir, 'gnomad-dup.sorted.bed')
hi_cds_file = os.path.join(result_data_dir, 'hi-cds.sorted.bed')
hi_exon_file = os.path.join(result_data_dir, 'hi-exon.sorted.bed')
hi_gene_file = os.path.join(result_data_dir, 'hi-gene.sorted.bed')
hi_region_file = os.path.join(result_data_dir, 'hi-region.sorted.bed')
ts_gene_file = os.path.join(result_data_dir, 'ts-gene.sorted.bed')
ts_region_file = os.path.join(result_data_dir, 'ts-region.sorted.bed')
uhi_gene_file = os.path.join(result_data_dir, 'uhi-gene.sorted.bed')
uhi_region_file = os.path.join(result_data_dir, 'uhi-region.sorted.bed')
uts_gene_file = os.path.join(result_data_dir, 'uts-gene.sorted.bed')
uts_region_file = os.path.join(result_data_dir, 'uts-region.sorted.bed')
hgnc_gene_fam_file = os.path.join(result_data_dir, 'family.csv')

#########################################

# 数据库制作
def cytoband():
    cyto_band_df = pd.read_csv(
        cyto_band_ori_file, sep='\t',
        names=['#chrom', 'start', 'end', 'name', 'gieStain']
    )
    cyto_band_df.to_csv(os.path.join(result_data_dir, "cyto-band.bed"), sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {result_data_dir}/cyto-band.bed > {result_data_dir}/cyto-band.bed.gz
        {tabix} -fp bed {result_data_dir}/cyto-band.bed.gz
        rm {result_data_dir}/cyto-band.bed
    """
    os.system(cmd)


def fetch_variants(region):
    need_fields = [
        'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
        'variants/AF_ESP', 'variants/AF_EXAC', 'variants/AF_TGP', 'variants/CLNSIG'
    ]
    fields, samples, headers, it = allel.iter_vcf_chunks(
        clinvar_ori_vcf_file, fields=need_fields, alt_number=1, region=region
    )
    for variants, *_ in it:
        esp_filter = np.isnan(variants['variants/AF_ESP'])
        esp_filter[~esp_filter] |= variants['variants/AF_ESP'][~esp_filter] < 0.01

        exac_filter = np.isnan(variants['variants/AF_EXAC'])
        exac_filter[~exac_filter] |= variants['variants/AF_EXAC'][~exac_filter] < 0.01

        tgp_filter = np.isnan(variants['variants/AF_TGP'])
        tgp_filter[~tgp_filter] |= variants['variants/AF_TGP'][~tgp_filter] < 0.01

        pathogenic_filter = np.isin(
            variants['variants/CLNSIG'], ['Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic']
        )

        af_filter = esp_filter & exac_filter & tgp_filter & pathogenic_filter

        filtered_variants = {k: v[af_filter] for k, v in variants.items()}

        filtered_variants['variants/CHROM'] = 'chr' + filtered_variants['variants/CHROM']

        return allel.normalize_callset(filtered_variants)

def fetch_gene(gene, chrom, start, end):
    return ','.join(gene.loc[
        (gene['chrom'] == chrom) & (gene['txEnd'] >= start) & (gene['txStart'] <= end), 'name2'
    ])

#
def gene_info():
    cols = [
        'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
        'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'ExonFrames']
    refgene = pd.read_csv(ref_gene_ori_file, sep='\t', names=cols)
    refgene = refgene[~refgene['chrom'].str.match(r'.*fix$')]
    refgene['length'] = refgene['cdsEnd'] - refgene['cdsStart']
    refgene = refgene.sort_values('length', ascending=False)
    refgene = refgene.drop_duplicates('name2', keep='first')
    gene_info = pd.read_csv(gene_info_ori_file, sep='\t')
    refgene_info = refgene.merge(gene_info, left_on='name2', right_on='Symbol')
    gene_fam = pd.read_csv(hgnc_gene_fam_file, sep=',')
    gene_fam = gene_fam[gene_fam['typical_gene'].notna()].rename(columns={'id': 'fam'})
    gene_fam_grp = gene_fam.groupby('typical_gene').agg({'fam': set})
    gene = refgene_info.loc[
        refgene_info['type_of_gene'] == 'protein-coding',
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    gene.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(gene_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {gene_file} > {gene_file}.gz
        {tabix} -fp bed {gene_file}.gz
        rm {gene_file}
    """
    os.system(cmd)


    omim_gene_list = []
    with open(omim_gene_list_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith('#'):
                lines = line.rstrip().split('\t')
                if lines[1] != "gene":
                    continue
                if len(lines) < 4:
                    continue
                gene_check = lines[3]
                if not gene_check == "":
                    omim_gene_list.append(gene_check)

    omim_gene = set(omim_gene_list)
    omim_gene = refgene_info.loc[
        refgene_info['name2'].isin(omim_gene),
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    omim_gene.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(omim_gene_file, index=False, sep='\t')

    cmd = f"""
        {bgzip} -cf {omim_gene_file} > {omim_gene_file}.gz
        {tabix} -fp bed {omim_gene_file}.gz
        rm {omim_gene_file}
    """
    os.system(cmd)

    curation_gene = pd.read_csv(clingen_gene_ori_file, sep='\t', dtype=str, skiprows=5)
    hi_genes = set(
    curation_gene.loc[
        curation_gene['Haploinsufficiency Score'] == '3', '#Gene Symbol'
    ])
    hi_gene = refgene_info.loc[
        refgene_info['name2'].isin(hi_genes),
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    hi_gene.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(hi_gene_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {hi_gene_file} > {hi_gene_file}.gz
        {tabix} -fp bed {hi_gene_file}.gz
        rm {hi_gene_file}
    """
    os.system(cmd)

    hi_cds = refgene_info.loc[
        (refgene_info['name2'].isin(hi_genes)) & (refgene_info['length'] != 0),
        ['chrom', 'cdsStart', 'cdsEnd', 'GeneID', 'name2', 'name']
    ].sort_values(['chrom', 'cdsStart', 'cdsEnd'])
    hi_cds.rename(columns={
        'chrom': '#chrom', 'cdsStart': 'start', 'cdsEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(hi_cds_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {hi_cds_file} > {hi_cds_file}.gz
        {tabix} -fp bed {hi_cds_file}.gz
        rm {hi_cds_file}
    """
    os.system(cmd)

    hi_exons = refgene_info.loc[
        refgene_info['name2'].isin(hi_genes), ['chrom', 'exonStarts', 'exonEnds', 'GeneID', 'name2', 'name', 'strand']
    ].copy()
    hi_exons['exonStarts'] = hi_exons['exonStarts'].str.replace(r',$', '')
    hi_exons['exonEnds'] = hi_exons['exonEnds'].str.replace(r',$', '')
    start = hi_exons['exonStarts'].str.split(',').apply(pd.Series).stack().reset_index()
    start = start.rename(columns={'level_0': 'row', 0: 'start'})[['row', 'start']]
    start['start'] = start['start'].astype(int)
    end = hi_exons['exonEnds'].str.split(',').apply(pd.Series).stack().reset_index()
    end = end.rename(columns={0: 'end'})['end'].astype(int)
    position = start.join(end)
    exon = position.merge(
        hi_exons[['chrom', 'GeneID', 'name2', 'name', 'strand']], how='left', left_on='row', right_index=True
    )
    exon = exon.sort_values(['chrom', 'start', 'end'])
    exon['+'] = exon.groupby(['name2', 'name'])['start'].rank('first', ascending=True).astype(int)
    exon['-'] = exon.groupby(['name2', 'name'])['start'].rank('first', ascending=False).astype(int)
    exon['exon'] = pd.concat([exon.loc[exon['strand'] == '+', '+'], exon.loc[exon['strand'] == '-', '-']])
    exon['last_exon'] = exon.groupby(['name2', 'name'])['exon'].transform('max') == exon['exon']
    exon = exon[
        ['chrom', 'start', 'end', 'GeneID', 'name2', 'name', 'exon', 'last_exon']
    ].sort_values(['chrom', 'start', 'end'])
    exon.rename(columns={
        'chrom': '#chrom',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(hi_exon_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {hi_exon_file} > {hi_exon_file}.gz
        {tabix} -fp bed {hi_exon_file}.gz
        rm {hi_exon_file}
    """
    os.system(cmd)

    last_exon = exon[exon['last_exon']]
    last_exon_region = last_exon['chrom'] + ':' + last_exon['start'].astype(str) + '-' + last_exon['end'].astype(str)
    last_exon_region = last_exon_region.str.replace('chr', '')
    need_fields = [
        'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
        'variants/AF_ESP', 'variants/AF_EXAC', 'variants/AF_TGP', 'variants/CLNSIG'
    ]
    cmd = f"""
        {tabix} -f {clinvar_ori_vcf_file}
    """
    os.system(cmd)
    with open(clinvar_file, 'w') as f:
        headers = allel.read_vcf_headers(clinvar_ori_vcf_file)
        f.write(''.join(headers.headers))
                
        with Pool(processes=8) as pool:
            variants = pool.map(fetch_variants, last_exon_region)
        
        for names, callset in tqdm(filter(lambda x: x is not None, variants)):
            allel.write_vcf_data(f, names, callset, None, {field: np.nan for field in need_fields})
    
    cmd = f"""
        sed -i 's/AF_.\+=nan;//' {clinvar_file}
        {bgzip} -cf {clinvar_file} > {clinvar_file}.gz
        {tabix} -fp vcf {clinvar_file}.gz
        rm {clinvar_file}
    """
    os.system(cmd)

    uhi_genes = set(
        curation_gene.loc[curation_gene['Haploinsufficiency Score'] == '40', '#Gene Symbol']
    )
    uhi_gene = refgene_info.loc[
        refgene_info['name2'].isin(uhi_genes),
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    uhi_gene.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(uhi_gene_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {uhi_gene_file} > {uhi_gene_file}.gz
        {tabix} -fp bed {uhi_gene_file}.gz
        rm {uhi_gene_file}
    """
    os.system(cmd)

    region = pd.read_csv(clingen_region_ori_file, sep='\t', skiprows=5, dtype=str)
    position = region['Genomic Location'].str.extract(r'(?P<chrom>chr\w+)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)')
    position = position.dropna()
    position['start'] = position['start'].astype(int)
    position['end'] = position['end'].astype(int)
    region_pos = region.merge(position, how='left', left_index=True, right_index=True)
    func_region = region_pos.loc[
        (region_pos['Haploinsufficiency Score'].isin(['1', '2','3']))
        | (region_pos['Triplosensitivity Score'].isin(['1', '2','3'])),
        ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
    ].sort_values(['chrom', 'start', 'end'])
    func_region = func_region.dropna(subset=['start', 'end'])
    func_region['start'] = func_region['start'].astype(int)
    func_region['end'] = func_region['end'].astype(int)
    
    # 清理数据：去除多余的空格，确保染色体名称格式正确
    func_region['#ISCA ID'] = func_region['#ISCA ID'].str.strip()
    func_region['ISCA Region Name'] = func_region['ISCA Region Name'].str.strip()
    func_region['chrom'] = func_region['chrom'].str.strip()
    
    func_region.rename(columns={
        'chrom': '#chrom',
        '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
    }).to_csv(func_region_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {func_region_file} > {func_region_file}.gz  
        {tabix} -fp bed {func_region_file}.gz
        rm {func_region_file}
    """
    os.system(cmd)

    hi_region = region_pos.loc[
        region_pos['Haploinsufficiency Score'] == '3',
        ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
    ].sort_values(['chrom', 'start', 'end'])
    hi_region['omim_genes'] = hi_region.apply(
        lambda row: ','.join(omim_gene.loc[
            (omim_gene['chrom'] == row['chrom'])
            & (omim_gene['txEnd'] >= row['start'])
            & (omim_gene['txStart'] <= row['end']),
            'name2'
        ]),
        axis=1
    )
    hi_region = hi_region.dropna(subset=['start', 'end'])
    hi_region['start'] = hi_region['start'].astype(int)
    hi_region['end'] = hi_region['end'].astype(int)
    hi_region.rename(columns={
        'chrom': '#chrom',
        '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
    }).to_csv(hi_region_file, sep='\t', index=False)

    cmd = f"""
        {bgzip} -cf {hi_region_file} > {hi_region_file}.gz
        {tabix} -fp bed {hi_region_file}.gz
        rm {hi_region_file}
    """
    os.system(cmd)

    uhi_region = region_pos.loc[
        region_pos['Haploinsufficiency Score'] == '40',
        ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
    ].sort_values(['chrom', 'start', 'end'])
    uhi_region['genes'] = uhi_region.apply(
        lambda row: ','.join(gene.loc[
            (gene['chrom'] == row['chrom'])
            & (gene['txEnd'] >= row['start'])
            & (gene['txStart'] <= row['end']),
            'name2'
        ]),
        axis=1
    )
    uhi_region['start'] = uhi_region['start'].astype(int)
    uhi_region['end'] = uhi_region['end'].astype(int)
    uhi_region.rename(columns={
        'chrom': '#chrom',
        '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
    }).to_csv(uhi_region_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {uhi_region_file} > {uhi_region_file}.gz
        {tabix} -fp bed {uhi_region_file}.gz
        rm {uhi_region_file}
    """
    os.system(cmd)

    decipher = pd.read_csv(hi_pred_ori_file, sep='\t',skiprows=1, header=None, usecols=[3,])
    decipher = decipher[3].str.split('|', expand=True).rename(columns={0: 'symbol', 1: 'hi_score', 2: 'hi_index'})
    decipher['hi_index'] = decipher['hi_index'].str.replace('%', '').astype(float)
    decipher = decipher.merge(gene, left_on='symbol', right_on='name2')
    gnomad = pd.read_csv(gnomad_lof_ori_file, sep='\t', index_col=0, compression='gzip')
    decipher = decipher.join(gnomad['pLI'], on='name2')
    decipher = decipher.join(gnomad['oe_lof_upper'], on='name2')
    decipher = decipher.loc[
        (decipher['pLI'] >= 0.9) & (decipher['hi_index'] < 10) & (decipher['oe_lof_upper'] < 0.35),
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'pLI', 'hi_score']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    decipher.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(decipher_gene_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {decipher_gene_file} > {decipher_gene_file}.gz  
        {tabix} -fp bed {decipher_gene_file}.gz
        rm {decipher_gene_file}
    """
    os.system(cmd)

    ts_genes = set(
        curation_gene.loc[
            curation_gene['Triplosensitivity Score'] == '3', '#Gene Symbol'
        ]
    )
    ts_gene = refgene_info.loc[
        refgene_info['name2'].isin(ts_genes),
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    ts_gene.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(ts_gene_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {ts_gene_file} > {ts_gene_file}.gz
        {tabix} -fp bed {ts_gene_file}.gz
        rm {ts_gene_file}
    """
    os.system(cmd)

    uts_genes = set(
        curation_gene.loc[curation_gene['Triplosensitivity Score'] == '40', '#Gene Symbol']
    )
    uts_gene = refgene_info.loc[
        refgene_info['name2'].isin(uts_genes),
        ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
    ].sort_values(['chrom', 'txStart', 'txEnd'])
    uts_gene.rename(columns={
        'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
        'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
    }).to_csv(uts_gene_file, index=False, sep='\t')
    cmd = f"""
        {bgzip} -cf {uts_gene_file} > {uts_gene_file}.gz
        {tabix} -fp bed {uts_gene_file}.gz
        rm {uts_gene_file}
    """
    os.system(cmd)

    ts_region = region_pos.loc[
        region_pos['Triplosensitivity Score'] == '3',
        ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
    ].sort_values(['chrom', 'start', 'end'])
    ts_region['omim_genes'] = ts_region.apply(
        lambda row: ','.join(omim_gene.loc[
            (omim_gene['chrom'] == row['chrom'])
            & (omim_gene['txEnd'] >= row['start'])
            & (omim_gene['txStart'] <= row['end']),
            'name2'
        ]),
        axis=1
    )
    ts_region['start'] = ts_region['start'].astype(int)
    ts_region['end'] = ts_region['end'].astype(int)
    ts_region.rename(columns={
        'chrom': '#chrom',
        '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
    }).to_csv(ts_region_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {ts_region_file} > {ts_region_file}.gz
        {tabix} -fp bed {ts_region_file}.gz
        rm {ts_region_file}
    """
    os.system(cmd)

    uts_region = region_pos.loc[
        region_pos['Triplosensitivity Score'] == '40',
        ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
    ].sort_values(['chrom', 'start', 'end'])
    uts_region['genes'] = uts_region.apply(
        lambda row: ','.join(gene.loc[
            (gene['chrom'] == row['chrom'])
            & (gene['txEnd'] >= row['start'])
            & (gene['txStart'] <= row['end']),
            'name2'
        ]),
        axis=1
    )
    uts_region['start'] = uts_region['start'].astype(int)
    uts_region['end'] = uts_region['end'].astype(int)
    uts_region.rename(columns={
        'chrom': '#chrom',
        '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
    }).to_csv(uts_region_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {uts_region_file} > {uts_region_file}.gz
        {tabix} -fp bed {uts_region_file}.gz
        rm {uts_region_file}
    """
    os.system(cmd)

    dgv = pd.read_csv(
        dgv_ori_file,
        sep='\t', names=['chrom', 'info'], usecols=[0, 8]
    ).drop_duplicates('info')
    info = dgv['info'].str.extract(
        r'ID=(?P<id>[^;]+).*variant_sub_type=(?P<type>[^;]+).*inner_start=(?P<start>[^;]+).*inner_end=(?P<end>[^;]+).*Frequency=(?P<freq>\S+?)%;.*num_unique_samples_tested=(?P<sample>[^;]+)'
    )
    
    # 先将数值列转换为float类型，这样可以保留NaN值
    numeric_cols = ['start', 'end', 'sample']
    info[numeric_cols] = info[numeric_cols].apply(pd.to_numeric, errors='coerce')
    
    # 去除包含NaN的行
    info = info.dropna(subset=numeric_cols)
    
    # 然后将需要的列转换为整数类型
    info[numeric_cols] = info[numeric_cols].astype(int)
    info['freq'] = pd.to_numeric(info['freq'], errors='coerce')
    dgv = dgv.merge(info, left_index=True, right_index=True)
    dgv['af'] = dgv['freq'] / 100
    dgv = dgv[dgv['sample'] >= 1000].sort_values(['chrom', 'start', 'end'])
    # def fetch_gene(gene, chrom, start, end):
    #     return ','.join(gene.loc[
    #         (gene['chrom'] == chrom) & (gene['txEnd'] >= start) & (gene['txStart'] <= end), 'name2'
    #     ])
    with Pool(processes=7) as pool:
        dgv['genes'] = pool.starmap(fetch_gene, (
            (gene, row['chrom'], row['start'], row['end']) for _, row in dgv.iterrows()
        ), chunksize=70)
    dgv.loc[
        dgv['type'] == 'Gain', ['chrom', 'start', 'end', 'id', 'genes', 'af']
    ].rename(columns={'chrom': '#chrom'}).to_csv(dgv_gain_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {dgv_gain_file} > {dgv_gain_file}.gz
        {tabix} -fp bed {dgv_gain_file}.gz
        rm {dgv_gain_file}
    """
    os.system(cmd)

    dgv.loc[
        dgv['type'] == 'Loss', ['chrom', 'start', 'end', 'id', 'genes', 'af']
    ].rename(columns={'chrom': '#chrom'}).to_csv(dgv_loss_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {dgv_loss_file} > {dgv_loss_file}.gz
        {tabix} -fp bed {dgv_loss_file}.gz
        rm {dgv_loss_file}
    """
    os.system(cmd)

    gnomad = pd.read_csv(
        gnomad_control_ori_file, sep='\t',  dtype=str,
        usecols=[0, 1, 2, 3, 4, 37, 38, 73, 74, 107, 108, 141, 142, 175, 176, 241]
    )
    gnomad = gnomad[
        (gnomad['FILTER'] == 'PASS') & gnomad['svtype'].isin(['DEL', 'DUP'])
    ]
    gnomad = gnomad[
        (gnomad['N_BI_GENOS'].astype(int) >= 1000) |
        (gnomad['AFR_N_BI_GENOS'].astype(int) >= 1000) |
        (gnomad['AMR_N_BI_GENOS'].astype(int) >= 1000) |
        (gnomad['EAS_N_BI_GENOS'].astype(int) >= 1000) |
        (gnomad['EUR_N_BI_GENOS'].astype(int) >= 1000)
    ]
    gnomad['#chrom'] = 'chr' + gnomad['#chrom']
    gnomad['start'] = gnomad['start'].astype(int)
    gnomad['end'] = gnomad['end'].astype(int)

    with Pool(processes=8) as pool:
        gnomad['genes'] = pool.starmap(fetch_gene, (
            (gene, row['#chrom'], row['start'], row['end']) for _, row in gnomad.iterrows()
        ), chunksize=70)
    gnomad['gene'] = gnomad.apply(lambda row: fetch_gene(gene, row['#chrom'], row['start'], row['end']), axis=1)
    gnomad[
        gnomad['svtype'] == 'DEL'
    ][['#chrom', 'start', 'end', 'name', 'genes', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF']].rename(columns={
        'AF': 'af', 'AFR_AF': 'af_afr', 'AMR_AF': 'af_amr', 'EAS_AF': 'af_eas', 'EUR_AF': 'af_eur'
    }).to_csv(gnomad_del_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {gnomad_del_file} > {gnomad_del_file}.gz
        {tabix} -fp bed {gnomad_del_file}.gz
        rm {gnomad_del_file}
    """
    os.system(cmd)

    gnomad.loc[
        gnomad['svtype'] == 'DUP',
        ['#chrom', 'start', 'end', 'name', 'genes', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF']
    ].rename(columns={
        'AF': 'af', 'AFR_AF': 'af_afr', 'AMR_AF': 'af_amr', 'EAS_AF': 'af_eas', 'EUR_AF': 'af_eur'
    }).to_csv(gnomad_dup_file, sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf {gnomad_dup_file} > {gnomad_dup_file}.gz
        {tabix} -fp bed {gnomad_dup_file}.gz
        rm {gnomad_dup_file}
    """
    os.system(cmd)

    refGene_df = read_gtf(exon_ori_file)
    refGene_df = refGene_df[refGene_df['feature'] == 'exon']
    gene_db = pd.read_csv('data/gene.sorted.bed.gz', sep='\t')
    merge_df = gene_db.merge(refGene_df, left_on='transcript', right_on='transcript_id')[[
        '#chrom', 'start_y', 'end_y', 'gene_id_x', 'symbol', 'exon_number', 'transcript'
    ]]
    merge_df.columns = [x.strip('_x').strip('_y') for x in merge_df.columns]
    merge_df = merge_df.sort_values(['#chrom', 'start', 'end'])
    merge_df.to_csv('data/exon.sorted.bed', sep='\t', index=False)
    cmd = f"""
        {bgzip} -cf data/exon.sorted.bed > data/exon.sorted.bed.gz
        {tabix} -fp bed data/exon.sorted.bed.gz
        rm  data/exon.sorted.bed
    """
    os.system(cmd)

############
cytoband()
gene_info()

