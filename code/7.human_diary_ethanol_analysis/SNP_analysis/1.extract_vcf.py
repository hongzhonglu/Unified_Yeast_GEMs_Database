# -*- coding: utf-8 -*-
# date : 2024/11/14
'''Extract human,bioethanol,diary and wild type strains from 1011Matrix.gvcf, ouput as vcf and tsv'''
import vcfpy
import pandas as pd
import tqdm


def ids_match(key_list, id_list):
    '''check if the key_list is in id_list'''
    id_match_dict=dict()
    for key in key_list:
        value=key[:3]
        if value in id_list:
            id_match_dict[key]=value
        elif 'SACE_'+value in id_list:
            id_match_dict[key]='SACE_'+value
        else:
            print('Can not find %s' % key)
    return id_match_dict


def get_key_by_value(d, target_value):
    for k, v in d.items():
        if v == target_value:
            return k
    return None

def calculate_record_genotype_freq(record):
    record_genotype=dict()
    for call in record.calls:
        record_genotype[call.sample]=call.data['GT']
    df_genotype=pd.Series(record_genotype)
    record_genotype_freq=df_genotype.value_counts(normalize=True)
    return record_genotype_freq


reader = vcfpy.Reader.from_path('/mnt/e/data/sce1011/1011Matrix.gvcf/1011Matrix.gvcf')
na1011_list=list(reader.header.samples.names)

# load strain info
df_strain_info=pd.read_excel(r'data/1897_strains_info.xlsx', index_col=0)
wildtypelist=['14. CHNIII ','20. CHN V ', '15. CHNII ','17. Taiwanese ', '24. Asian islands ', '18. Far East Asia ', '19. Malaysian ', '22. Far East Russian ']
wt_strainList=df_strain_info[(df_strain_info['nature_clade'].isin(wildtypelist)) & (df_strain_info['type']=='Wild')].index.tolist()
bioethanol_strainList=df_strain_info[(df_strain_info['nature_clade']=='3. Brazilian bioethanol ') & (df_strain_info['type']=='Industry')].index.tolist()
human_strainList=df_strain_info[(df_strain_info['nature_clade']=='10. French Guiana human ')&(df_strain_info['type']=='Human')].index.tolist()
dairy_strainList=df_strain_info[(df_strain_info['nature_clade']=='5. French dairy ')&(df_strain_info['type']=='Fermentation')].index.tolist()
select_strainList=human_strainList+bioethanol_strainList+dairy_strainList+wt_strainList

sce1807id_to_na1011id_dict=ids_match(select_strainList,na1011_list)

# update header
new_header = reader.header.copy()
new_header.samples = vcfpy.SamplesInfos(
    [get_key_by_value(sce1807id_to_na1011id_dict,s) for s in reader.header.samples.names if s in sce1807id_to_na1011id_dict.values()]
)
# write to vcf
writer = vcfpy.Writer.from_path('/mnt/e/data/sce_paper/human_bioethanol_dairy_wt_SNP_indel_filter.gvcf', new_header)
for record in tqdm.tqdm(reader):
    # update id
    for call in record.calls:
        if call.sample in sce1807id_to_na1011id_dict.values():
            new_id=get_key_by_value(sce1807id_to_na1011id_dict,call.sample)
            call.sample=new_id
            record.call_for_sample[new_id]=call
    record.calls=[call for call in record.calls if call.sample in select_strainList]

    # filter according to genotype freq: None<0.5, 0 and 1 both > 0.1
    record_genotype_freq=calculate_record_genotype_freq(record)
    # make sure 0/0 and 1/1 are in included
    if '0/0' in record_genotype_freq.index and '1/1' in record_genotype_freq.index:
        if record_genotype_freq['0/0']>0.1 and record_genotype_freq['1/1']>0.1:
            if './.' in record_genotype_freq.index and record_genotype_freq['./.']<0.5:
                writer.write_record(record)
                print('write valuable record %s' % record.CHROM+':'+str(record.POS))
            elif './.' not in record_genotype_freq.index:
                writer.write_record(record)
                print('write valuable record %s' % record.CHROM+':'+str(record.POS))
writer.close()


# extract all SNPs and indels respectively and output as tsv
reader = vcfpy.Reader.from_path('/mnt/e/data/sce_paper/human_bioethanol_dairy_wt_SNP_indel.gvcf')
cols = ['chromosome', 'position','QUAL','type', 'ref', 'alt'] + reader.header.samples.names
df = pd.DataFrame(columns=cols)
for record in reader:
    chrosome = record.CHROM
    position = record.POS
    QUAL = record.QUAL
    type = record.ALT[0].type
    ref = record.REF
    alt = record.ALT[0].value
    line = [chrosome, position, QUAL, type, ref, alt]
    for sample in record.calls:
        gt=sample.data['GT']
        if gt =='0/0':
            gt=0
        elif gt=='1/1':
            gt=1
        else:
            gt=None
        line.append(gt)

