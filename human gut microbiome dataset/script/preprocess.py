import pandas as pd
from collections import Counter
import os
os.chdir("D:/Unifrac-functional distance/13-preprocess/V49")
# read the protein identification table.
protein_group = pd.read_table("proteinGroups.txt", sep = "\t")
# First, filter any protein group that includes REV_ in the id names
protein_group = protein_group[~protein_group['Protein IDs'].str.contains('REV_')]
# Second, expand protein groups by semicolon
# Split the values in the "Protein Group accession" column by semicolons and then explode the dataframe
protein_group['Protein IDs'] = protein_group['Protein IDs'].str.split(';')  # Split values by semicolons
protein_group = protein_group.explode('Protein IDs')  # Explode the dataframe
protein_group_filter = protein_group.filter(regex='^(id|Protein IDs|.*ntensity.*)$')
protein_group_filter=protein_group_filter.drop("Intensity",axis=1)
protein_group_filter.columns = [col.replace(' ', '.') for col in protein_group_filter.columns]

# The original file contains protein id and taxonomic annotations
original = pd.read_table("gtdbtk.bac120.summary.tsv", sep = "\t")
original=original[["user_genome","classification"]]
# We need to match to the ids first
cog_function=pd.read_csv("cog_V49.csv", sep = ',')
#delete protein lines with nan value
cog_function = cog_function.dropna(subset=['Protein ID'])
cog_function['genome_ID'] = cog_function['Protein ID'].str.extract(r'(.*?bin.*?)_')

# First, perform taxonomic annotation. 
# generate a table containing ProteinGroupIDs, and taxonomic annotation of each protein in the protein groups.
merged_taxon_cog= cog_function.merge(original, left_on='genome_ID', right_on='user_genome', how='inner')
merged_taxon_cog = merged_taxon_cog.drop("user_genome",axis=1)

split_columns = merged_taxon_cog['classification'].str.split(';', expand=True)  # Split values by semicolons
merged_taxon_cog_split = pd.concat([merged_taxon_cog, split_columns], axis=1)
merged_taxon_cog_split.columns = merged_taxon_cog_split.columns.tolist()[:-7] + ['d','p','c','o','f', 'g', 's']

# Now we merge taxon_cog table with the protein data.
merge_proteingroup=protein_group_filter.merge(merged_taxon_cog_split, left_on='Protein.IDs', right_on='Protein ID', how='inner')
merge_proteingroup_tidy=merge_proteingroup.drop(["classification","Protein.IDs"],axis=1)
merge_proteingroup_tidy=merge_proteingroup_tidy.rename(columns={'id':'ProteinGroupIDs'})
protein_group_LFQ_intensity=merge_proteingroup_tidy.filter(regex='^(?!.*Intensity).*$', axis=1)
protein_group_Intensity=merge_proteingroup_tidy.filter(regex='^(?!.*LFQ).*$', axis=1)

#change column name of sample from table meta_data.csv
meta_data=pd.read_csv("meta_data.csv", sep = ',')
meta_data_1=meta_data[["Name","Peptide.Name"]]
meta_data_2=meta_data[["Name","Protein.Name"]]
meta_data_1.rename(columns={"Peptide.Name": "sample_name"}, inplace=True)
meta_data_2.rename(columns={"Protein.Name": "sample_name"}, inplace=True)
sample_name= pd.concat([meta_data_1, meta_data_2], ignore_index=True)


for column in protein_group_Intensity.columns:
    new_name=sample_name.loc[sample_name["sample_name"]==column,"Name"].values
    if len(new_name) > 0:
        protein_group_Intensity.rename(columns={column:new_name[0]}, inplace=True)

for column in protein_group_LFQ_intensity.columns:
    new_name=sample_name.loc[sample_name["sample_name"]==column,"Name"].values
    if len(new_name) > 0:
        protein_group_LFQ_intensity.rename(columns={column:new_name[0]}, inplace=True)
        
## Compute LCA on the species level use LFQintensity
Group_ID_species =  protein_group_LFQ_intensity[['ProteinGroupIDs', 's']].reset_index(inplace=False).drop(columns='index')
LCA_species = Group_ID_species.drop_duplicates() # remove if an id matches to same species
LCA_species_key = pd.DataFrame.from_dict(Counter(LCA_species['ProteinGroupIDs']).keys(), orient='columns')
LCA_species_number = pd.DataFrame.from_dict(Counter(LCA_species['ProteinGroupIDs']).values(), orient='columns')
LCA_species_count_list = pd.concat([LCA_species_key, LCA_species_number], axis=1)
LCA_species_count_list.columns = ['Group_ID', 'Count']
Unique_LCA_species_level = LCA_species_count_list.loc[LCA_species_count_list['Count'] == 1]
Protein_species_unique = protein_group_LFQ_intensity.merge(Unique_LCA_species_level, left_on='ProteinGroupIDs', right_on='Group_ID')
Protein_species_unique = Protein_species_unique[pd.notnull(Protein_species_unique['s'])]
Protein_species_unique=Protein_species_unique[Protein_species_unique['s']!='s__']
Protein_species_unique['LCA'] = "species"
Protein_species_unique

Protein_species_unique_index = Protein_species_unique['ProteinGroupIDs']
merge_table_noSpeciesLCA = protein_group_LFQ_intensity[~protein_group_LFQ_intensity['ProteinGroupIDs'].isin(Protein_species_unique_index)]
merge_table_noSpeciesLCA

Group_ID_genus =  merge_table_noSpeciesLCA[['ProteinGroupIDs', 'g']].reset_index(inplace=False).drop(columns='index')
LCA_genus = Group_ID_genus.drop_duplicates()  # remove if an id matches to the same genus
LCA_genus_key = pd.DataFrame.from_dict(Counter(LCA_genus['ProteinGroupIDs']).keys(), orient='columns')
LCA_genus_number = pd.DataFrame.from_dict(Counter(LCA_genus['ProteinGroupIDs']).values(), orient='columns')
LCA_genus_count_list = pd.concat([LCA_genus_key, LCA_genus_number], axis=1)
LCA_genus_count_list.columns = ['Group_ID', 'Count']
Unique_LCA_genus_level = LCA_genus_count_list.loc[LCA_genus_count_list['Count'] == 1]
Protein_genus_unique = merge_table_noSpeciesLCA.merge(Unique_LCA_genus_level, left_on='ProteinGroupIDs', right_on='Group_ID')
Protein_genus_unique = Protein_genus_unique[pd.notnull(Protein_genus_unique['g'])]
Protein_genus_unique=Protein_genus_unique[Protein_genus_unique['g']!='g__']
Protein_genus_unique['LCA'] = "genus"
Protein_genus_unique

len(Protein_genus_unique['ProteinGroupIDs'].drop_duplicates())
Protein_genus_unique_index = Protein_genus_unique['ProteinGroupIDs']
merge_table_noGenusLCA = merge_table_noSpeciesLCA[~merge_table_noSpeciesLCA['ProteinGroupIDs'].isin(Protein_genus_unique_index)]
len(merge_table_noGenusLCA)

Group_ID_family =  merge_table_noGenusLCA[['ProteinGroupIDs', 'f']].reset_index(inplace=False).drop(columns='index')
LCA_family = Group_ID_family.drop_duplicates()  # remove if an id matches to the same family
LCA_family_key = pd.DataFrame.from_dict(Counter(LCA_family['ProteinGroupIDs']).keys(), orient='columns')
LCA_family_number = pd.DataFrame.from_dict(Counter(LCA_family['ProteinGroupIDs']).values(), orient='columns')
LCA_family_count_list = pd.concat([LCA_family_key, LCA_family_number], axis=1)
LCA_family_count_list.columns = ['Group_ID', 'Count']
Unique_LCA_family_level = LCA_family_count_list.loc[LCA_family_count_list['Count'] == 1]
Protein_family_unique = merge_table_noGenusLCA.merge(Unique_LCA_family_level, left_on='ProteinGroupIDs', right_on='Group_ID')
Protein_family_unique = Protein_family_unique[pd.notnull(Protein_family_unique['f'])]
Protein_family_unique=Protein_family_unique[Protein_family_unique['f']!='f__']
Protein_family_unique['LCA'] = "family"
Protein_family_unique

len(Protein_family_unique['ProteinGroupIDs'].drop_duplicates())
Protein_family_unique_index = Protein_family_unique['ProteinGroupIDs']
merge_table_noFamilyLCA = merge_table_noGenusLCA[~merge_table_noGenusLCA['ProteinGroupIDs'].isin(Protein_family_unique_index)]
merge_table_noFamilyLCA

Group_ID_order =  merge_table_noFamilyLCA[['ProteinGroupIDs', 'o']].reset_index(inplace=False).drop(columns='index')
LCA_order = Group_ID_order.drop_duplicates()  # remove if an id matches to the same order
LCA_order_key = pd.DataFrame.from_dict(Counter(LCA_order['ProteinGroupIDs']).keys(), orient='columns')
LCA_order_number = pd.DataFrame.from_dict(Counter(LCA_order['ProteinGroupIDs']).values(), orient='columns')
LCA_order_count_list = pd.concat([LCA_order_key, LCA_order_number], axis=1)
LCA_order_count_list.columns = ['Group_ID', 'Count']
Unique_LCA_order_level = LCA_order_count_list.loc[LCA_order_count_list['Count'] == 1]
Protein_order_unique = merge_table_noFamilyLCA.merge(Unique_LCA_order_level, left_on='ProteinGroupIDs', right_on='Group_ID')
Protein_order_unique = Protein_order_unique[pd.notnull(Protein_order_unique['o'])]
Protein_order_unique=Protein_order_unique[Protein_order_unique['o']!='o__']
Protein_order_unique['LCA'] = "order"
Protein_order_unique

len(Protein_order_unique['ProteinGroupIDs'].drop_duplicates())
Protein_order_unique_index = Protein_order_unique['ProteinGroupIDs']
merge_table_noOrderLCA = merge_table_noFamilyLCA[~merge_table_noFamilyLCA['ProteinGroupIDs'].isin(Protein_order_unique_index)]
merge_table_noOrderLCA

Group_ID_class =  merge_table_noOrderLCA[['ProteinGroupIDs', 'c']].reset_index(inplace=False).drop(columns='index')
LCA_class = Group_ID_class.drop_duplicates()  # remove if an id matches to the same class
LCA_class_key = pd.DataFrame.from_dict(Counter(LCA_class['ProteinGroupIDs']).keys(), orient='columns')
LCA_class_number = pd.DataFrame.from_dict(Counter(LCA_class['ProteinGroupIDs']).values(), orient='columns')
LCA_class_count_list = pd.concat([LCA_class_key, LCA_class_number], axis=1)
LCA_class_count_list.columns = ['Group_ID', 'Count']
Unique_LCA_class_level = LCA_class_count_list.loc[LCA_class_count_list['Count'] == 1]
Protein_class_unique = merge_table_noOrderLCA.merge(Unique_LCA_class_level, left_on='ProteinGroupIDs', right_on='Group_ID')
Protein_class_unique = Protein_class_unique[pd.notnull(Protein_class_unique['c'])]
Protein_class_unique=Protein_class_unique[Protein_class_unique['c']!='c__']
Protein_class_unique['LCA'] = "class"
Protein_class_unique

len(Protein_class_unique['ProteinGroupIDs'].drop_duplicates())
Protein_class_unique_index = Protein_class_unique['ProteinGroupIDs']
merge_table_noClassLCA = merge_table_noOrderLCA[~merge_table_noOrderLCA['ProteinGroupIDs'].isin(Protein_class_unique_index)]
merge_table_noClassLCA

Group_ID_phylum =  merge_table_noClassLCA[['ProteinGroupIDs', 'p']].reset_index(inplace=False).drop(columns='index')
Group_ID_phylum['p'] = Group_ID_phylum['p'].str.split('_').str.get(2)
Group_ID_phylum['p'] = Group_ID_phylum['p'].apply(lambda x: 'p__' + str(x))
LCA_phylum = Group_ID_phylum.drop_duplicates() # remove if an id matches to same phylum
LCA_phylum_key = pd.DataFrame.from_dict(Counter(LCA_phylum['ProteinGroupIDs']).keys(), orient='columns')
LCA_phylum_number = pd.DataFrame.from_dict(Counter(LCA_phylum['ProteinGroupIDs']).values(), orient='columns')
LCA_phylum_count_list = pd.concat([LCA_phylum_key, LCA_phylum_number], axis=1)
LCA_phylum_count_list.columns = ['Group_ID', 'Count']
Unique_LCA_phylum_level = LCA_phylum_count_list.loc[LCA_phylum_count_list['Count'] == 1]
Protein_phylum_unique = merge_table_noClassLCA.merge(Unique_LCA_phylum_level, left_on='ProteinGroupIDs', right_on='Group_ID')
Protein_phylum_unique = Protein_phylum_unique[pd.notnull(Protein_phylum_unique['p'])]
Protein_phylum_unique=Protein_phylum_unique[Protein_phylum_unique['p']!='p__']
Protein_phylum_unique['LCA'] = "phylum"
Protein_phylum_unique

len(Protein_phylum_unique['ProteinGroupIDs'].drop_duplicates())
Protein_phylum_unique_index = Protein_phylum_unique['ProteinGroupIDs']
merge_table_nophylumLCA = merge_table_noClassLCA[~merge_table_noClassLCA['ProteinGroupIDs'].isin(Protein_phylum_unique_index)]
merge_table_nophylumLCA

Protein_phylum_unique_index = Protein_phylum_unique['ProteinGroupIDs']
Protein_LCA_table = pd.concat([Protein_species_unique, Protein_genus_unique, Protein_family_unique, 
                               Protein_order_unique, Protein_class_unique, Protein_phylum_unique])
Protein_LCA_table=Protein_LCA_table.drop(columns=["Protein ID"])
Protein_LCA_table=Protein_LCA_table.drop_duplicates()
Protein_LCA_table.to_csv("Protein_LCA_table_V49_LFQ.csv")
Protein_LCA_species = Protein_LCA_table[(Protein_LCA_table['LCA']=="species")]
Protein_LCA_species_tidy=Protein_LCA_species.drop(columns=["COG_category","Group_ID","Count",'d','p','c','o','f', 'g','LCA','genome_ID'])
Protein_LCA_species_tidy.to_csv("Protein_LCA_table_V49_species_LFQ.csv")
Protein_LCA_species_tidy=Protein_LCA_species_tidy.drop(columns=["ProteinGroupIDs"])
Protein_LCA_species_COG_sum = Protein_LCA_species_tidy.groupby(['s', 'COG_number']).sum().reset_index()
Protein_LCA_species_COG_sum['s'] = Protein_LCA_species_COG_sum['s'].str.replace('s__', 's_')
Protein_LCA_species_COG_sum['s'] = Protein_LCA_species_COG_sum['s'].str.replace(' ', '_')
Protein_LCA_species_COG_sum = Protein_LCA_species_COG_sum.rename(columns={"s": "taxon"})
Protein_LCA_species_COG_sum.to_csv("Protein_LCA_species_COG_sum_V49_LFQ.csv",index=False)
