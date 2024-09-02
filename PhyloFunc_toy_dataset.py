### This code is for generating PhyloFunc distace for a toy dataset

### About original data files in the folder \1_database_search_results\1_toy_dataset:
    ## tree of toy_dataset.nwk, Taxon-Function Table.csv are written directly.
    #tree_branches,tax_composition,weighted_function_composition.csv,weighted_function_composition_percentage.csv, extend_taxon_composition_merge_all_nodes.csv are automatically generated during program execution, making it convenient to verify the calculation process of the program
from Bio import Phylo
import pandas as pd
import csv
#step1 read tree file and creat a table of branches
#read tree file
tree = Phylo.read("tree of toy_dataset.nwk", "newick")
# name internal nodes because the tree only has the name of leaves nodes
def assign_names(tree):
    node_count = 0
    for clade in tree.find_clades(order='postorder'):  # traverse from leaves to root node
        if not clade.name:  # generate a unique name of node if it doesn't have name
            node_count += 1
            clade.name = f"Node{node_count}"
    return tree
tree_with_names = assign_names(tree)
        
#Define a function to recursively traverse each branch of the tree and generate the names of the antecedent and consequent nodes for that branch
def generate_branches(clade):
    branches = []
    for child in clade.clades:
        if child.branch_length is not None:
            if clade.name:
                branches.append(clade.name)  
            else:
                branches.append(clade.clades[0].name)  
            if child.name:
                branches.append(child.name)  
            else:
                branches.append(child.clades[0].name)  
            branches.append(child.count_terminals())  
            branches.append(child.branch_length) 
    return branches
# creat all of branches and write to tree_branches.csv
def collect_branches(clade):
    if clade.is_terminal():
        return 
    branches = generate_branches(clade)
    all_branches.append(branches)
    for child_clade in clade.clades:
        collect_branches(child_clade)
all_branches = []
collect_branches(tree_with_names.root)  
with open("tree_branches.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Precedent", "consequent", "num_childnode", "branch_length"])
    for branches in all_branches:
        writer.writerow(branches[0:4])
        writer.writerow(branches[4:])
branch=pd.read_csv("tree_branches.csv", sep = ',')
#Step2 calculate tax_composition, function composition and weighted_function_composition according to Taxon-Function Table
#read Taxon-Function Table.csv
example_data=pd.read_csv("Taxon-Function Table.csv", sep = ',')
#calculate tax_composition
column_A = 'Taxon'  
column_Bs = example_data.columns.tolist()[2:]
tax_composition = pd.DataFrame()  
for column_B in column_Bs:
   grouped_sum = example_data.groupby(column_A)[column_B].sum()
   total_sum = example_data[column_B].sum()
   column_B_percentage = grouped_sum / total_sum
   tax_composition[column_B] = column_B_percentage
pd.DataFrame(tax_composition).to_csv('tax_composition.csv', sep=',', encoding='utf-8', index=True)  

#calculate function_composition
function_composition = pd.DataFrame() 
for column_B in column_Bs:
    grouped_sum = example_data.groupby(['Taxon', 'Function'])[column_B].sum()
    total_sum = example_data.groupby(column_A)[column_B].sum()
    column_B_percentage = grouped_sum / total_sum
    function_composition[column_B] = column_B_percentage
#pd.DataFrame(function_composition).to_csv('function_composition.csv', sep=',', encoding='utf-8', index=True)  
#calculate weighted_function_composition
weighted_function_composition = pd.DataFrame() 
for column_B in column_Bs:
   grouped = example_data.groupby(['Taxon', 'Function'])[column_B].sum()
   total_sum = example_data[column_B].sum()
   weighted_function_composition[column_B] = grouped / total_sum
pd.DataFrame(weighted_function_composition).to_csv('weighted_function_composition.csv', sep=',', encoding='utf-8', index=True) 

#Step3 expand to represent all nodes up to the root of the phylogeny through summing up each node,
#calculate extend_weighted_function_composition and weighted_function_composition_merge_all_nodes
#change weighted function composition to percentage and creat table weighted_function_composition_percentage.csv
weighted_function_composition=pd.read_csv("weighted_function_composition.csv", sep = ',')
tree=tree_with_names
#define a function to merge all of internal nodes
def merge_weighted_function_composition_for_inner_nodes(clade):
    leaf_nodes = pd.DataFrame()
    if not clade.is_terminal():  
        for sub_clade in clade.clades:
            if sub_clade.is_terminal() :  
                inner_node_results = pd.DataFrame(weighted_function_composition[weighted_function_composition['Taxon'] == sub_clade.name][:])
                inner_node_results['Taxon'] = clade.name  
                leaf_nodes = pd.concat([leaf_nodes, inner_node_results])
            else:
                leaf_nodes = pd.concat([leaf_nodes, merge_weighted_function_composition_for_inner_nodes(sub_clade)])   
        if not leaf_nodes.empty:
            leaf_nodes.columns = weighted_function_composition.columns
    return leaf_nodes

result_leaf_nodes = merge_weighted_function_composition_for_inner_nodes(tree.root) 
column_Bs = example_data.columns.tolist()[2:]
extend_weighted_function_composition=pd.DataFrame()

for clade in tree.find_clades():
    if not clade.is_terminal():
        leaf_set = merge_weighted_function_composition_for_inner_nodes(clade)
        if not leaf_set.empty: 
            group=pd.DataFrame()
            for column_B in column_Bs:
               group[column_B] = leaf_set.groupby(['Function'])[column_B].sum()
            group.insert(0, 'Taxon', clade.name)
            extend_weighted_function_composition = pd.concat([extend_weighted_function_composition, group])

extend_weighted_function_composition.reset_index(inplace=True)
cols = extend_weighted_function_composition.columns.tolist()
cols = cols[1:2] + cols[0:1] + cols[2:]
extend_weighted_function_composition = extend_weighted_function_composition[cols]           
weighted_function_composition_all=pd.concat([weighted_function_composition, extend_weighted_function_composition])
#pd.DataFrame(extend_weighted_function_composition).to_csv('extend_weighted_function_composition.csv', sep=',', encoding='utf-8', index=True) 
#pd.DataFrame(weighted_function_composition_all).to_csv('weighted_function_composition_merge_all_nodes.csv', sep=',', encoding='utf-8', index=True) 

#change weighted_function_composition to percentage
weighted_function_composition_percentage = pd.DataFrame()  
for column_B in column_Bs:
    grouped_sum = weighted_function_composition_all.groupby(['Taxon', 'Function'])[column_B].sum()
    total_sum = weighted_function_composition_all.groupby(column_A)[column_B].sum()
    column_B_percentage = grouped_sum / total_sum
    weighted_function_composition_percentage[column_B] = column_B_percentage
pd.DataFrame(weighted_function_composition_percentage).to_csv('weighted_function_composition_percentage.csv', sep=',', encoding='utf-8', index=True)  

#Sept4. expand to represent all nodes up to the root of the phylogeny through summing up each node,
#calculate extend_taxon_composition and extend_taxon_composition_merge_all_nodes
tree=tree_with_names
taxon_composition=pd.read_csv("tax_composition.csv", sep = ',')
def merge_taxon_composition_for_inner_nodes(clade):
    leaf_nodes = pd.DataFrame()
    if not clade.is_terminal():
        for sub_clade in clade.clades:
            if sub_clade.is_terminal():  # if it is the leaf node, add to set
                inner_node_results = pd.DataFrame(taxon_composition[taxon_composition['Taxon'] == sub_clade.name][:])
                leaf_nodes = pd.concat([leaf_nodes,inner_node_results])
            else:
                leaf_nodes = pd.concat([leaf_nodes, merge_taxon_composition_for_inner_nodes(sub_clade)])
        if not leaf_nodes.empty:  # if leaf_nodes is not null, set the name of column
            leaf_nodes.columns = taxon_composition.columns
            leaf_nodes['Taxon'] = clade.name  
    return leaf_nodes

#generate extend_taxon_composition_merge_all_nodes
result_leaf_nodes = merge_taxon_composition_for_inner_nodes(tree.root) 
extend_taxon_composition=pd.DataFrame()
for clade in tree.find_clades():
    if not clade.is_terminal():
        leaf_set = merge_taxon_composition_for_inner_nodes(clade)
        if not leaf_set.empty: 
            group = leaf_set.sum().to_frame().transpose()
            group["Taxon"]=clade.name
            extend_taxon_composition = pd.concat([extend_taxon_composition, group], ignore_index=True)
#extend_taxon_composition.to_csv('extend_taxon_composition.csv', sep=',', encoding='utf-8', index=True) 
#taxon_composition=pd.read_csv("tax_composition.csv", sep = ',')
extend_taxon_composition_merge_all_nodes=pd.concat([taxon_composition, extend_taxon_composition])
pd.DataFrame(extend_taxon_composition_merge_all_nodes).to_csv('extend_taxon_composition_merge_all_nodes.csv', sep=',', encoding='utf-8', index=True) 

# Step5 calculate PhyloFunc distance basd on tables "weighted_function_composition_percentage" and "extend_taxon_composition_merge_all_nodes"
weighted_function_composition_percentage=pd.read_csv('weighted_function_composition_percentage.csv', sep=',')  
weight=pd.read_table('tree_branches.csv', sep = ',')
dab_matrix_norm = pd.DataFrame()
taxon_len=len(weighted_function_composition_percentage.index)
column_len=len(weighted_function_composition_percentage.columns)
for a in range(2, column_len):
    Ga_function = weighted_function_composition_percentage.iloc[:, [0,1,a]]
    Ga_taxon = extend_taxon_composition_merge_all_nodes.iloc[:, [0,a-1]]
    for b in range(2, column_len):
        Gb_function = weighted_function_composition_percentage.iloc[:, b]
        Gb_taxon = extend_taxon_composition_merge_all_nodes.iloc[:, b-1]
        Gab_function = pd.concat([Ga_function, Gb_function], axis=1)
        Gab_taxon = pd.concat([Ga_taxon, Gb_taxon], axis=1)
        Phylofunc=0
        first_column_value = weighted_function_composition_percentage['Taxon'].unique()
        for t in first_column_value:
            if t=='Node2':# for the root node
                weight_taxon=1
            else:
                weight_taxon=weight.loc[weight['consequent']==t]['branch_length'].iloc[0]
            data_cog_tax=Gab_function[(Gab_function["Taxon"]==t)]
            origin_data_norm = data_cog_tax.iloc[:,2:].apply(lambda x: x/sum(x), axis = 0).fillna(0) # normalized COG within each sample
            dab = 1 - sum(origin_data_norm.apply(lambda x: min(x), axis=1)) / sum(origin_data_norm.apply(lambda x: max(x), axis=1))
            a_abundance=Gab_taxon[(Gab_taxon["Taxon"]==t)].iloc[0, 1] 
            b_abundance=Gab_taxon[(Gab_taxon["Taxon"]==t)].iloc[0, 2] 
            Phylofunc=Phylofunc+dab*weight_taxon*a_abundance*b_abundance
        dab_matrix_norm.at[Gab_function.columns.values[2], Gab_function.columns.values[3]] = Phylofunc
        dab_matrix_norm.to_csv('PhyloFunc_distance.csv')
