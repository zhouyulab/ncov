#!/usr/bin/env python
"""
Get genes belonging to specific pathways
"""
import os
import sys
import json

def read_json(json_file):
    """
    Load in the ".json" file with information of pathway.
    Save as a dict.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data

def search_target_pathway(data, kegg_id):
    """
    Traversal the dict, seach for the kegg_id.
    Return the list of gene information list
    """
    for a in data["children"]:
        for b in a["children"]:
            for c in b["children"]:
                if kegg_id in c["name"]:
                    # some kegg_id doesn't have key named 'children'
                    try:
                        return(c["children"])
                    except:
                        return None

def get_gene_list(genelist):
    """
    From a line of gene information list, extract gene symbol.
    The info before the first ";" is "id-space-symbol"
    """
    gene_symbol_list = []
    for gene in genelist:
        id_and_gene = gene["name"].split(';')[1]
        gene_symbol = id_and_gene.split(' ')[-1]
        gene_symbol_list.append(gene_symbol)
    return gene_symbol_list

def output_gene_list(kegg_id_list, foldname, data, json_path):
    """
    Input a txt file of kegg id list and a json file of pathwy.
    Create a new directory with the same name of the txt.
    For each line fo kegg id list, output a txt file with gene symbol list.
    """
    # creat a new directory of the input kegg id list file
    out_path = os.path.join(json_path, foldname)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # read each line of kegg id list file
    for kegg_id in kegg_id_list:
        # set outfile name
        out_name = kegg_id + ".txt"
        out_full_path = os.path.join(out_path, out_name)
        # search for the given kegg id
        genelist = search_target_pathway(data, kegg_id)
        if genelist:
            # get the gene symbol list of given kegg id
            # output the txt file of gene symbol list 
            gene_symbol_list = get_gene_list(genelist)
            with open(out_full_path, 'w') as f:
                print(*gene_symbol_list, sep='\n', file=f)

def get_info_from_txt(json_file, highlight_id_list):
    """
    Input json file and kegg id list file,
    and split the file name, and get the directory name.
    """
    fold_name = os.path.splitext(os.path.basename(highlight_id_list))[0]
    json_path = os.path.dirname(json_file)    
    data = read_json(json_file)
    with open(highlight_id_list, 'r') as f:
        kegg_id_list = f.read().strip().split('\n')
    output_gene_list(kegg_id_list, fold_name, data, json_path)


if __name__=="__main__":
    json_file = r"hsa00001.json"
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: program query.txt\n")
        sys.exit(1)

    f_query = sys.argv[1]
    assert os.path.exists(f_query)
    get_info_from_txt(json_file, f_query)

