import pandas as pd

meta = pd.read_csv("/home/nolson/Documents/CCQM/CCQM_Analysis_taxa_id/ccqm_metadata_rdp_analysis.txt", \
                   names = ["file","org","method","lab","rep","platform","template","trim"])
ecoli_taxa_levels = ["Proteobacteria","Gammaproteobacteria","Enterobacteriales","Enterobacteriaceae","Escherichia/Shigella"]
lmono_taxa_levels = ["Firmicutes","Bacilli","Bacillales","Listeriaceae","Listeria"]

taxa_all = pd.DataFrame()

for id_file in meta['file']:
    print id_file
    file_name = "/home/nolson/Documents/CCQM/CCQM_Analysis_taxa_id/" + id_file
    taxa_file = pd.read_csv(file_name, sep = "\t", names = ["Read","taxonomy","confidence"])
    taxa_file['file'] = id_file
    taxa_split = pd.DataFrame(taxa_file.taxonomy.str.split(';').tolist(), \
                                   columns = ["Domain","Phylum","Class","Order","Family", "Genus"])
    taxa = pd.concat([taxa_file, taxa_split], axis=1)
    taxa = pd.merge(taxa, meta, on = 'file')
    
    taxonomy = taxa.taxonomy
    tax = []
    level = []
    for i in taxonomy:
        i_split = i.split(';')
        level.append(len(i_split))
        tax.append(i_split[-1])  
    taxa['tax'] = tax
    taxa['level'] = level
    taxa['id_type'] = "contam"

    
    if taxa['org'][0] == "Ecoli":
        for i in range(0,5):
            taxa['id_type'][taxa['tax'] == ecoli_taxa_levels[i]] = "expected"
        group = [taxa.groupby(['Class','level','id_type']), \
            taxa.groupby(['Order','level','id_type']), \
            taxa.groupby(['Family','level','id_type'])]
        for i in range(0,3):
            if (ecoli_taxa_levels[i + 1], i + 4, 'contam') in group[i].groups.keys():
                taxa['id_type'][group[i].groups[(ecoli_taxa_levels[i + 1], i + 4, 'contam')]] = "neighbor"
    else:
        for i in range(0,5):
            taxa['id_type'][taxa['tax'] == lmono_taxa_levels[i]] = "expected"
        group = [taxa.groupby(['Class','level','id_type']), \
            taxa.groupby(['Order','level','id_type']), \
            taxa.groupby(['Family','level','id_type'])]
        for i in range(0,3):
            if (lmono_taxa_levels[i + 1], i + 4, 'contam') in group[i].groups.keys():
                taxa['id_type'][group[i].groups[(lmono_taxa_levels[i + 1], i + 4, 'contam')]] = "neighbor"
    taxa_all = taxa_all.append(taxa)

taxa_all.to_csv("contam.csv", sep=',', na_rep='NA',header=True)

