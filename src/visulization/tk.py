"""
utilities that help auto-mated generate DNA block plot

"""



fullcurated_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_curated_status.tsv',sep='\t',index_col=0).T

kos = ['K00845', 'K01810', 'K06859', 'K16370', 'K01623', 'K01624', 'K11645',
       'K01803', 'K00134', 'K00927', 'K15633', 'K01689', 'K00873']
gids = ['F6', 'B4', 'FI4', 'GJ9', 'BA3', 'AD3', 'L7', 'AG8', 'GNM000011965']
sublocus_df = kegg_df.loc[gids,kos]
substats_df = fullcurated_df.loc[gids,kos]

for ko,col in substats_df.iteritems():
    for gid,s in col.to_dict().items():
        if 'intact' not in s:
            print(ko,gid)

