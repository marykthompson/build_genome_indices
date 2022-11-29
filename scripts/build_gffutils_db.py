'''
build_gffutils_db.py
- Create a gffutils database to use in downstream scripts and figure creation
- Output a file containing all rRNA type genes in order to allow their removal if needed
'''

import gffutils
import pandas as pd

gtf_file = snakemake.input['combo_gtf']
db_out = snakemake.output['gffutils_db']
rrna_gene_file = snakemake.output['rrna_gene_file']
rrna_gene_file_nopg = snakemake.output['rrna_gene_file_nopg']
len_by_gene_med_out = snakemake.output['len_by_gene_med']
len_by_gene_mean_out = snakemake.output['len_by_gene_mean']
len_by_txt_out = snakemake.output['len_by_txt']

#write gffutils db
db = gffutils.create_db(gtf_file, db_out, disable_infer_genes = True, disable_infer_transcripts = True,
force = True, merge_strategy = 'create_unique')
introns = list(db.create_introns())
#Add introns back to the database
db.update(introns, disable_infer_genes = True, disable_infer_transcripts = True,
          merge_strategy = 'create_unique', make_backup = False)

# Write lengths of intronic and exonic regions by gene and by transcript
all_genes = db.features_of_type('gene')
d = {}
for gene in all_genes:
    this_gene = gene.id
    d[this_gene] = {}
    txts = db.children(this_gene, featuretype='transcript')
    for t in txts:
        d[this_gene][t.id] = {'intron_bp':0, 'exon_bp':0}
        introns = [i for i in db.children(t.id, featuretype='intron')]
        exons = [i for i in db.children(t.id, featuretype='exon')]
        for i in exons:
            d[this_gene][t.id]['exon_bp'] += i.stop - i.start + 1
        for i in introns:
            d[this_gene][t.id]['intron_bp'] += i.stop - i.start + 1

# https://stackoverflow.com/questions/13575090/construct-pandas-dataframe-from-items-in-nested-dictionary
len_df = pd.DataFrame.from_dict({(i,j): d[i][j]
                         for i in d.keys()
                         for j in d[i].keys()}, orient='index')
len_df.index.set_names(['gene', 'txt'], inplace=True)

intron_lens_med = len_df.groupby('gene')['intron_bp'].median()
exon_lens_med = len_df.groupby('gene')['exon_bp'].median()
intron_lens_mean = len_df.groupby('gene')['intron_bp'].mean()
exon_lens_mean = len_df.groupby('gene')['exon_bp'].mean()
intron_lens_med.name = 'intron_bp'
intron_lens_mean.name = 'intron_bp'
exon_lens_med.name = 'exon_bp'
exon_lens_mean.name = 'exon_bp'
len_by_gene_med = pd.concat([intron_lens_med, exon_lens_med], axis=1)
len_by_gene_mean = pd.concat([intron_lens_mean, exon_lens_mean], axis=1)
len_by_txt = len_df.reset_index('gene', drop=True)
len_by_gene_med.to_csv(len_by_gene_med_out)
len_by_gene_mean.to_csv(len_by_gene_mean_out)
len_by_txt.to_csv(len_by_txt_out)

#write file containing all the rRNA gene IDs
allgenes = db.features_of_type('gene')
rRNA_genes_all = set()
#and also a file with no rRNA pseudogenes
rRNA_genes_nopg = set()
for gene in allgenes:
    try:
        biotype = db[gene].attributes['gene_biotype'][0]
        genename = db[gene].attributes['gene_name'][0]
        if biotype == 'rRNA':
            rRNA_genes_all.add(gene.id)
            rRNA_genes_nopg.add(gene.id)
        if 'rRNA' in genename:
            rRNA_genes_all.add(gene.id)
    except KeyError:
        continue

pd.Series(list(rRNA_genes_all)).to_csv(rrna_gene_file, index=False, header=None)
pd.Series(list(rRNA_genes_nopg)).to_csv(rrna_gene_file_nopg, index=False, header=None)
