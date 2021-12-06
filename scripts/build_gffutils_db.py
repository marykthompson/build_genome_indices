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

#write gffutils db
db = gffutils.create_db(gtf_file, db_out, disable_infer_genes = True, disable_infer_transcripts = True,
force = True, merge_strategy = 'create_unique')
introns = list(db.create_introns())
#Add introns back to the database
db.update(introns, disable_infer_genes = True, disable_infer_transcripts = True,
          merge_strategy = 'create_unique', make_backup = False)

#write file containing all the rRNA gene IDs
allgenes = db.features_of_type('gene')
rRNA_genes = []
for gene in allgenes:
    try:
        biotype = db[gene].attributes['gene_biotype'][0]
        if biotype == 'rRNA':
            rRNA_genes.append(gene.id)
    except KeyError:
        continue

pd.Series(rRNA_genes).to_csv(rrna_gene_file, index=False, header=None)
