'''
Download and unzip the files.
Create combined gtf and fasta files from the genome and spike-ins.
'''

rule download_spikes:
    input:
    output:
        'indices/spike_files/{}.fa'.format(config['spike_name']),
        'indices/spike_files/{}.gtf'.format(config['spike_name'])
    params:
        sirv_zip_address = config['sirv_zip_address'],
        ercc_fasta_address = config['ercc_fasta_address'],
        spike_name = config['spike_name'],
        outdir = 'indices/spike_files',
        include_sirvs = config['include_sirvs'],
        include_erccs = config['include_erccs'],
        ercc_subset = config['ercc_subset']
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_spike_seqs.py'

rule download_genomic_seqs:
    input:
    output:
        fasta = 'indices/genome_files/{}_genome.fa'.format(config['genome_name']),
        gtf = 'indices/genome_files/{}.gtf'.format(config['genome_name'])
    params:
        fasta_address = config['genome_fasta_address'],
        gtf_address = config['genome_gtf_address'],
        outdir = 'indices/genome_files'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_ensembl_seqs.py'

rule combine_genome_and_spike:
    input:
        fastas = ['indices/genome_files/{}_genome.fa'.format(config['genome_name']),
        'indices/spike_files/{}.fa'.format(config['spike_name'])],
        gtfs = ['indices/genome_files/{}.gtf'.format(config['genome_name']),
        'indices/spike_files/{}.gtf'.format(config['spike_name'])]
    output:
        combo_fasta = 'indices/combo_files/{}.fa'.format(config['index_name']),
        combo_gtf = 'indices/combo_files/{}.gtf'.format(config['index_name'])
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/combine_genomes.py'
