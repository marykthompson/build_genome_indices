rule build_star_index:
    input:
        fasta = 'indices/combo_files/{}.fa'.format(config['index_name']),
        gtf = 'indices/combo_files/{}.gtf'.format(config['index_name'])
    output:
        outdir = directory('indices/star_index_{index_name}'.format(index_name = config['index_name']))
    params:
        outdir = 'indices/star_index_{index_name}'.format(index_name = config['index_name']),
        extra='--genomeSAsparseD {} --genomeSAindexNbases {} --limitGenomeGenerateRAM {}'
        .format(config['build_star_params']['genomeSAsparseD'],
        config['build_star_params']['genomeSAindexNbases'],
        config['build_star_params']['limitGenomeGenerateRAM'])
    log:
        'indices/log_files/star_build.log'
    threads: 24
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/build_star_index.py'

rule build_kallisto_index:
    input:
        fasta = 'indices/combo_files/{}.fa'.format(config['index_name']),
        gtf = 'indices/combo_files/{}.gtf'.format(config['index_name'])
    output:
        index = 'indices/kallisto_index/{}.idx'.format(config['index_name']),
        t2g = 'indices/{}_t2g.txt'.format(config['index_name']),
        cdna_fasta = 'indices/{}_cdna.fa'.format(config['index_name']),
        intron_fasta = 'indices/{}_intron.fa'.format(config['index_name']),
        cdna_t2c = 'indices/{}_cdna_t2c.txt'.format(config['index_name']),
        intron_t2c = 'indices/{}_intron_t2c.txt'.format(config['index_name'])
    conda:
        '../envs/main.yaml'
    log:
        'indices/log_files/kallisto_build.log'
    shell:
        '''
        kb ref -i {output.index} -g {output.t2g} -f1 {output.cdna_fasta} \
        -f2 {output.intron_fasta} -c1 {output.cdna_t2c} -c2 {output.intron_t2c} \
        --lamanno {input.fasta} {input.gtf}
        '''
