import os
from snakemake.utils import validate, min_version

#need to use workflow.basedir to get the relative directory
snake_dir = workflow.basedir

##### set minimum snakemake version #####
min_version('5.1.2')

##### load config and sample sheets #####
configfile: 'config.yaml'

##### target rules #####
rule all:
    input:
        #'indices/star_index_{index_name}'.format(index_name = config['index_name']),
        'indices/kallisto_index/{index_name}.idx'.format(index_name = config['index_name'])

##### setup report #####
report: 'report/workflow.rst'

##### load rules #####
include: 'rules/prepare_seqs.smk'
include: 'rules/build_indices.smk'
