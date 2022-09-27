from src.utils.yaml_utils import load_yaml
from src.utils.configuration import snakemake_validate_config, TECHNOLOGIES_CONFIG_FN, CONFIG_FN, get_number_of_indices
import os


PHASEDANCER_DATA_DIR = os.environ['PHASEDANCER_DATA_DIR']
OUTPUT_DIR = PHASEDANCER_DATA_DIR + '/'

config_technologies = load_yaml(TECHNOLOGIES_CONFIG_FN)

#config_fn = CONFIG_FN

#errors = snakemake_validate_config(config_fn, check_files = False)
errors = []

configfile: OUTPUT_DIR + 'config.yaml'

def aggregate_input(wildcards):

    checkpoint_output = checkpoints.fastq2fasta.get(**wildcards).output[0]

    return expand(OUTPUT_DIR + 'data/{sample}/index/{sample}_{number}.{ext}', **wildcards, number=range(get_number_of_indices(wildcards, config)))

if not errors:

    rule main:
        input:
            expand(OUTPUT_DIR +  'data/{sample}/{sample}.{ext}.done', sample = config['sample'], ext = ['mmi', 'fasta.fai'])

    checkpoint fastq2fasta:
        input:
            OUTPUT_DIR + 'data/{sample}/{sample}.fastq'
        output:
            clusters=directory(OUTPUT_DIR + 'data/{sample}/index')
        params:
            number_of_indices = lambda wildcards: get_number_of_indices(wildcards, config),
            output_file_pattern = lambda wildcards: OUTPUT_DIR + 'data/{sample}/index/{sample}_%d.fasta'.format(**wildcards)
        shell:
            """ mkdir -p {output} ; echo {params.number_of_indices} ;awk '{{ if (NR%4 == 1) {{ print ">" substr($0,2) }} else if ( NR%4 == 2)  {{ print $0 }} }}' {input} | """
	    """ awk 'BEGIN {{ n_seq=0; }} """
	    """ {{ if (NR%2 == 1) """
	    """        {{ file=sprintf("{params.output_file_pattern}",n_seq%({params.number_of_indices})); print ">" substr($0,2)  >> file; n_seq++;  }}  """
	    """    else  """
	    """        {{ print $0 >> file; }} """
	    """ }}' """

    rule generate_mmi_part:
        input:
            ref = OUTPUT_DIR + 'data/{sample}/index/{sample}_{number}.fasta'
        output:
            index = OUTPUT_DIR + 'data/{sample}/index/{sample}_{number}.mmi'
        params:
            preset = lambda wildcards: config_technologies['technology'][config['samples'][wildcards['sample']]['technology']]['minimap2-preset']
        threads:
            4
        log:
            err = 'data/{sample}/index/mmi-log/{sample}_{number}.mmi.log'
        shell:
           'minimap2 --idx-no-seq -I 450GB -t {threads} -x {params.preset} -d {output.index} {input.ref} 2> {log.err}'

    rule generate_fai_part:
        input:
            reads = OUTPUT_DIR + 'data/{sample}/index/{sample}_{number}.fasta'
        output:
            OUTPUT_DIR + 'data/{sample}/index/{sample}_{number}.fasta.fai'
        shell:
            'samtools faidx {input.reads}'

    rule aggregate:
        input:
           aggregate_input
        output:
            index = OUTPUT_DIR + 'data/{sample}/{sample}.{ext}.done'
        shell:
            'touch {output}'
