import pandas as pd
from workflow.lib.utils import get_path

RAW_DATA = get_path(config['data'], 'raw')
TAXONOMY_DATA = get_path(config['data'], 'taxonomy')

SUBUNIT_TO_LINK = {
    "SSU": config['data']['ssu']['fasta'],
    "LSU": config['data']['lsu']['fasta']
}

SUBUNIT_TO_TAXONOMY = {
    "SSU": config['data']['ssu']['tax'],
    "LSU": config['data']['lsu']['tax']
}

rule fetch_raw:
    output:
        temp(f"{RAW_DATA}/{{subunit}}.fasta.gz")
    wildcard_constraints:
        subunit="SSU|LSU"
    threads: 1
    conda: '../envs/raw_data.yml'
    log:
        f"logs/raw_data/fetch_raw/{{subunit}}.log"
    params:
        link=lambda wildcards: SUBUNIT_TO_LINK[wildcards.subunit]
    shell:
        r"""
        wget -O {output} \
            -o {log:q} \
            {params.link}
        """

rule grep_eukaryota:
    input:
        f"{RAW_DATA}/{{subunit}}.fasta.gz"
    output:
        temp(f"{RAW_DATA}/{{subunit}}_eukaryota.fasta.gz")
    wildcard_constraints:
        subunit="SSU|LSU"
    threads: max(1, config['max_threads'] // 2)
    conda: '../envs/raw_data.yml'
    log:
        f"logs/raw_data/grep_eukaryota/{{subunit}}.log"
    shell:
        r"""
        seqkit grep \
            --by-name \
            --use-regexp \
            --pattern 'Eukaryota;' \
            -j {threads} \
            {input:q} \
            --out-file {output:q} \
            2> {log:q}
        """

rule fetch_taxonomy:
    output:
        temp(f"{TAXONOMY_DATA}/{{subunit}}.txt.gz")
    wildcard_constraints:
        subunit="SSU|LSU"
    threads: 1
    conda: '../envs/raw_data.yml'
    log:
        f"logs/raw_data/fetch_taxonomy/{{subunit}}.log"
    params:
        link=lambda wildcards: SUBUNIT_TO_TAXONOMY[wildcards.subunit]
    shell:
        r"""
        wget -O {output} \
            -o {log:q} \
            {params.link}
        """

rule grep_eukaryota_taxonomy:
    input:
        f"{TAXONOMY_DATA}/{{subunit}}.txt.gz"
    output:
        temp(f"{TAXONOMY_DATA}/{{subunit}}_eukaryota.txt.gz")
    wildcard_constraints:
        subunit="SSU|LSU"
    threads: max(1, config['max_threads'] // 2)
    conda: '../envs/raw_data.yml'
    log:
        f"logs/raw_data/grep_eukaryota_taxonomy/{{subunit}}.log"
    shell:
        r"""
        set -euo pipefail

        zcat {input:q} \
            | awk -F '\t' 'NR == 1 || $4 ~ /^Eukaryota;/' \
            | gzip -c > {output:q} \
            2> {log:q}
        """