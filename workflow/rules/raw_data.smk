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

ORGANISM_TO_QUERY = {
    "fungi": "(^|;)Fungi;",
    "trees": "^(" + "|".join(config['data']['trees']) + ") "
}

ORGANISM_TO_COLUMN = {
    "fungi": 4,
    "trees": 5
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

rule grep_organism_taxonomy_id:
    input:
        f"{TAXONOMY_DATA}/{{subunit}}_eukaryota.txt.gz"
    output:
        f"{TAXONOMY_DATA}/{{subunit}}_{{organism}}_ids.txt.gz"
    wildcard_constraints:
        subunit="SSU|LSU",
        organism="fungi|trees"
    threads: 1
    conda: '../envs/raw_data.yml'
    log:
        f"logs/raw_data/grep_organism_taxonomy_id/{{subunit}}_{{organism}}.log"
    params:
        query=lambda wc: ORGANISM_TO_QUERY[wc.organism],
        column=lambda wc: ORGANISM_TO_COLUMN[wc.organism]
    shell:
        r"""
        set -euo pipefail

        zcat {input:q} \
            | awk -F '\t' \
                -v query={params.query:q} \
                -v col={params.column} \
                'NR > 1 && $col ~ query {{
                    print $1 "." $2 "." $3
                }}' \
            | sort -u \
            | gzip -c > {output:q} \
            2> {log:q}
        """

rule extract_organism_sequences:
    input:
        fasta=f"{RAW_DATA}/{{subunit}}_eukaryota.fasta.gz",
        ids=f"{TAXONOMY_DATA}/{{subunit}}_{{organism}}_ids.txt.gz"
    output:
        f"{RAW_DATA}/{{subunit}}_{{organism}}.fasta.gz"
    wildcard_constraints:
        subunit="SSU|LSU",
        organism="fungi|trees"
    threads: max(1, config["max_threads"] // 2)
    conda:
        "../envs/raw_data.yml"
    log:
        f"logs/raw_data/extract_organism_sequences/{{subunit}}_{{organism}}.log"
    shell:
        r"""
        set -euo pipefail

        seqkit grep \
            --threads {threads} \
            --pattern-file <(gzip -dc {input.ids:q}) \
            {input.fasta:q} \
            --out-file {output:q} \
            2> {log:q}
        """