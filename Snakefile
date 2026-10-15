import pandas as pd
from snakemake.utils import min_version
min_version("6.0")

from workflow.lib.utils import get_path

configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.csv")
TREES = set(samples["tree"].tolist())


SUBUNITS_DATA = expand(
    "{data_dir}/{subunit}_eukaryota.fasta.gz",
    data_dir = get_path(config['data'], 'raw'),
    subunit = ["SSU", "LSU"]
)

SUBUNITS_TAXONOMY = expand(
    "{data_dir}/{subunit}_eukaryota.txt.gz",
    data_dir = get_path(config['data'], 'taxonomy'),
    subunit = ["SSU", "LSU"]
)


rule all:
    input:
        SUBUNITS_DATA,
        SUBUNITS_TAXONOMY
        
RULES_DIR = get_path(config['workflow'], "rules")


module raw_data:
    snakefile: f"{RULES_DIR}/raw_data.smk"
    config: config
use rule * from raw_data