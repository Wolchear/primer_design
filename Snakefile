import pandas as pd
from snakemake.utils import min_version
min_version("6.0")

from workflow.lib.utils import get_path

configfile: "config/config.yaml"

samples = pd.read_csv("config/samples.csv")
TREES = set(samples["tree"].tolist())


SUBUNITS_DATA = expand(
    "{raw_fata}/{subunit}_eukaryota.fasta.gz",
    raw_fata = get_path(config['data'], 'raw'),
    subunit = ["SSU", "LSU"]
)

rule all:
    input:
        SUBUNITS_DATA
        
RULES_DIR = get_path(config['workflow'], "rules")


module raw_data:
    snakefile: f"{RULES_DIR}/raw_data.smk"
    config: config
use rule * from raw_data