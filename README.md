# MGPGtools
Tools for microbial graph pangenomes

## Contents

## Installation
#### Pre-installation
+ python3
+ odgi(0.8.3)
+ seqkit(v2.5.1)
+ prettytable
+ gfapy
+ Biopython
+ clusterW
+ toytree

### Arguments
    usage: python main.py stat -db DATABASE -rank [domain|phylum|class|order|family|genus] -name NAME -outdir OUTDIR [-h Help]
       eg. python main.py stat -db path_to_database -rank order -name Lactobacillales -outdir path_to_dir

           python main.py search -db DATABASE -name NAME [-gene|-cds] [-h]
       eg. python main.py search -db path_to_database -name Listeria [-gene LSE_RS00145] [-cds WP_003721652.1]

           python main.py viz -db DATABASE -name NAME -outdir OUTDIR [-outName OUTNAME] [-t Thread] [-x WIDTH] [-y HEIGTH] [-r PATH_RANGE] [-h Help]
       eg. python main.py viz -db path_to_database -name Listeria -outdir path_to_dir -outName Listeria -t 16 -x 300 -y 200 [-r GCF_000027145#1#NC_013891.1:31414-31752]

           python main.py describe -db DATABASE -name NAME -outdir OUTDIR [-t Thread] [-h Help]
       eg. python main.py describe -db path_to_database -name Listeria -outdir path_to_dir -t 16

           python main.py tree [-gfa GFA] [-sampleTxt sampleTxt] [-label LABEL] -db DATABASE -name NAME -gene GENE -outdir OUTDIR [-t Thread] [-h Help]
       eg. python main.py tree -db path_to_database -name Listeria -gene LSE_RS00145 -outdir path_to_dir -t 16
       eg. python main.py tree -gfa test/sample.gfa -sampleTxt test/sample.txt -label sample -db path_to_database -name Listeria -gene LSE_RS00145 -outdir path_to_dir -t 16
