# panTools

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
+ pyvcf

### Arguments
    usage: python main.py stat -db DATABASE -rank [domain|phylum|class|order|family|genus] -name NAME -outdir OUTDIR [-h Help]
           python main.py search -db DATABASE -name NAME [-gene|-cds] [-h]
           python main.py viz -db DATABASE -name NAME -outdir OUTDIR [-outName OUTNAME] [-t Thread] [-x WIDTH] [-y HEIGTH] [-r PATH_RANGE] [-h Help]
           python main.py describe -db DATABASE -name NAME -outdir OUTDIR [-t Thread] [-h Help]
           python main.py tree [-vcf VCF] -db DATABASE -name NAME -gene GENE -outdir OUTDIR [-t Thread] [-h Help]
