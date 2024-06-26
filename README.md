# SARS-CoV-2 Genome Analysis

## Overview
The SARS-CoV-2 Genome Assembly and Annotation Service provides a streamlined **“meta-service”** that accepts raw reads and performs genome assembly, annotation, and variation analysis for SARS-CoV-2 genome reads. The figure below provides an overview of the workflows of the service.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

There are three distinct application service specifications defined here:

1. [ComprehensiveSARS2Analysis](app_specs/ComprehensiveSARS2Analysis.md): Service that provides the backend for the BV-BRC web interface; it takes reads as input.
2. [SARS2Assembly](app_specs/SARS2Assembly.md): This is curently only available `p3-submit-sars2-assembly` command-line script; it takes reads as input. SARS2Assembly only preforms assembly where as Comprehensive SARS2 Analysis preforms assembly and annotation.


The code in this module provides the BV-BRC application service wrapper scripts for the SARS-CoV-2 Genome Analysis service as well
as some backend utilities:

| Script name | Purpose |
| ----------- | ------- |
| [App-ComprehensiveSARS2Analysis.pl](service-scripts/App-ComprehensiveSARS2Analysis.pl) | App script for the [SARS-CoV-2 Genome Assembly and Annotation Service](https://www.bv-brc.org/docs/quick_references/services/sars_cov_2_assembly_annotation_service.html) |
| [App-SARS2Assembly.pl](service-scripts/App-SARS2Assembly.pl) | App script for the [SARS-CoV-2 Assembly service](https://www.bv-brc.org/docs/quick_references/services/sars_cov_2_assembly_annotation_service.html) |


## See also

* [SARS-CoV-2 Genome Assembly and Annotation Service](https://www.bv-brc.org/app/ComprehensiveSARS2Analysis)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/sars_cov_2_assembly_annotation_service.html)
* [SARS-CoV-2 Genome Assembly and Annotation Service Tutorial](https://www.bv-brc.org/docs/tutorial/sars_cov_2_assembly_annotation/sars_cov_2_assembly_annotation.html)



## References
Etherington, G.J., R.H. Ramirez-Gonzalez, and D. MacLean, bio-samtools 2: a package for analysis and visualization of sequence and alignment data with SAMtools in Ruby. Bioinformatics, 2015. 31(15): p. 2565-2567.

Langmead, B. and S. Salzberg, Langmead. 2013. Bowtie2. Nature Methods, 2013. 9: p. 357-359.

Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018. 34(18): p. 3094-3100.

Martin, M., Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 2011. 17(1): p. 10-12.

