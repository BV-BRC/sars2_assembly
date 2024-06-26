
# Application specification: ComprehensiveSARS2Analysis

This is the application specification for service with identifier ComprehensiveSARS2Analysis.

The backend script implementing the application is [App-ComprehensiveSARS2Analysis.pl](../service-scripts/App-ComprehensiveSARS2Analysis.pl).

The raw JSON file for this specification is [ComprehensiveSARS2Analysis.json](ComprehensiveSARS2Analysis.json).

This service performs the following task:   Analyze a genome from reads or contigs, generating a detailed analysis report.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| input_type | Input Type | enum  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_ids | SRR ID | string  |  |  |
| recipe | Assembly recipe | enum  |  | auto |
| primers | Primer set to use for assembly | enum  | :heavy_check_mark: | ARTIC |
| primer_version | Version number for primer | string  |  |  |
| min_depth | Minimum depth to add to consensus | int  |  | 100 |
| keep_intermediates | Keep all intermediate output from the pipeline | int  |  | 0 |
| genbank_file | Genbank file | WS: genbank_file  |  |  |
| contigs | Contig file | WS: Contigs  |  |  |
| scientific_name | Scientific Name | string  | :heavy_check_mark: |  |
| taxonomy_id | NCBI Taxonomy ID | int  | :heavy_check_mark: |  |
| code | Genetic Code | enum  | :heavy_check_mark: | 1 |
| domain | Domain | enum  | :heavy_check_mark: | Viruses |
| public | Public | bool  |  | 0 |
| queue_nowait | Don't wait on indexing queue | bool  |  | 0 |
| skip_indexing | Don't index genome | bool  |  | 0 |
| reference_genome_id | Reference genome ID | string  |  |  |
| reference_virus_name | Reference virus name | string  |  |  |
| container_id | (Internal) Container to use for this run | string  |  |  |
| _parent_job | (Internal) Parent job for this annotation | string  |  |  |
| workflow | Custom workflow | string  |  |  |
| analyze_quality | Enable quality analysis of genome | bool  |  |  |
| debug_level | Debug level | int  |  | 0 |

