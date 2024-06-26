
# Application specification: SARS2Assembly

This is the application specification for service with identifier SARS2Assembly.

The backend script implementing the application is [App-SARS2Assembly.pl](../service-scripts/App-SARS2Assembly.pl).

The raw JSON file for this specification is [SARS2Assembly.json](SARS2Assembly.json).

This service performs the following task:   Assemble SARS2 reads into a consensus sequence

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_ids | SRR ID | string  |  |  |
| recipe | Assembly recipe | enum  |  | auto |
| primers | Primer set to use for assembly | enum  | :heavy_check_mark: | ARTIC |
| primer_version | Version number for primer | string  |  |  |
| min_depth | Minimum depth to add to consensus | int  |  | 100 |
| max_depth | Maximum depth to add to consensus | int  |  | 8000 |
| keep_intermediates | Keep all intermediate output from the pipeline | int  |  | 0 |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| debug_level | Debug level | int  |  | 0 |

