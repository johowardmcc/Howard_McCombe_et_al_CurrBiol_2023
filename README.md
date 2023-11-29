# Howard_McCombe_et_al_CurrBiol_2023

Scripts for assembly, genotyping and analyses described in Howard-McCombe et al. 2023 (https://www.cell.com/current-biology/fulltext/S0960-9822(23)01424-0)

Note: The dataset for this study was produced as part of a wider collaboration (see also: Jamieson et al. 2023 - https://www.cell.com/current-biology/fulltext/S0960-9822(23)01073-4).  As such, data for the following samples are available under the NCBI BioSample accessions:

  * FBI4 - SAMEA112127400
  * FLI3 - SAMN20192906 
  * FLI4 - SAMN20192907
  * FMA8 - SAMEA112127401
  * FSI204 - SAMEA112127402
  * FSX360 - SAMEA112127406
  * AJ324 - SAMEA112127580
  * AJ419 - SAMEA112127581

### Genome assembly
Scripts 1-8 detail the pipeline for QC, mapping and genotyping of sampels. These scripts were written for a PBS job scheduler, most are set up to run as array jobs.  Information about the tools needed to run each script are given in the 'Load programs' section.

### MOSAIC analyses
Setup for MOSAIC analyses is completed by running mosiac_input.sh and mosaic_rates_file.sh.  The MOSAIC run is controlled by run_mosaic.sh

### Masking domestic introgression
Domestic ancestry is masked in the VCF file as detailed in masking_domstic_haplotypes.R (using output from MOSAIC for local ancestry estimates).  Haploid masked data (e.g., for GONE analyses) is generated with make_haploid_data.R

### Estimating recent effective population size
Recent effective population size was estimated with GONE, run using a leave-one-out jackknife approach (leaving one chromosome out of the analyses per run) (run_gone.sh).  Jackknife bias can be estimated by comparing this output to that generated using all the data.  For this study GONE was run using the domestic and wildcat ancestry separtely, and using the complete (admixed) data.  Results from this analysis were presented as detailed in gone_plot.R

### Tests for selection
Code for selection analyses can be accessed at: https://github.com/danjlawson/localselection
