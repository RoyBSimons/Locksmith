# Locksmith Docs

## Installing Locksmith
Clone the repository with the following command:
```
$ git clone https://github.com/RoySimons96/Locksmith.git
```
Enter the repository by `cd Locksmith`.

Create a new conda environment.
```
conda create --name Locksmith
```
```
conda activate Locksmith
```
Install the following packages: biopython=1.79 .
```
conda install -c conda-forge biopython=1.79
```

### Preprocessing
Before running Locksmith several files are needed:
-	Config file

This file is supplied in the Locksmith directory and needs to be completed with the correct local paths. For the other parameters a default value is in place.

-	Chr2acc file: 

A tab delimited file which shows the chromosome number and the corresponding accession number. The header consists of two columns. 
Header: #Chromosome	Accession.version
This file can be obtained from  the NCBI FTP server e.g. https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc

-	Full genome.fa with right headers (create with generalize_genome_headers.py)

Obtain the Fasta sequence of the full reference genome.

Run the following command in order to change the headers into the chromosome numbers:
```
python scripts/preprocessing/generalize_genome_headers.py –g <Path to genome directory> -i <Path to chr2acc file>
```
-	Split genome: chr_1.fa, chr_2.fa etc. (needed to create CpG bed files with Create_CpG_IDs_from_Reference_genome.py). 
	
The files must be end with the following name format: ‘_[nr/X/Y].[fasta extension]’

Create CpG bed files for all chromosomes with the following command:
```
python scripts/preprocessing/generalize_genome_headers.py –g <Path to chromosome fasta files with wildcard> -c <Path to chr2acc file>
```
-	Target_list.bed

The target list is a bedfile with the following format:
```
Chromosome  Start End CpG_ID  Target strand
chr1  1232656 1232657 CpG23295  -
```

## Dependencies
- Python 3.9.0
- Snakemake v4.3.1
- Bedtools v2.29.1
- Primer3 v2.5.0
- mfeprimer V3.0
- tabix (htslib) 1.10.2

### Python packages
biopython-1.79
## Running Locksmith
Go into the Locksmith directory by:
```
cd <Path to Locksmith>
```
Make sure to activate the conda environment by: 
```
conda activate Locksmith
```
```
snakemake --cores [amount of cores] --configfile [path to config file]
```
### Test run
Check whether the local paths of the Locksmith/test_data/config_test.json file are all correct.
In order to perform a testrun of Locksmith go into the Locksmith directory and execute the following command:
```
snakemake --cores [amount of cores] --configfile [path to Locksmith directory]/test_data/config_test.json
```

### Configuration
The configuration file which is needed to run Locksmith includes all initial parameters needed to create probes and select a probe panel.

"path_to_config_file": The path in which the configuration file can be found.

"path_to_scripts": The directory in which the Lockmith python scripts are situated.

"target_list_bedfile_path": The path to the bedfile which states the targets for probe design.

"acc_2_chrom_list_path": The path to the chr2acc file from the used reference genome. This file is needed in order to find frequent SNPs in the probe arms through the NCBI FTP server, without needing to download all of the SNP database.

"output_directory":  The path of the directory in which the output files will be created.

"cpg_bed_path": The path (including wildcard) to the bedfiles of all CpG loci per chromosome.

"PATH_to_reference_genome_fasta": The Path to the reference genome Fasta file.

"ftp_path_snp_database": Path to the FTP server of the NCBI SNP database e.g. "ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz",

"target_range": Integer value of amount of basepairs to include on either side of the target CpG for probe design.

"snp_frequency_threshold": The threshold on the frequency of a SNP in order to count the SNP as frequent. e.g. 0.01

"dimer_exclusion_factor": The factor which determines how much dimers from the panel will be regarded as conflicting each iteration. This is a float between 0.0 and 1.0. To indicate that only one dimer must be chosen each iteration use the value 2.

"seed": The seed to use in random selection of the probes based on their cost. The seed makes sure that a Locksmith run with the exact same configuration yields the same results.

"backbone_sequence": Backbone sequence used for each probe in the panel. This consists of three parts, which are all 5’-3’ and will be appended to each other in the following order:
1.	"reverse_complement_universal_forward_primer": Sequence for annealing the forward primer during amplification of captured target.
2.	"common_region": Spacer sequence.
3.	"universal_reverse_primer": Sequence for annealing the reverse primer during amplification of captured target.
    "probe_specifics": 

"probe_specifics": A list of parameters which are used for probe design.
- "min_arm_length": Minimum length of the upstream and downstream arms.
- "max_arm_length": Maximum length of the upstream and downstream arms.
- "min_target_length": Minimum length of the captured target sequence.
- "max_target_length":  Maximum length of the captured target sequence
- "min_cg_percentage": Minimum CG percentage of the probe arms.
- "max_cg_percentage": Maximum CG percentage of the probe arms.
- "cpg_flanks": Amount of basepairs which must be flanking the target locus on both sides before annealing sites for the probe arms.
- "max_cpgs_in_arms": Maximum combined amount of CpGs in the probe arms which will be tolerated during probe design.

"mfeprimer_hairpin_parameters": Parameters to use for hairpin check.
- "score_cutoff": Hairpin score threshold for reporting hairpins.
- "tm_cutoff": Hairpin Tm threshold for reporting hairpins.

"mfeprimer_dimer_parameters": Parameters to use for dimer check.
- "score_cutoff": Dimer score threshold for reporting dimers.
- "tm_cutoff": Dimer Tm threshold for reporting dimers.

"permutations": Maximum amount of iterations to be performed before choosing the panel with the least amount of dimers.

"max_threads": Maximum amount of threads that will be used by Locksmith
