# Locksmith Docs
## Installing Locksmith
Clone the repository with the following command:
```
$ git clone https://github.com/RoySimons96/Locksmith.git
```
Enter the repository by `cd Locksmith`.

## Dependencies
- [Bedtools v2.29.1](https://github.com/arq5x/bedtools2/releases/tag/v2.30.0)

Download the right version of bedtools
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
```
Rename the file
```
mv bedtools.static.binary bedtools
```
Make sure the file has the executable permission
```
chmod a+x bedtools
```
Export the executable to PATH
```
export PATH=<Path to bedtools executable>:$PATH
```
- [mfeprimer V3.2.4](https://github.com/quwubin/MFEprimer-3.0/releases/tag/v3.2.4)

Download the right version of mfeprimer
```
wget -c https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.2.4/mfeprimer-3.2.4-linux-amd64.gz
```
Uncompress the file
```
gunzip mfeprimer-3.2.4-linux-amd64.gz
```
Rename the file
```
mv mfeprimer-3.2.4-linux-amd64 mfeprimer
```
Make sure the file has the executable permission
```
chmod a+x mfeprimer
```
Export the executable to PATH
```
export PATH=<Path to mfeprimer executable>:$PATH
```
- [tabix (htslib) 1.10.2](https://github.com/samtools/htslib/releases/tag/1.10.2)

Download tabix 
```
wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2
```
Untar the tabix file
```
tar -xvf htslib-1.10.tar.bz2
```
Change working directory
```
cd htslib-1.10
```
Install tabix
```
make
```
```
make install
```
Export the executable to PATH
```
export PATH=<Path to tabix executable>:$PATH
```
### Set up conda environment
- Download [Conda](https://www.anaconda.com/products/individual) with Python 3.9

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
```

Install Conda
```
bash Anaconda3-2021.11-Linux-x86_64.sh
```

- [Snakemake](https://snakemake.readthedocs.io/) (at least v4.3.1) and [biopython-1.79](https://biopython.org/docs/1.79/api/Bio.html)

Install Mamba into your Conda-based python distribution
```
conda install -n base -c conda-forge mamba
```
Activate the Conda base environment (which now includes Mamba).
```
conda activate base
```
Create a new conda environment called ```Locksmith``` with snakemake in it.
```
mamba create -c conda-forge -c bioconda -n snakemake Locksmith
```
Activate the ```Locksmith``` conda environment.
```
conda activate Locksmith
```
Check whether Snakemake is succesfully installed by running the following command:
```
snakemake --help
```
Install the following packages: biopython=1.79 .
```
conda install -c conda-forge biopython=1.79
```
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
In order to perform a testrun of Locksmith go into the Locksmith directory and execute the following command:
```
snakemake --cores [amount of cores] --configfile test_data/config_test.json
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
This file can be created with the following command:
```
python scripts/preprocessing/cg_illumina_to_target_bedfile.py -c <Path to Illumina cg identifier file> -i <Path to Infinium annotation file> -o <output_file_name> -p <Path to directory with CpG bedfiles>*.bed -e .bed
```
The Infinium annotation file can be obtained at e.g. https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.tsv.gz

The Illumina cg identifier file is a headerless file with only the Illumina cg identifier for each target per row.

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

"scoring_weights": Weights for calculating the probe cost.
- "hairpin": Hairpin weight
- "cpg": CpG weight
- "snp": SNP weight
- "tm": Tm weight

Formula: Hairpin_score * hairpin_weight + nr_of_cpgs_in_arms * cpg_weight + nr_of_snps_in_arms * snp_weight + Tm_difference_between_arms * tm_weight

"permutations": Maximum amount of iterations to be performed before choosing the panel with the least amount of dimers.

"max_threads": Maximum amount of threads that will be used by Locksmith
