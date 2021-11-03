configfile: "config.json"
rule all:
	input:
		fasta=config['output_directory']+"/chosen_probes.fasta",
                tms=config['output_directory']+"/tms.csv",
                cpg=config['output_directory']+"/CpG_conflicts.csv",
                probes=config['output_directory']+"/possible_probes_all_targets.csv",
                hairpins=config['output_directory']+"/hairpin_scores.csv",
                snp=config['output_directory']+"/SNP_conflicts.csv",
                scores=config['output_directory']+"/probe_scores.csv"

#---------------------------------------------------------------------------------------------------
#STEP 1: GET TARGET SEQUENCES
rule create_bed_file_range:
	input:
		"data/target_list.bed"
	output:
		config['output_directory']+"/target_list_range.bed"
	shell:
		"python scripts/target_2_target_range_list.py -f {input} -o {output} -r {config[target_range]}"

rule extract_target_sequences:
	input:
		genome="data/hg19.fa", #change this to ftp?
		bed=config['output_directory']+"/target_list_range.bed"
	output:
		config['output_directory']+"/target_sequences.fa"
	shell:
		"bedtools getfasta -fi {input.genome} -bed {input.bed}  -fo {output}"

#---------------------------------------------------------------------------------------------------
rule keep_high_freq_SNP: #Keep the SNPs that have a frequency that is equal or higher than the threshold in the config file:SNP_frequency_threshold.
	input:
		config['output_directory']+"/targets.vcf"
	output:
		config['output_directory']+"/high_targets.vcf"
	shell:
		"python scripts/extract_high_frequent_SNPs.py -i {input} -o {output} -f {config[SNP_frequency_threshold]}"
#---------------------------------------------------------------------------------------------------
#STEP 3
rule create_all_arm_combinations:
	input:
		targets=config['output_directory']+"/target_sequences.fa",
		config_file="config.json",
		bed=config['output_directory']+"/target_list_range.bed"
	output:
		fasta=config['output_directory']+"/chosen_probes.fasta",
		tms=config['output_directory']+"/tms.csv",
                cpg=config['output_directory']+"/CpG_conflicts.csv",
                probes=config['output_directory']+"/possible_probes_all_targets.csv",
                hairpins=config['output_directory']+"/hairpin_scores.csv",
                snp=config['output_directory']+"/SNP_conflicts.csv",
                scores=config['output_directory']+"/probe_scores.csv"

	threads: 20
	shell:
		"python scripts/create_all_arm_combinations.py -i {input.targets} -o {config[output_directory]} -c {input.config_file} -b {input.bed} -t {threads}"

