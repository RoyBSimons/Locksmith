configfile: "config.json"
rule all:
	input:
		fasta=config['output_directory']+"/chosen_panel.fasta"
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
		genome=config['PATH_to_reference_genome_fasta'],
		bed=config['output_directory']+"/target_list_range.bed"
	output:
		config['output_directory']+"/target_sequences.fa"
	shell:
		"bedtools getfasta -fi {input.genome} -bed {input.bed}  -fo {output}"

#---------------------------------------------------------------------------------------------------
#STEP 2
rule create_probes:
	input:
		targets=config['output_directory']+"/target_sequences.fa",
		config_file="config.json",
		bed=config['output_directory']+"/target_list_range.bed"
	output:
		tms=config['output_directory']+"/tms.csv",
                cpg=config['output_directory']+"/CpG_conflicts.csv",
                probes=config['output_directory']+"/possible_probes_all_targets.csv",
                hairpins=config['output_directory']+"/hairpin_scores.csv",
                snp=config['output_directory']+"/SNP_conflicts.csv",
                scores=config['output_directory']+"/probe_scores.csv"
	shell:
		"python scripts/create_all_arm_combinations.py -i {input.targets} -o {config[output_directory]} -c {input.config_file} -b {input.bed} -t {config[max_threads]}"
#--------------------------------------------------------------------------------------------------
#STEP 3
rule create_and_choose_panel:
	input:
		config_file="config.json",
		tms=config['output_directory']+"/tms.csv",
		cpg=config['output_directory']+"/CpG_conflicts.csv",
		probes=config['output_directory']+"/possible_probes_all_targets.csv",
		hairpins=config['output_directory']+"/hairpin_scores.csv",
		snp=config['output_directory']+"/SNP_conflicts.csv",
		scores=config['output_directory']+"/probe_scores.csv"
	output:
		fasta=config['output_directory']+"/chosen_panel.fasta"
	shell:
		"python scripts/create_and_choose_panel.py -o {output.fasta} -c {input.config_file} -m {input.tms} -g {input.cpg} -p {input.probes} -a {input.hairpins} -n {input.snp} -s {input.scores} -t {config[max_threads]}"
