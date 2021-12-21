configfile:	config['path_to_config_file']
rule all:
	input:
		fasta=config['output_directory']+"/chosen_panel_iterative.fasta"
#---------------------------------------------------------------------------------------------------
rule copy_config_to_output_directory:
	input:
		 config['path_to_config_file']
	output:
		config['output_directory']+"/config.json"
	shell:
		"cp {input} {output}"
#STEP 1: GET TARGET SEQUENCES
rule create_bed_file_range:
	input:
		target=config['target_list_bedfile_path'],
		c=config['output_directory']+"/config.json"
	output:
		config['output_directory']+"/target_list_range.bed"
	shell:
		"python {config[path_to_scripts]}target_2_target_range_list.py -f {input.target} -o {output} -r {config[target_range]}"

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
		config_file=config['output_directory']+"/config.json",
		bed=config['output_directory']+"/target_list_range.bed"
	output:
		tms=config['output_directory']+"/tms.csv",
                cpg=config['output_directory']+"/CpG_conflicts.csv",
                probes=config['output_directory']+"/possible_probes_all_targets.csv",
                hairpins=config['output_directory']+"/hairpin_scores.csv",
                snp=config['output_directory']+"/SNP_conflicts.csv",
		arms=config['output_directory']+"/probe_arms.csv"

	log:
		out = config['output_directory']+"/logs/create_probes_stdout.log",
		err = config['output_directory']+"/logs/create_probes_stderr.err"

	shell:
		"python {config[path_to_scripts]}create_all_arm_combinations.py -i {input.targets} -o {config[output_directory]} -c {input.config_file} -b {input.bed} -t {config[max_threads]} 2> {log.err} 1> {log.out}"
#--------------------------------------------------------------------------------------------------
#STEP 3
rule choose_panel_iteratively:
        input:
                config_file=config['output_directory']+"/config.json",
                tms=config['output_directory']+"/tms.csv",
                cpg=config['output_directory']+"/CpG_conflicts.csv",
                probes=config['output_directory']+"/possible_probes_all_targets.csv",
                hairpins=config['output_directory']+"/hairpin_scores.csv",
                snp=config['output_directory']+"/SNP_conflicts.csv",
		arms=config['output_directory']+"/probe_arms.csv",
                target=config['output_directory']+"/target_list_range.bed"

	output:
                fasta=config['output_directory']+"/chosen_panel_iterative.fasta",
		conflicts=config['output_directory']+"/conflicting_probes_per_iteration.csv"

	log:
		out = config['output_directory']+"/logs/choose_panel_iteratively_stdout.log",
		err = config['output_directory']+"/logs/choose_panel_iteratively_stderr.err"

	shell:
                "python {config[path_to_scripts]}choose_panel_iterative.py -o {output.fasta} -c {input.config_file} -m {input.tms} -g {input.cpg} -p {input.probes} -a {input.hairpins} -n {input.snp} -r {input.target} -b {input.arms} -t {config[max_threads]} -u {config[output_directory]} -f {output.conflicts} 2> {log.err} 1> {log.out}"

