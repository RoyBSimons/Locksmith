configfile: "config.json"
	rule all:
	input:
		expand(config['output_directory']+"/selection_{selection_round}/conflicting_probes_hairpins.tsv",selection_round=range(0,int(config["Selection_rounds_for_QC"])))
	shell:
		"cp config.json "+config['output_directory']+"/config.json"
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
#STEP 2: FIND CpG AND SNP IN TARGET SEQUENCES

rule convert_bed_file: #Convert the chromosome names from chr1 to NC_000001.10 format to be compatible with the ftp assembly in the find_overlapping_VCF rule
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
        input:
                config['output_directory']+"/target_list_range.bed"
        output:
                config['output_directory']+"/target_list_range_conv.bed"
        shell:
                "sed -e 's/chr1\t/NC_000001.10\t/' {input} | sed -e 's/chr2\t/NC_000002.11\t/' | "
                "sed -e 's/chr3\t/NC_000003.11\t/' |sed -e 's/chr4\t/NC_000004.11\t/' | "
                "sed -e 's/chr5\t/NC_000005.9\t/' |sed -e 's/chr6\t/NC_000006.11\t/' | "
                "sed -e 's/chr7\t/NC_000007.13\t/' |sed -e 's/chr8\t/NC_000008.10\t/' | "
                "sed -e 's/chr9\t/NC_000009.11\t/' |sed -e 's/chr10\t/NC_000010.10\t/' | "
                "sed -e 's/chr11\t/NC_000011.9\t/' |sed -e 's/chr12\t/NC_000012.11\t/' | "
                "sed -e 's/chr13\t/NC_000013.10\t/' |sed -e 's/chr14\t/NC_000014.8\t/' | "
                "sed -e 's/chr15\t/NC_000015.9\t/' |sed -e 's/chr16\t/NC_000016.9\t/' | "
                "sed -e 's/chr17\t/NC_000017.10\t/' |sed -e 's/chr18\t/NC_000018.9\t/' | "
                "sed -e 's/chr19\t/NC_000019.9\t/' |sed -e 's/chr20\t/NC_000020.10\t/' | "
                "sed -e 's/chr21\t/NC_000021.8\t/' |sed -e 's/chr22\t/NC_000022.10\t/' | "
                "sed -e 's/chrX\t/NC_000023.10\t/' |sed -e 's/chrY\t/NC_000024.9\t/' "
                "> {output}"

rule find_overlapping_VCF:#Obtain all known SNPs 
	input:
		config['output_directory']+"/target_list_range_conv.bed"
	output:
		config['output_directory']+"/targets.vcf"	
	shell:
		"tabix ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz -R {input} > {output}"

rule keep_high_freq_SNP: #Keep the SNPs that have a frequency that is equal or higher than the threshold in the config file:SNP_frequency_threshold.
	input:
		config['output_directory']+"/targets.vcf"
	output:
		config['output_directory']+"/high_targets.vcf"
	shell:
		"python scripts/extract_high_frequent_SNPs.py -i {input} -o {output} -f {config[SNP_frequency_threshold]}"

rule find_CpGs_in_targets:
	input:
		config['output_directory']+"/target_list_range.bed"
	output:
		config['output_directory']+"/CpG_in_targets.bed"
	shell:
		"intersectBed -a {input} -b data/CpG_bed_nomenclature/CpG_hg19_chr* -wa -wb> {output}"

rule convert_vcf_file:
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
        input:
                config['output_directory']+"/high_targets.vcf"
        output:
                tmp_snp=config['output_directory']+"/SNP_in_targets.vcf",
		snp=config['output_directory']+"/SNP_in_targets.bed"
        shell:
                "sed -e 's/NC_000001.10\t/chr1\t/' {input} | sed -e 's/NC_000002.11\t/chr2\t/' | "
                "sed -e 's/NC_000003.11\t/chr3\t/' |sed -e 's/NC_000004.11\t/chr4\t/' | "
                "sed -e 's/NC_000005.9\t/chr5\t/' |sed -e 's/NC_000006.11\t/chr6\t/' | "
                "sed -e 's/NC_000007.13\t/chr7\t/' |sed -e 's/NC_000008.10\t/chr8\t/' | "
                "sed -e 's/NC_000009.11\t/chr9\t/' |sed -e 's/NC_000010.10\t/chr10\t/' | "
                "sed -e 's/NC_000011.9\t/chr11\t/' |sed -e 's/NC_000012.11\t/chr12\t/' | "
                "sed -e 's/NC_000013.10\t/chr13\t/' |sed -e 's/NC_000014.8\t/chr14\t/' | "
                "sed -e 's/NC_000015.9\t/chr15\t/' |sed -e 's/NC_000016.9\t/chr16\t/' | "
                "sed -e 's/NC_000017.10\t/chr17\t/' |sed -e 's/NC_000018.9\t/chr18\t/' | "
                "sed -e 's/NC_000019.9\t/chr19\t/' |sed -e 's/NC_000020.10\t/chr20\t/' | "
                "sed -e 's/NC_000021.8\t/chr21\t/' |sed -e 's/NC_000022.10\t/chr22\t/' | "
                "sed -e 's/NC_000023.10\t/chrX\t/' |sed -e 's/NC_000024.9\t/chrY\t/' "
                "> {output.tmp_snp} && "
		"""cat {output.tmp_snp} | awk -F"\\t" "BEGIN {{OFS = FS}}{{print \$1,\$2-1,\$2,\$0}}" | cut -f -3,6- > {output.snp}"""

#---------------------------------------------------------------------------------------------------
#STEP 3
rule create_all_arm_combinations:
	input:
		targets=config['output_directory']+"/target_sequences.fa",
		config_file="config.json",
		bed=config['output_directory']+"/target_list_range.bed"
	output:
		tsv=config['output_directory']+"/all_arms.tsv", 
		up_tmp=temp(config['output_directory']+"/all_arms_upstream_tmp.bed"),
                up=config['output_directory']+"/all_arms_upstream.bed",
                down_tmp=temp(config['output_directory']+"/all_arms_downstream_tmp.bed"),
		down=config['output_directory']+"/all_arms_downstream.bed"
	shell:
		"python scripts/create_all_arm_combinations.py -i {input.targets} -o {output.tsv} -c {input.config_file} -u {output.up_tmp} -d {output.down_tmp} -b {input.bed} && "
		"cat {output.up_tmp} | tr -d '\r' > {output.up} && "
                "cat {output.down_tmp} | tr -d '\r' > {output.down}"
#---------------------------------------------------------------------------------------------------
#STEP 4
rule report_SNP_and_CpG_in_arms:
	input:
		snp=config['output_directory']+"/SNP_in_targets.bed",
		down=config['output_directory']+"/all_arms_downstream.bed",
		up=config['output_directory']+"/all_arms_upstream.bed",
		arms=config['output_directory']+"/all_arms.tsv"
	output:
		cpg_up=config['output_directory']+"/all_arms_upstream_cpg.bed",
		cpg_down=config['output_directory']+"/all_arms_downstream_cpg.bed",
		snp_up=config['output_directory']+"/all_arms_upstream_snp.bed",
		snp_down=config['output_directory']+"/all_arms_downstream_snp.bed",
		arms=config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms.tsv"
	params:
		value='{iteration}',
		iteration='{iteration}'

	shell:
		"bedtools intersect -a {input.up} -b {config[cpg_bed_path]} -c > {output.cpg_up} && "
		"bedtools intersect -a {input.down} -b {config[cpg_bed_path]} -c > {output.cpg_down} && "
                "bedtools intersect -a {input.up} -b {input.snp} -c > {output.snp_up} && "
                "bedtools intersect -a {input.down} -b {input.snp} -c > {output.snp_down} && "
		"python scripts/remove_CpG_SNP_conflicts.py -i {input.arms} -u {output.cpg_up} -d {output.cpg_down} -s {output.snp_up} -t {output.snp_down} -o {output.arms} -l {params.value}"
		
def get_correct_wildcard(wildcards):
        if int(wildcards.iteration)== 0:
                output_value=config['output_directory']+'/Arm_selection_initialization_file'
        else:
                output_value=config['output_directory']+"/iteration_"+str(int(wildcards.iteration)-1)+"/not_selected_arms_"+str(int(wildcards.iteration)-1)+"_0.tsv"
        return output_value

rule Initialize_arm_selection:
        input:
        output:
                iteration_initiation=config['output_directory']+"/Arm_selection_initialization_file"
        shell:
                "echo This is an empty file, which fills the place of input.non in the first iteration of arm selection > {output.iteration_initiation}"
rule Initialize_probe_selection:
        input:
        output:
                selection_initiation=config['output_directory']+"/probe_selection_initialization_file"
        shell:
                "touch {output.selection_initiation}"

rule Select_arms_iterative: #This rule only keeps the arms from "output/iteration_{iteration}/cpg_snp_{iteration}_arms.tsv" which are not yet selected by previous iterations
        input:
                non=lambda wildcards: get_correct_wildcard(wildcards),
                arms=config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms.tsv"
        output:
                config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms_chosen.tsv"
        shell:
                "python scripts/choose_arms_for_iteration.py -i {input.arms} -n {input.non} -o {output}"

#---------------------------------------------------------------------------------------------------
#STEP 5
rule Obtain_Tm_arms:
	input:
		config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms_chosen.tsv"
	output:
		up=config['output_directory']+"/iteration_{iteration}/arms_tm_upstream_{iteration}.txt",
		down=config['output_directory']+"/iteration_{iteration}/arms_tm_downstream_{iteration}.txt"
	shell:
		"""seq_list_down="" &&"""
		""" seq_list_up="" &&"""
		""" while IFS= read line; do seq_list_down=$seq_list_down$(awk "{{print \$1}}"); done < "{input}" &&"""
		""" while IFS= read line; do seq_list_up=$seq_list_up$(awk "{{print \$2}}"); done < "{input}" &&"""
		""" lines_in_file=$(wc -l <{input}) &&"""
		"""if [ $lines_in_file -gt 1 ]; then echo $seq_list_down | xargs parallel -k -j 400 oligotm ::: >> {output.down} &&  echo $seq_list_up | xargs parallel -k -j 400 oligotm ::: >> {output.up}; else echo No Tms are obtained as there are no arms in the input file && touch {output.down} && touch {output.up}; fi; """

rule Add_Tm_arms:
	input:
		up=config['output_directory']+"/iteration_{iteration}/arms_tm_upstream_{iteration}.txt",
		down=config['output_directory']+"/iteration_{iteration}/arms_tm_downstream_{iteration}.txt",
		arms=config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms_chosen.tsv"
	output:
		config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms_tm.tsv"
	shell:
		"python scripts/add_tm_to_arms.py -i {input.arms} -o {output} -d {input.down} -u {input.up}"
		
#---------------------------------------------------------------------------------------------------
#STEP 6
def get_correct_wildcard_selection(name):
        selection_nr=int(name[2:-2])
        if selection_nr== 0:
                output_value=config['output_directory']+"/probe_selection_initialization_file"
        else:
                output_value=config['output_directory']+"/selection_"+str(selection_nr-1)+"/conflicting_probes_hairpins_dimers.txt"
        return output_value

def get_correct(wildcards):
	if int(wildcards.iteration)== 0:
		output_value=config['output_directory']+"/Arm_selection_initialization_file"
	else:
		output_value=config['output_directory']+"/iteration_"+str(int(wildcards.iteration)-1)+"/selected_arms_"+str(int(wildcards.iteration)-1)+"_"+str(int(wildcards.selection_round))+".tsv"
	return output_value
def get_correct_conflict_hairpins(wildcards):
	if int(wildcards.selection_round)== 0:
		output_value=config['output_directory']+"/Arm_selection_initialization_file"
	else:
		output_value=config['output_directory']+"/selection_"+str(int(wildcards.selection_round)-1)+"/conflicting_probes_hairpins_combined.tsv"
	return output_value
def get_correct_conflict_dimers(wildcards):
        if int(wildcards.selection_round)== 0:
                output_value=config['output_directory']+"/Arm_selection_initialization_file"
        else:
                output_value=config['output_directory']+"/selection_"+str(int(wildcards.selection_round)-1)+"/conflicting_probes_dimers_combined.tsv"
        return output_value

rule Select_arms:
	input:
		arms=config['output_directory']+"/iteration_{iteration}/cpg_snp_{iteration}_arms_tm.tsv",
		config_file="config.json",
		target="data/target_list.bed",
		conflicts_hairpins=lambda wildcards: get_correct_conflict_hairpins(wildcards),
		conflicts_dimers=lambda wildcards: get_correct_conflict_dimers(wildcards),
		previous_selected=lambda wildcards: get_correct(wildcards),
		up=config['output_directory']+"/iteration_0/arms_tm_upstream_0.txt",
		down=config['output_directory']+"/iteration_0/arms_tm_downstream_0.txt"

	output:
		selected=config['output_directory']+"/iteration_{iteration}/selected_arms_{iteration}_{selection_round}.tsv",
		not_selected=config['output_directory']+"/iteration_{iteration}/not_selected_arms_{iteration}_{selection_round}.tsv"
	params:
		iteration='{iteration}',
		selection_round='{selection_round}'
	shell:
		"python scripts/select_arms.py -i {input.arms} -o {output.selected} -n {output.not_selected} -t {input.target} -y {params.iteration} -d {input.conflicts_dimers} -c {input.conflicts_hairpins} -s {params.selection_round} -a {input.up} -b {input.down}"

rule Combine_selected_arms:
	input: 
		expand(config['output_directory']+"/iteration_{iteration}/selected_arms_{iteration}_{{selection_round}}.tsv",iteration=range(0,1+int(config['probe_specifics'][0]['max_cpgs_in_arms'])))
	output:
		config['output_directory']+"/selection_{selection_round}/selected_arms_combined.tsv"
	shell:
		"awk 'FNR==1 && NR!=1{{next;}}{{print}}' "+config['output_directory']+"/iteration_*/selected*_{wildcards.selection_round}.tsv > {output}"
#---------------------------------------------------------------------------------------------------
#STEP 7
rule Create_probes:
	input:
                config['output_directory']+"/selection_{selection_round}/selected_arms_combined.tsv"
	output:
		config['output_directory']+"/selection_{selection_round}/probes.fasta"
	shell:
		"python scripts/Add_backbone.py -i {input} -c {config[Backbone_sequence]} -o {output}"
		
#---------------------------------------------------------------------------------------------------
#STEP 8
rule probe_QC_by_MFEprimer:
	input:
		config['output_directory']+"/selection_{selection_round}/probes.fasta"
	output:
		hairpin=config['output_directory']+"/selection_{selection_round}/probe_QC_hairpin.json",
                dimer=config['output_directory']+"/selection_{selection_round}/probe_QC_dimer.json"
	shell:
#		"db_caller='' &&"
#		"databases={config[PATH_to_genome_database_directory]}*.fa &&"
#		"for file in $databases; do db_caller=$db_caller'-d '$file' ' ; done &&"
#		"mfeprimer -i {input} -o {output} --json $db_caller"
		"outputname_hairpin={output.hairpin} && outputname_without_suffix_hairpin=${{outputname_hairpin:0: -5}} &&"
                "mfeprimer hairpin -i {input} -o $outputname_without_suffix_hairpin --json &&"
                "outputname_dimer={output.dimer} && outputname_without_suffix_dimer=${{outputname_dimer:0: -5}} &&"
                "mfeprimer dimer -i {input} -o $outputname_without_suffix_dimer --json"


rule Rerun_probes_with_hairpins_or_dimers:
	input:
		dimer=config['output_directory']+"/selection_{selection_round}/probe_QC_dimer.json",
		hairpin=config['output_directory']+"/selection_{selection_round}/probe_QC_hairpin.json",
		config="config.json"
	output:
		conflicts_hairpin=config['output_directory']+"/selection_{selection_round}/conflicting_probes_hairpins.tsv",
		conflicts_dimer=config['output_directory']+"/selection_{selection_round}/conflicting_probes_dimers.tsv",
		combined_hairpins=config['output_directory']+"/selection_{selection_round}/conflicting_probes_hairpins_combined.tsv",
                combined_dimers=config['output_directory']+"/selection_{selection_round}/conflicting_probes_dimers_combined.tsv"
	shell:
		"python scripts/rerun_targets_from_hairpin_and_dimers.py -d {input.dimer} -p {input.hairpin} -o {output.conflicts_hairpin} -c {output.conflicts_dimer} -i {input.config} &&"
		"cat "config['output_directory']+"/selection_*/conflicting_probes_hairpins.tsv > {output.combined_hairpins} &&"
		"cat "config['output_directory']+"/selection_*/conflicting_probes_dimers.tsv > {output.combined_dimers}"
