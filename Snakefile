configfile: "config.json"
rule all:
	input: 
		"data/high_targets_conv.vcf",
		"data/target_sequences.fa",
		"data/CpG_in_targets.bed"


#---------------------------------------------------------------------------------------------------
#STEP 1: GET TARGET SEQUENCES
rule create_bed_file_range:
	input:
		"data/target_list.bed"
	output:
		"data/target_list_range.bed"
	shell:
		"python scripts/target_2_target_range_list.py -f {input} -o {output} -r {config[target_range]}"

rule extract_target_sequences:
	input:
		genome="data/hg19.fa",
		bed="data/target_list_range.bed"
	output:
		"data/target_sequences.fa"
	shell:
		"bedtools getfasta -fi {input.genome} -bed {input.bed}  -fo {output}"

#---------------------------------------------------------------------------------------------------
#STEP 2: FIND CpG AND SNP IN TARGET SEQUENCES

rule convert_bed_file:
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
        input:
                "data/target_list_range.bed"
        output:
                "data/target_list_range_conv.bed"
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

rule find_overlapping_VCF:
	input:
		"data/target_list_range_conv.bed"
	output:
		"data/targets.vcf"	
	shell:
		"tabix -f ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz -R {input} > {output}"

rule keep_high_freq_SNP:
	input:
		"data/targets.vcf"
	output:
		"data/high_targets.vcf"
	shell:
		"python scripts/extract_high_frequent_SNPs.py -i {input} -o {output} -f {config[SNP_frequency_threshold]}"

rule find_CpGs_in_targets:
	input:
		"data/target_list_range.bed"
	output:
		"data/CpG_in_targets.bed"
	shell:
		"intersectBed -a {input} -b ~/genomes/hg19/CpG_bed/no_header/CpG_hg19_chr* > {output}"

rule convert_vcf_file:
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
        input:
                "data/high_targets.vcf"
        output:
                "data/high_targets_conv.vcf"
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
                "sed -e 's/NC_000023.10\t/chr1\X/' |sed -e 's/NC_000024.9\t/chrY\t/' "
                "> {output}"
