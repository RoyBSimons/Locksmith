rule extract_target_sequences:
	input:
	"data/genome.fa"
	"data/target_list"
	output:
	"data/target_sequences.fa"
