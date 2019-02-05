###########################################################################################################################
# Genomics England VCF filtering
###########################################################################################################################

config_location = "config/config_vm.yaml"

configfile: config_location

###########################################################################################################################
# Utility Functions for getting files
###########################################################################################################################

def get_family(wildcards):

	file_list = []

	family_members = config['families'][wildcards.family_id]

	for member in family_members:

		file_list.append("output/roi/{wildcards.family_id}/{member}_roi.vcf.gz".format(wildcards=wildcards, member=member))

	return file_list
	
###########################################################################################################################
# Main program
###########################################################################################################################

# Take the region of interest bed file and filter each WGS bed file
rule filter_by_roi:
	input:
		vcf = "input/{family_id}/{sample}.vcf.gz",
		index = "input/{family_id}/{sample}.vcf.gz.tbi"
	output:
		vcf = "output/roi/{family_id}/{sample}_roi.vcf"
	params:
		roi_bed = config["roi_bed"]
	shell:
		"bcftools view -Ov -T {params.roi_bed} "
		"{input.vcf} > {output.vcf} "


# Take the region of interest bed file and filter each WGS bed file
rule compress_and_index_vcf:
	input:
		vcf = "output/roi/{family_id}/{sample}_roi.vcf"
	output:
		vcf = "output/roi/{family_id}/{sample}_roi.vcf.gz",
		index = "output/roi/{family_id}/{sample}_roi.vcf.gz.tbi"
	params:
		roi_bed = config["roi_bed"]
	shell:
		"bgzip {input.vcf} && tabix {output.vcf}"

# Merge vcf by family
rule merge_vcfs:
	input:
		lambda wildcards: get_family(wildcards)
	output:
		"output/merged/{family_id}_merged.vcf.gz"
	params:
		vcfs = lambda wildcards, input: " ".join(input)
	shell:
		"bcftools merge {params.vcfs} | bgzip -c > {output} && tabix {output}"


# Create a PED file for downstream analysis
rule create_ped_file:
	input:
		config_location
	output:
		"output/config/all_families.ped"
	shell:
		"python scripts/make_ped.py --config {input} > {output}"


# Check the ped file makes sense using peddy
rule ped_check:
	input:
		"output/merged/{family_id}_merged.vcf.gz"
	output:
		pass





# annotate with vep:
rule annotate_with_vep:
	input:
		"output/merged/{family_id}.vcf.gz"
	output:
		"output/merged_vep/{family_id}_merged_vep.vcf"
	params:
		vep_cache = config["vep_cache_location"],
		ref = config["reference"],
		gnomad_genomes = config["gnomad_genomes"],
		gnomad_exomes = config["gnomad_exomes"]
	threads:
		config["vep_threads"]
	shell:
		"vep --verbose "
		"--format vcf "
		"--everything "
		"--fork {threads} "
		"--species homo_sapiens "
		"--assembly GRCh38  "
		"--input_file {input}  "
		"--output_file {output} "
		"--force_overwrite "
		"--cache "
		"--dir  {params.vep_cache} "
		"--fasta {params.ref} "
		"--offline "
		"--cache_version 94 "
		"--no_escape "
		"--shift_hgvs 1 "
		"--vcf "
		"--refseq "
		"--exclude_predicted "
		"--flag_pick "
		"--pick_order biotype,canonical,appris,tsl,ccds,rank,length "
		"--custom {params.gnomad_genomes},gnomADg,vcf,exact,0,AF_POPMAX "
		"--custom {params.gnomad_exomes},gnomADe,vcf,exact,0,AF_POPMAX "
		"--custom {params.ccrs},ccrs,bed,overlap,1 "
		"--custom {params.spliceai},SpliceAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL,SYMBOL "

# Convert to CSV using gatk
rule vcf_to_table:
	input:
		"output/merged_vep/{family_id}_merged_vep.vcf"
	output:
		"output/csv/{family_id}_merged_vep.csv"
	shell:
		"gatk VariantsToTable "
		"-V {input} -O {output} "
		"-F CHROM "
		"-F POS "
		"-F REF "
		"-F ALT "
		"-F ID "
		"-F QUAL "
		"-F FILTER "
		"-F CSQ "
		"-F TC "
		"-F AC "
		"-F TYPE "
		"-GF GT "
		"-GF GQ "
		"-GF NR "
		"-GF NV"


# Run germline_filtering_script


