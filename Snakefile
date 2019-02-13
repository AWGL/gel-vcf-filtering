###########################################################################################################################
# Genomics England VCF filtering
###########################################################################################################################

config_location = "config/config_vm_validation.yaml"

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


def get_family_index(wildcards):

	file_list = []

	family_members = config['families'][wildcards.family_id]

	for member in family_members:

		file_list.append("output/roi/{wildcards.family_id}/{member}_roi.vcf.gz.tbi".format(wildcards=wildcards, member=member))

	return file_list

###########################################################################################################################
# End Target Rule
###########################################################################################################################

families = config["families"]

rule all:
	input:
		expand("output/filtered_csv/{family_id}/{family_id}_finished.txt", family_id = families),
		expand("output/peddy/{family_id}/{family_id}.html", family_id = families)

###########################################################################################################################
# Main program
###########################################################################################################################


# remove '-' from sample names and replace with '_'
rule fix_sample_names:
	input:
		vcf = "input/{family_id}/{sample}.vcf.gz",
		index = "input/{family_id}/{sample}.vcf.gz.tbi"
	output:
		vcf = temp("output/fixed_name_vcfs/{family_id}/{sample}_fixed_names.vcf.gz"),
		index = temp("output/fixed_name_vcfs/{family_id}/{sample}_fixed_names.vcf.gz.tbi"),
		original_name = "output/temp/{family_id}/{sample}_original_name.txt",
		new_name = "output/temp/{family_id}/{sample}_new_name.txt"
	shell:
		"""
		bcftools query --list-samples {input.vcf} > {output.original_name}

		cat {output.original_name} | tr "-" "_" | tr "." "_" > {output.new_name}
		
		bcftools reheader --samples {output.new_name} {input.vcf} | bcftools view -O z > {output.vcf}

		tabix {output.vcf}

		"""

# Take the region of interest bed file and filter each WGS bed file
rule filter_by_roi:
	input:
		vcf = "output/fixed_name_vcfs/{family_id}/{sample}_fixed_names.vcf.gz",
		index = "output/fixed_name_vcfs/{family_id}/{sample}_fixed_names.vcf.gz.tbi"
	output:
		vcf = "output/roi/{family_id}/{sample}_roi.vcf"
	params:
		roi_bed = config["roi_bed"]
	shell:
		"bcftools view -Ov -T {params.roi_bed} "
		"{input.vcf} > {output.vcf} "


# Compress and index the vcf
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
		vcf = lambda wildcards: get_family(wildcards),
		index = lambda wildcards: get_family_index(wildcards)
	output:
		"output/merged/{family_id}/{family_id}_merged.vcf.gz"
	params:
		vcfs = lambda wildcards, input: " ".join(input.vcf)
	shell:
		"bcftools merge {params.vcfs} | bgzip -c > {output} && tabix {output}"

# Create a PED file for downstream analysis
rule create_ped_file:
	input:
		config_location
	output:
		"output/all_ped/all_families.ped"
	shell:
		"python scripts/make_ped.py --config {input} > {output}"

# Create a PED file for each family
rule create_per_family_ped:
	input:
		ped = "output/all_ped/all_families.ped",
		family = "output/merged/{family_id}/{family_id}_merged.vcf.gz"
	output:
		"output/family_ped/{family_id}/{family_id}.ped"
	shell:
		"grep {wildcards.family_id} {input.ped} > {output}"

# Run peddy to check for problems with ped file. Spot errors when setting up config
rule run_peddy:
	input:
		vcf = "output/merged/{family_id}/{family_id}_merged.vcf.gz",
		ped = "output/family_ped/{family_id}/{family_id}.ped"
	output:
		"output/peddy/{family_id}/{family_id}.background_pca.json",
		"output/peddy/{family_id}/{family_id}.het_check.csv",
		"output/peddy/{family_id}/{family_id}.het_check.png",
		"output/peddy/{family_id}/{family_id}.html",
		"output/peddy/{family_id}/{family_id}.pca_check.png",
		"output/peddy/{family_id}/{family_id}.ped_check.csv",
		"output/peddy/{family_id}/{family_id}.ped_check.png",
		"output/peddy/{family_id}/{family_id}.ped_check.rel-difference.csv",
		"output/peddy/{family_id}/{family_id}.peddy.ped",
		"output/peddy/{family_id}/{family_id}.sex_check.csv",
		"output/peddy/{family_id}/{family_id}.sex_check.png",
		"output/peddy/{family_id}/{family_id}.vs.html"
	threads:
		config["peddy_threads"]
	conda:
		"envs/python2.yaml"
	shell:
		"peddy --plot -p {threads} {input.vcf} {input.ped} --sites hg38 --prefix output/peddy/{wildcards.family_id}/{wildcards.family_id} "

# annotate with vep:
rule annotate_with_vep:
	input:
		"output/merged/{family_id}/{family_id}_merged.vcf.gz"
	output:
		"output/merged_vep/{family_id}/{family_id}_merged_vep.vcf"
	params:
		vep_cache = config["vep_cache_location"],
		ref = config["reference"],
		gnomad_genomes = config["gnomad_genomes"],
		gnomad_exomes = config["gnomad_exomes"],
		ccrs = config["ccrs"],
		spliceai = config["spliceai"],
		cadd_snvs = config["cadd_snvs"],
		cadd_indels = config["cadd_snvs"]
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
		"--custom {params.ccrs},ccrs,bed,overlap,0 "
		"--custom {params.spliceai},SpliceAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL,SYMBOL "
		"--plugin CADD,{params.cadd_snvs},{params.cadd_indels}"

# Convert to CSV using gatk
rule vcf_to_table:
	input:
		"output/merged_vep/{family_id}/{family_id}_merged_vep.vcf"
	output:
		temp("output/csv/{family_id}/{family_id}_merged_vep.csv")
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

# remove 'chr' string from chromosomes e.g. chr1 > 1
rule fix_chromsome_notation:
	input:
		"output/csv/{family_id}/{family_id}_merged_vep.csv"
	output:
		temp("output/csv_fixed/{family_id}/{family_id}_merged_vep_fixed.csv")
	shell:
		"awk '{{gsub(/^chr/,\"\"); print}}' {input} > {output}   "


# Rule get CSQ string
rule get_csq_string:
	input:
		"output/merged_vep/{family_id}/{family_id}_merged_vep.vcf"
	output:
		"output/csq/{family_id}/{family_id}_csq.txt"
	shell:
		"grep \"^##INFO=<ID=CSQ\" {input} | awk 'BEGIN {{ FS = \":\" }} ; {{ print $2 }}' | tr -d '>\" ' > {output}     "


# Run germline_filtering_script
rule filter_variants:
	input:
		csv = "output/csv_fixed/{family_id}/{family_id}_merged_vep_fixed.csv",
		ped = "output/family_ped/{family_id}/{family_id}.ped",
		csq = "output/csq/{family_id}/{family_id}_csq.txt"
	output:
		"output/filtered_csv/{family_id}/{family_id}_finished.txt"
	params:
		germline_variant_filter = config["germline_variant_filter"],
		filter_config = config["filter_config"],
		local_panel_app_dump = config["local_panel_app_dump"],
		hpo_file = config["hpo_file"]
	shell:
		"python {params.germline_variant_filter} "
		"--config {params.filter_config} "
		"--ped {input.ped} "
		"--input {input.csv} "
		"--panelapp "
		"--local-panel-app-dump {params.local_panel_app_dump} "
		"--spliceai "
		"--smart-synonymous "
		"--gnomad-constraint-scores "
		"--worksheet {wildcards.family_id} "
		"--results-dir output/filtered_csv/{wildcards.family_id}/ "
		"--csq $(cat {input.csq}) "
		"--patient-hpos {params.hpo_file} "
		"--add-ccrs && touch {output}"





























