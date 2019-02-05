import yaml
import argparse
"""
Script for parsing the YAML config file for the pipeline and creating a PED file.

Usage:

python make_ped.py config.yaml > myped.ped

Requires yaml - will be installed as part of smakemake conda environment

"""

# Get parser arguments
parser = argparse.ArgumentParser(description='Make a PED file from the variable files')
parser.add_argument('--config', type=str, nargs=1, required=True,
					help='The location of the YAML config')
args = parser.parse_args()
config = args.config[0]

# Parse YAML
config_dict = yaml.load(open(config))


# loop through samples and create ped file.
for sample in config_dict['sample_info'].keys():

	family_id = config_dict['sample_info'][sample]['familyId']
	sample_id = sample
	paternal_id = config_dict['sample_info'][sample]['paternalId']
	maternal_id = config_dict['sample_info'][sample]['maternalId']
	sex = config_dict['sample_info'][sample]['sex']
	phenotype = config_dict['sample_info'][sample]['phenotype']


	# Set values if null in config
	if family_id == None:
		family_id = 0

	if paternal_id == None:
		paternal_id = 0

	if maternal_id == None:
		maternal_id = 0

	if sex == None:
		sex = 0

	if phenotype == None:
		phenotype = 2

	print (f'{family_id}\t{sample_id}\t{paternal_id}\t{maternal_id}\t{sex}\t{phenotype}')