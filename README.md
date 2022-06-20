# Cobra functions
This is a repository of COnstraint-Based Reconstruction (COnstraint-Based Reconstruction and Analysis) scripts and functions used to build organ-specific personalized flux maps from genetically imputed transcript abundances. First, a reference flux map is computed for each organ using the GIM3E algorithm to integrate average transcript abundances and organ-specific metabolic functions into organ-specific Genome-Scale Metabolic Models (GSMMs). In parallel, transcript abundances imputed from genotype data are mapped to reactions in the metabolic network. Then, the quadratic metabolic transformation algorithm (qMTA) is used to integrate the organ-specific transcript abundances mapped to reactions and reference flux map and compute personalized organ-specific metabolic flux maps. The resulting flux maps can be used to perform fluxome-wide association analysis to complex traits or diseases. 

Genetically personalized flux maps can be computed running three scripts in succession:
 - gim3e_and_sampling.py : Runs GIM3E  to compute an organ-specific reference flux distribution from average organ transcript abundance.
 - map_expression_to_reactions.py :  Maps individual-level gene expression data as reaction activity fold changes relative to average transcript abundance. 
 - run_qMTA.py: Runs qMTA to compute the personalized flux distributions most consistent with the reaction activity fold changes starting from the reference flux distribution 
## Installation
The python scripts and underlying functions can be run in either python 2 or python 3. Source code is provided for either distribution in the directories python2 and  python3. Running “python setup.py install” will install the scripts and all their open source dependencies. Additionally, run_qMTA.py also requires the solver  CPLEX and its associated python package which must be installed separately. Cplex is freely available for academic use as part of the [IBM academic initiative](https://www.ibm.com/support/pages/ibm-ilog-optimization-academic-initiative) . Creating a virtualenv is advised to prevent compatibility issues with existing python installations.  
## Required inputs
 - Average organ gene expression: A CSV or XLSX file defining the average gene expression (in TPM, FPKM or equivalent units) in a given organ or tissue. Such data can be obtained from the GTEx database. The first column must define entrez ID, and subsequent columns should contain the average gene expression in each organ. 
 - Organ specific genome-scale metabolic models: Systems Biology Markup Language (SBML) models defining a metabolic network representing the metabolic potential of each organ of interest. Each model must be constrained to fulfil at least some of the metabolic function of the organ (e.g. ATP production and synthesis of neurotransmitter for Brain tissue) or else GIM3E will not generate meaningful results. The model should have a “.SBML” or “.XML” extension and be named either “organ_specific_metabolic_network_ORGAN_NAME” or “ORGAN_NAME” where ORGAN_NAME must match a column name in the average organ gene expression file. 
 - Organ-specific transcript abundance patterns imputed from genotype data. Such data can be derived by applying PredictDB or equivalent models on genotype data with appropriate tools (e.g. PLINK2 or Predict.py). Genes must be in rows and individuals in columns. The first column should be the Entrez ID.

Examples of inputs are provided in the directory inputs.
## gim3e_and_sampling.py
Average organ transcript abundances are mapped to the organ-specific subnetworks using the gene reaction annotations in each network. More in detail, transcripts abundances of isoenzymes and enzymes subunits catalysing each reaction or transport processes are added and, subsequently, LOG2 transformed. The resulting values are used as input to apply the GIM3E algorithm. The GIM3E algorithm applies a flux minimization weighted by transcript abundance allowing to identify solutions that are both enzymatically efficient and consistent with gene expression data. Subsequently, the flux ranges within 99% of the GIM3E optimal solution are identified and the resulting solution space is sampled using the Artificially Centered hit-and-run (ACHR) algorithm. The average of these flux samples can be used as the reference or average flux map of each organ.

### **Usage: 	gim3e_and_sampling.py [INPUTS...]** 

**INPUTS:**

*-m, --organ_specific_model_directory*  : Directory containing organ-specific genome-scale metabolic models in SBML format. All SBML models will be analyzed in succession. 

*-r, --reference_transcript_abundance* :  Path to the CSV or XLSX file defining the average organ gene expression in TPM or FPKM.

*-o, --output_directory* : Output directory. Will be created if it does not exist. 

**OUTPUTS:**
 
 - Reference_fluxes_ORGAN_NAME.csv : A CSV file with the GIM3E optimal solution, GIM3E flux ranges and average flux distributions for all reactions in all analysed organs. 
 - Reference_fluxes_ORGAN_NAME_dict.json : A json file with the average flux distributions for all reactions in all analysed organs. It is used as input for run_qMTA.py
 - ORGAN_NAME_gim3e__constrained_model.sbml : Organ-specific model constrained to the Gim3e flux ranges. Used as input for run_qMTA.py. 
## map_expression_to_reactions.py
This script maps genetically imputed patient-specific expression patterns to organ-specific models using the gene reaction annotations in these models. Imputed values are expressed as Log2 fold changes relative to average gene expression in a given organ and then mapped to reactions in the organ-specific model considering the relative transcript abundance of isoenzymes and enzyme subunits (e.g. in a reaction catalysed by multiple isoenzymes genetic variation on the isoenzyme with the highest expression will have a stronger effect on putative reaction activity). The script must be run for each organ under study. 
### Usage: map_expression_to_reactions.py [INPUTS...] 
**INPUTS:**

*-i, --imputed_transcript_abundance* : Path to the CSV or XLSX file with the Organ-specific transcript abundance patterns imputed from genotype data. Genes must be in rows and individuals/samples in columns. 

*-m, --organ_specific_model* : path to the organ-specific model in SBML format. Can also take the gim3e__constrained_model.sbml model as input. 

*-r, --reference_transcript_abundance* :  Path to the CSV or XLSX file defining the average organ gene expression in TPM or FPKM.

*-o, --output_directory* : Working and output directory. Will be created if it does not exist. 

*-s, --sample_list* : Optional, list of sample/individual IDs that should be analysed. If not provided all samples will be analysed. Each row should contain a sample/individuals ID.

*-t, --organ_name* :  Optional, Organ or tissue to be analysed. Has to match a column in the reference_transcript_abundance file. If not provided it will take organ name from the  organ_specific_model file name.

*-g, --gene_id_column_name* : Optional, defines the column name in imputed_transcript_abundance that defines the gene identifiers used in the model. If it is not provided, it will be set to "NCBI.gene..formerly.Entrezgene..ID". If it is not present it will be assumed to be the first column in the file. 


**OUTPUTS:**

- reaction_expression: CSV file containing putative reaction activity fold changes for each individual. Used as input for run_qMTA.py.
## run_qMTA.py
This script runs the quadratic metabolic transformation algorithm (qMTA). qMTA seeks to minimize the difference between the simulated flux distribution and the product of the putative fold changes by the reference flux distribution (target flux) while also minimizing the flux variation from the reference flux distribution in reactions without gene expression fold change. Additionally, both terms are scaled by the difference between the reference flux distribution and the target flux and the reference flux distribution, respectively, to prevent biases towards reactions with high reference flux values. Thus, qMTA can identify the flux map most consistent with gene expression fold changes starting from a reference flux distribution and compute personalized flux maps.
Usage run_qMTA.py [INPUTS...] 

**INPUTS:**

*-i, --imputed_reaction_fold_change* : Path to the reaction_expression CSV file containing the putative reaction activity fold changes for each individual generated by the map_expression_to_reactions.py script.

*-f, --reference_flux_json* : Path to the json file with the average flux distributions for all reactions in all analysed organs. It is generated by the gim3e_and_sampling.py script.

*-m, --organ_specific_gim3e_model* : Organ-specific model constrained to the GIM3E flux ranges generated by the gim3e_and_sampling.py script.

*-o, --working_directory* : Output and Working directory. Will be created if it does not exist.

*-a, --output_prefix* : Optional, prefix added to the output.

*-s, --sample_list* : Optional, list of sample/individual IDs that should be analysed. If not provided all samples will be analysed. Each row should contain a sample/individuals ID.

*-w, --min_flux_weight* : Optional, minimum value allowed for the flux normalization. Default is 1e-6. Lower values give more weight to reactions with low reference flux values.

*-c, --min_flux_fold_change* : Optional, minimum flux to consider reaction activity fold changes relevant. Putative reaction activity fold changes for reactions under this threshold will be ignored. Default is 1e-6.  

*-r, --unchanged_reaction_weight* : Optional, weight given to the minimization of reactions without imputed gene expression data. Default is 1. Lower values will allow more variation in such reactions.  

*-t, --organ_name* : Optional, Organ or tissue to be analysed. It should match  the organ name used in both gim3e_and_sampling.py and map_expression_to_reactions.py. If it is not provided, it will take the organ name from the organ_specific_gim3e_model file name. 

**OUTPUTS:**
 - Personalized fluxes file: CSV with the organ-specific personalized flux values computed for each individual. 
