# Agentic AI Integrated with Scientific Knowledge: Laboratory Validation in Systems Biology

## Install and setup

### Install Python and R dependencies

Using Python (3.10.13) and R (4.2.3) the dependencies can be installed from the requirements.txt file, e.g. using conda and the following commands:
```
$ conda create --name genExp python=3.10.13 r-base=4.2.3 && \\
    conda activate genExp && \\
    pip install -r requirements.txt
```
### Install SWI-Prolog

Follow download and install instructions [here](https://www.swi-prolog.org/download/stable). It can also be installed using package managers such as apt, snap, and brew.

### Install RMLMapper

Run `scripts/install_rml.sh` to install the RMLMapper JAR in `/opt/tools`. If you decide to install this somewhere else in your filesystem, you will need to change the `RMLMAPPER_JAR` variable in `map_protocols.sh`.

### Setting up a ChEBI SPARQL endpoint

Currently the database creation programmes rely on a privately hosted SPARQL endpoint for ChEBI to look up compounds for inclusion in the graph database.

We recommend these steps.

1. Follow the instructions for our (recently tested) [containerised Fuseki server implementation](https://github.com/TW-Genesis/genesis-database-system) to build the Docker image locally.
2. Download the `chebi.ttl` file from the Zenodo store and save in the `data/` directory.
3. Again following the instructions, run ```DATA=chebi.ttl docker compose up load-data```.
4. Once the data loading is complete, run ```docker compose up start-server -d```.

We are investigating publically hosted options that have the same functionality, which should simplify this process.

### Setting up an API-key

Add a key.txt file (see .gitignore) containing only the API key ("sk-proj---XXXXXXXX...") for OpenAI in the root folder.

## Usage

### Prompts

The template prompts used for the examples in the manuscript are available in `/context`. Investigation-specific contexts/prompts (i.e. the ones actually used for the different investigations) can be found inside the experiments folder. e.g. `experiments/completed_experiments/arginine_202503131539/versions`.

### Generation and experiment execution parameters

Settings regarding the runs can be found in `scripts/config.py`.
Note that you will need to change paths to relevant executables in the config-file.

### Generating an hypothesis, experimental design and liquid handling scripts

From the `/script` folder, run the following command in the terminal (fill in the blanks):

```
python genexp.py --target <string> --N <integer> --alpha <float>
```

- `target` denotes the metabolite observable used for learning the association (an amino acid, in this case).
- `N` denotes the number of patterns passed to the hypothesis generation step. If more than 1, an LLM-agent will select the most reasonable one. Default at 1.
- `alpha` is a float between 0 and 1 that is used to penalize patterns not unique to the specific metabolite observable (a number closer to 1 will ensure that patterns that are only deemed important for your specific target will rank higher). Default at 0.0
- `override` is an optional parameter, if you want to manually override the selected logic program with one of your own.

This will create a folder in `/experiments` with all of the details regarding the experiments (e.g. hypothesis, protocol, liquid handling scripts, ...). Note that you will be prompted for stock concentrations (if compounds are not present in the library) during the run. When the scripts have been run on the Hamilton, EVE and via AutonoMS [2] and data has been aquired and saved in `data/growth/raw` and `data/metabolomics/raw`, run the following command to process and analyse all of the data. For details regarding data acquisition, see `protocol/hamilton/scripts`, `protocol/overlord/scripts` and `protocol/mass_spectrometry`. `output_folder` denotes the folder created by the prior step. 

```
python analysis.py --output_folder <string> --metabolomics_analysis <bool>
```
This will automatically run outlier curation, processing, normalization and statistical testing on the growth data and metabolomics data. It will also generate a basic result report.

## References:
1. https://www.sciencedirect.com/science/article/pii/S266731852300017X?via%3Dihub
2. https://pubs.acs.org/doi/10.1021/jasms.3c00396
