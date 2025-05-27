# GenExp

## Install Python dependencies

Using Python 3.10.13 the dependencies can be installed from the requirements.txt file, e.g. using conda and the following commands:
```
$ conda create --name genExp python=3.10.13 r-base=4.2.3 && \\
    conda activate genExp && \\
    pip install -r requirements.txt
```
## Install SWI-Prolog

Follow download and install instructions [here](https://www.swi-prolog.org/download/stable). It can also be installed using package managers such as apt, snap, and brew. For more instructions on how to generate the patterns used for the hypothesis generation steps, see the `/prolog` folder.

## Setting up an API-key

Add a key.txt file (see .gitignore) containing only the API key ("sk-proj---XXXXXXXX...") for GPT4 in the root folder. Easiest to generate personal one (I can also give you one if you prefer to use the same).

## Generation and experiment execution parameters

Settings regarding the runs can be found in `scripts/config.py`.
Note that you will need to change paths to relevant executables in the config-file.

## Generating an hypothesis, experimental design and liquid handling scripts

From the `/script` folder, run the following command in the terminal (fill in the blanks):

```
python genexp.py --target <string> --N <integer> --alpha <float>
```

- `target` denotes the metabolite observable used for the implication (an amino acid, in this case).
- `N` denotes the number of patterns passed to the hypothesis generation step (a higher number will allow for more variance, but lower ranked patterns are less likely to be true). Default at 10.
- `alpha` is a float between 0 and 1 that is used to penalize patterns not unique to the specific metabolite observable (a number closer to 1 will ensure that patterns that are only deemed important for your specific target will rank higher). Default at 0.0
- `override` is an optional parameter, if you want to manually override the selected logic program with one of your own.

This will create a folder in `/experiments` with all of the details regarding the experiments (e.g. hypothesis, protocol, liquid handling scripts, ...). Note that you will be prompted for stock concentrations (if compounds are not present in the library) during the run. 

When the scripts have been run on the Hamilton, EVE and via AutonoMS [2] and data has been aquired and saved in `data/growth/raw` and `data/metabolomics/raw`, run the following command to process and analyse all of the data. For details regarding data acquisition, see `protocol/hamilton/scripts`, `protocol/overlord/scripts` and `protocol/mass_spectrometry`.

```
python analysis.py --output_folder <string> --metabolomics_analysis <bool>
```
This will automatically run outlier curation, processing, normalization and statistical testing on the growth data and metabolomics data. 

## TODO:
- Database integration scripts (alec)
- Database querying/resue (alec)

### References:
1. https://www.sciencedirect.com/science/article/pii/S266731852300017X?via%3Dihub
2. https://pubs.acs.org/doi/10.1021/jasms.3c00396
