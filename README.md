# GenExp

## Install Python dependencies

Using Python 3.10.13 the dependencies can be installed from the requirements.txt file, e.g. using conda and the following commands:
```
$ conda create --name genExp python=3.10.13 && \\
    conda activate genExp && \\
    pip install -r requirements.txt
```
## Install SWI-Prolog

Follow download and install instructions [here](https://www.swi-prolog.org/download/stable). It can also be installed using package managers such as apt, snap, and brew. For more instructions on how to generate the patterns used for the hypothesis generation steps, see the `/prolog` folder.

## Setting up an API-key

Add a key.txt file (see .gitignore) containing only the API key ("sk-proj---XXXXXXXX...") for GPT4 in the root folder. Easiest to generate personal one (I can also give you one if you prefer to use the same).

## Generation and experiment execution parameters

Settings regarding the runs can be found in `scripts/config.py`.

## Generating an hypothesis, experimental design and liquid handling scripts

From the `/script` folder, run the following command in the terminal (fill in the blanks):

```
python genexp.py --target <string> --N <integer> --alpha <float>
```

- `target` denotes the metabolite observable used for the implication (an amino acid, in this case).
- `N` denotes the number of patterns passed to the hypothesis generation step (a higher number will allow for more variance, but lower ranked patterns are less likely to be true). Default at 10.
- `alpha` is a float between 0 and 1 that is used to penalize patterns not unique to the specific metabolite observable (a number closer to 1 will ensure that patterns that are only deemed important for your specific target will rank higher). Default at 0.0


This will create a folder in `/experiments` with all of the details regarding the experiments (e.g. hypothesis, protocol, liquid handling scripts, ...). Note that you will be prompted for stock concentrations (if applicable) during the run.

## TODO:
- Potentially generate the Overlord-protocol (e.g. replace some variables in the reference protocol such as sampling time and regenerate)
- Need to generate a SPE-IMS runlist using the plate-layout (connected to AutonoMS) [FINISHED]
- Still need the AutonoMS parts [2] (should be somewhat trivial, as AutonoMS outputs an annotated tsv with peak areas)
- Hypothesis formalization (alec & filip)
- Automated testing (daniel & alec & filip) 

### References:
1. https://www.sciencedirect.com/science/article/pii/S266731852300017X?via%3Dihub
2. https://pubs.acs.org/doi/10.1021/jasms.3c00396
