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

## Generating an hypothesis, experimental design and liquid handling scripts

From the `/script` folder, run the following command in the terminal (fill in the blanks):

```
python hgen.py --target <string> --N <integer> --alpha <float> --volume <integer> --T <float>
```

- `target` denotes the metabolite observable used for the implication (an amino acid, in this case).
- `N` denotes the number of patterns passed to the hypothesis generation step (a higher number will allow for more variance, but lower ranked patterns are less likely to be true). Default at 10.
- `alpha` is a float between 0 and 1 that is used to penalize patterns not unique to the specific metabolite observable (a number closer to 1 will ensure that patterns that are only deemed important for your specific target will rank higher). Default at 0.1
- `volume` denotes the final cultivation volume in microliters (e.g. 225 uL for a regular 96 well plate).
- `T` denotes the temperature in the prompting step. Keep between 0 and 1. (0 more deterministic, 1 less deterministic). Reasonable default around 0.3-0.7.


This will create a folder in `/experiments` with all of the details regarding the experiments (e.g. hypothesis, protocol, liquid handling scripts, ...). Note that you will be prompted for stock concentrations (if applicable) during the run.

<br/>

Alternatively you can run the steps in sequence:

```
python pattern_selector.py --target <string> --N <integer> --alpha <float>
```
This extracts several pattern (given the `target` and `N`), and passes it as an initial prompt to the LLM. The expected return is a list of feasible patterns, along with a short description of their relevance. It will also create the folder structure for the project. Additional outputs are the selection summaries (`hypothesis/feasibility/selection.json`) and the given prompt (`hypothesis/patterns/prompt.txt`)

```
python hypothesis_generation.py --T <string> --folder <path>
```
Generates the hypothesis, given the path to the previously generated experiment folder. `T` denotes the temperature of the LLM. Output is the generated hypothesis text (`hypothesis/generated_hypothesis.txt`)

```
python autoformalize_protocol.py --folder <path>
```
Autoformalizes parts of the hypothesis and rough experimental protocol into a JSON-file containing all the needed parameters for subsequent automation. Output is the formalized protocol in JSON-format (`protocol/protocol.json`).

```
python plate_layout.py --folder <path>
```
Reads the formalized protocol and generates a robust plate-layout with constraint programming, using PLAID [1]. Outputs are the plate layout in a tabular format (`protocol/plate_layout/plate_layout.tsv`) and the JSON file used to generate the layout in minizinc (`protocol/plate_layout/minzinc_reference.json`). Note that this currently only works with a 96 well format (might fix at some point).

```
python hamilton_protocol.py --folder <path> --volume <integer> --S <string> --Treatment <string>
```
Uses the plate layout and the formalized protocol design a runlist for a Hamilton Microlab Star. `S` denotes stock concentration of the media supplement (typically an amino acid, but can differ). `Treatment` denotes the stock concentration of the chemical treatment (if applicable). Any unit of concentration (e.g. mM (preferrable), mg/ml or % (v/v) should work). These, along with the final well `volume`, are needed to calculate the exact amount of volume to dispense in the plate preparation steps. Outputs include the the experimental layout along with a summary of each variable combination (`protocol/hamilton/pipetting_layout.xlsx`), a Hamilton Microlab Star runlist (`protocol/hamilton/runlist.xlsx`) and instructions for detailed liquid channel usage (`protocol/hamilton/channels/`).


## TODO:
- Potentially generate the Overlord-protocol (e.g. replace some variables in the reference protocol such as sampling time and regenerate)
- Need to generate a SPE-IMS runlist using the plate-layout (connected to AutonoMS)
- Still need the AutonoMS parts [2] (should be somewhat trivial, as AutonoMS outputs an annotated tsv with peak areas)
- Hypothesis formalization (alec & filip)
- Automated testing (daniel & alec & filip)

### References:
1. https://www.sciencedirect.com/science/article/pii/S266731852300017X?via%3Dihub
2. https://pubs.acs.org/doi/10.1021/jasms.3c00396
