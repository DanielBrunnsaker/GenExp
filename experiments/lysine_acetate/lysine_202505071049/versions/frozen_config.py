
# config.py

# Workflow settings
PRINT_HYPOTHESIS_DURING_WORKFLOW = True
ALLOW_FOR_CHOICE = True

# Logic program prioritization settings
REUSE_LOGIC_PROGRAM = False # Below logic can be changed if you want to bias the experiments differently
RELEVANCE_SCORES = {
        # predicates
        'compound_name': 5, # Reward patterns with a specific condition in mind (if there is one)
        'condition': 10, # Reward patterns with a specific chemical in mind (if there is one)
        'participates_in_metabolism': 5, 
        'compound_modulates_target': 5, 
        'metabolism of': 2
        #'production of': 2,
        #'consumption of': 2,
        #'metabolism of': 1,
    }



# Dispensing settings
WELL_VOLUME = 270 # in uL
WELL_VOLUME_UNIT = 'microliter'
DILUTION_FACTOR = 15 # Dilution factor for preculture

# Layout generation settings
MINIZINC_PATH = "/Applications/MiniZincIDE.app/Contents/Resources/minizinc"

# Generation settings
LLM_PATTERN_SELECTION_MODEL = 'gpt-4o'
LLM_VARIATION_MODEL = 'gpt-4o'
LLM_EXPERIMENT_SELECTION_MODEL = 'gpt-4o'
LLM_AUTOFORMALIZATION_MODEL = 'gpt-4o'

LLM_PATTERN_SELECTION_TEMPERATURE = 0.5
LLM_VARIATION_TEMPERATURE = 0.5
LLM_EXPERIMENT_SELECTION_TEMPERATURE = 0.5
LLM_AUTOFORMALIZATION_TEMPERATURE = 0.05
VARIATIONS = 5 # How many variations to generate from the pattern selection step

# Mass spectrometry settings
PERFORM_MASS_SPEC = True
SAMPLES_IN_SEQUENCE = 6 # How many samples to include in sample blocks
BLANK_BETWEEN_SAMPLE_BLOCKS = 4
TUNE_INTERVAL = 0 # How many sample-blocks (blank, qc, sample, qc, blank) to have inbetween reference tuning injections
RF_SETTINGS_FILE = '../data/rf_params/AutonoMS_template.xlsx'
RUNLIST_RANDOMIZATION = True # If one wants to randomize the runlit. Recommended!
MS_POLARITY = 'Positive'
METHOD_NAME = '2024-10-10_AminoAcids_test_1_Pos.m' # Mass spec method
COLUMN_TYPE = 'H'
PLATE_TYPE = 'P96'


# Growth analysis settings
BLANK_SUBTRACTION = True
N_CLOSEST_BLANKS = 3 # How many of the closest blanks to use for blank-normalization in OD processing
BLANK_FILLIN_VALUE = 0.01 # If negative value after blank subtraction, what value to set to instead

OUTLIER_METHOD = 'mad'
ROLLING_MEAN_WINDOWSIZE = 3 # Windowsize for smoothing
THRESHOLD = 3.5 # MAD-threshold to use for curve fitering. #3 for MAD, 1.5 for IQR?



