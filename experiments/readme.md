## Folder structure of the experiment folders

- target_datetime/
  
  - hypothesis/
    - generated_hypotheses/
      - initial_stage/  # Contains hypotheses and experimental plans from the pattern selection step (number will be equal to N)
      - second_stage/  # Contains all of the generated variants (after pattern selection, M variants will be generated, best one is selected)
    - selected_hypothesis/
      - hypothesis_details.json  # basic information regarding the selected hypothesis
      - hypothesis.txt  # hypothesis in plaintext
      - reasoning.txt  # selection critertion for variant selection
      - variants.txt  # plaintext of all variants (should be same as the ones in second_stage)

  - protocol/
    - EVE/  # Contains the generated OVERLORD-protocol (not done yet)
    - hamilton/ # contains the files needed for the hamilton automation
     - channels/  # specifies activation and deactivation of hamilton-channels to avoid cross-contamination
         
    - mass_spectrometry/ # contains generated runlist and experimental template for AutonoMS
    - plate_layout/ # contains metadata for PLAID, and the generated layout
    - protocol.json # All experimental details regarding the experiment, extracted from selected hypothesis
      
  - results/ 
    - growth/
      - processed/  # Contains all the processed growth curves
      - raw/  # raw polarstar output
      - tests/  # contains significance testing data and metadata (growth-based)
        
    - metabolomics/ # Contains resulting transitionlist from AutonoMS (not completed)
    - plots/  # Generated plots
    - 
  - versions/  # Frozen settings, contexts and LLM-versions
  
