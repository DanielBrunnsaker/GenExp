## Folder structure of the experiment folders

Alec & Filip, does this seem like a reasonable structure? Feel free to make subfolders in the experiments for your stuff

- target_datetime/
  
  - hypothesis/
    - generated_hypotheses/
      - initial_stage/ # Contains hypotheses and experimental plans from the pattern selection step (number will be equal to N)
      - second_stage/ # Contains all of the generated variants (after pattern selection, M variants will be generated, best one is selected)
    - selected_hypothesis/
      - hypothesis_details.json # basic information regarding the selected hypothesis
      - hypothesis.txt # hypothesis in plaintext
      - reasoning.txt # selection critertion for variant selection
      - variants.txt # plaintext of all variants (should be same as the ones in second_stage)

  - protocol/
    - EVE/
    - hamilton/
    - mass_spectrometry/
    - plate_layout/
    - protocol.json
      
  - results/
    - growth/
    - metabolomics/
    - plots/
    - 
  - versions/
  
