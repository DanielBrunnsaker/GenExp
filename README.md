# GenExp
 
Run h_generator.py from the scripts folder. e.g. 

python h_generator.py --target "alanine" --N 10 --T 0.5 --alpha 1.0 --beta 0.25

 - target: Amino acid to base the hypothesis around (i.e. the coefficients to use)
 - N: number of patterns to present to the LLM. Tested with 10, performance decreases with more, and a bit boring with less than five.
 - T: temperature of the hypothesis generation step. Best results come with the other steps (selection and protocol) fixed with a low temperature.
 - alpha/beta: alpha * target_values (normalized coefficients for selected target) - beta * penalty_values (where the penalty is the sum of the normalized coefficients of the row, e.g. all amino acids except target). Makes it so that the method will prioritize patterns that are more specific to one amino acid. Otherwise the conclusions will be quite unclear (can also be adapted so that we take more amino acids into account, but complicates the hypothesis, and i do not think it is needed for a proof of concept).


Outputs should then be saved to the experiments/generated_outputs folder.

 - selection (json): Filters down the given patterns to a set of five, based on feasibility and safety, given a rough explanation of our setting. The output here is a JSON-file with the selected pattern(s) and a brief summary.
 - hypothesis: Main text, contains the hypothesis and experimental plan. Output is a text file.
 - protocol: Summarized experimental protocol in JSON.


