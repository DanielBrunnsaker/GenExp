# GenExp
 
Run h_generator.py from the scripts folder. e.g. 

python h_generator.py --target "alanine" --N 10 --T 0.5 --alpha 1.0 --beta 0.25

Outputs should then be saved to the experiments/generated_outputs folder.

 - selection (json): Filters down the given patterns to a set of five, based on feasibility and safety, given a rough explanation of our setting. The output here is a JSON-file with the selected pattern(s) and a brief summary.
 - hypothesis: Main text, contains the hypothesis and experimental plan. Output is a text file.
 - protocol: Summarized experimental protocol in JSON.


