import json
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from utils.hgen_support import *
from utils.gpt_support import *
from utils.selection_utils import *
import config

# Ensure relative paths work correctly
BASE_DIR = Path(__file__).resolve().parent

# Load disallowed clauses
with open(BASE_DIR / '../experiments/forbidden_patterns.txt', "r") as file:
    disallowed_clauses = [line.strip() for line in file.readlines()]

with open(BASE_DIR / '../experiments/succesful_patterns.txt', "r") as file:
    successful_clauses = [line.strip() for line in file.readlines()]

disallowed_clauses += successful_clauses


def pattern_selection(target, alpha, N, output_folder, override, override_negative):

    # Load supporting data
    id_dict, orf_dict = load_supporting_data(BASE_DIR)

    # Load reference data
    aaSet = pd.read_excel(BASE_DIR / '../data/AA.xls',
                          sheet_name='intracellular_concentration_mM').set_index('ORF').iloc[:, 1:]

    # Path for coefficients (update based on latest results)
    coefficients_path = (
        Path(BASE_DIR) / Path('../results/coefficients/20241114')).resolve()
    hypotheses, all_hypotheses = load_and_filter_hypotheses(
        coefficients_path, target, aaSet, BASE_DIR)

    # Filter predictable amino acids
    metrics = pd.read_csv(
        BASE_DIR / '../results/metrics/eCV_results_20241114.csv').groupby('Amino acid').mean()
    all_hypotheses = all_hypotheses[metrics.index[metrics['R2'] > 0]]

    # Rank hypotheses
    filtered_sorted_specifity = rank_hypotheses(all_hypotheses, target, alpha)

    hypothesis_counter = 0
    for counter, hypothesis in enumerate(filtered_sorted_specifity.index):
        directionality = 'higher' if np.sign(
            hypotheses.loc[hypothesis]) == 1.0 else 'lower'

        with open(BASE_DIR / '../results/patterns/feature_dict.pickle', 'rb') as handle:
            duplicate_dict = pickle.load(handle)

        current_clause = pick_best_clause(hypothesis, duplicate_dict, BASE_DIR)
        clause = prepare_clause_text(current_clause, BASE_DIR)

        # Find negative control examples
        if override_negative:
            neg_examples = [override_negative]
        else:
            neg_examples = list(
                all_hypotheses.loc[hypothesis][all_hypotheses.loc[hypothesis] == 0].index)

        if neg_examples:
            add_negative_examples(
                BASE_DIR / '../context/hypgen_per_pattern.txt',
                BASE_DIR / '../context/hypgen_per_pattern_temp.txt',
                neg_examples
            )
        else:
            continue

        # Check for disallowed hypotheses
        if any(f'{target}:'+clause.replace("'", "") in s for s in disallowed_clauses):
            continue

        if override:
            # prompt_part = override
            ind_statement = override
        else:
            prompt_part = f"{counter}. Cells with {directionality} levels of intracellular {target} in standard conditions associate with:"
            ind_statement = prompt_part + clause
        # Generate hypothesis text
        hypothesis_counter += 1

        completion = prompt_gpt_for_hypothesis(
            BASE_DIR / '../context/hypgen_per_pattern_temp.txt',
            ind_statement,
            open(BASE_DIR / "../key.txt", "r").read(),
            config.LLM_PATTERN_SELECTION_TEMPERATURE,
            config.LLM_PATTERN_SELECTION_MODEL
        )

        hypothesis_text = completion.choices[0].message.content
        with open(Path(output_folder) / "versions/v_hypothesis_selection.txt", "w") as file:
            file.write(completion.model)

        save_hypothesis_details(clause.replace("'", ""), target, directionality,
                                filtered_sorted_specifity[hypothesis], ind_statement, output_folder, hypothesis_counter)

        with open(Path(output_folder).resolve() / 'hypothesis/generated_hypotheses/initial_stage' / f'hypothesis_{hypothesis_counter}.txt', 'w') as file:
            file.write(hypothesis_text)

        # Stop when reaching the desired number of hypotheses
        if hypothesis_counter == int(N):
            break


# Run as standalone script
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Pattern Selection")
    parser.add_argument("--target", required=True, type=str,
                        help="Target molecule (e.g., alanine)")
    parser.add_argument("--alpha", required=False, type=float,
                        default=0.0, help="Weight for pattern ranking")
    parser.add_argument("--N", required=True, type=int,
                        help="Number of hypotheses to select")
    parser.add_argument("--override", required=False, type=str,
                        default=None, help="Override clause.")
    parser.add_argument("--override_negative", required=False,
                        type=str, default=None, help="Override negative control.")

    args = parser.parse_args()

    # Call function with parsed arguments
    pattern_selection(args.target, args.alpha, args.N,
                      args.override, args.override_negative)
