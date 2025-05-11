#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:26:00 2024

@author: danbru
"""

import re
import json
from utils.hgen_support import *
from utils.gpt_support import *
import config
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def hypothesis_generation(output_folder,  N):

    # T = config.LLM_GENERATION_TEMPERATURE
    variations = config.VARIATIONS

    # Load OpenAI-key
    key = open(BASE_DIR / "../key.txt", "r")
    key = key.read()

    EXPERIMENT_DIR = Path(output_folder).resolve()

    '''
    Select a logic program (and its associated experimental plan and hypothesis) to work with
    '''

    if int(N) > 1:
        complete_prompt = ''
        for nr_hyp in range(1, N+1):
            complete_prompt += f'####### GENERATED HYPOTHESIS {nr_hyp} #######\n'
            hyp_text = open(
                EXPERIMENT_DIR / f"hypothesis/generated_hypotheses/initial_stage/hypothesis_{nr_hyp}.txt", "r").read()
            complete_prompt += hyp_text + '\n\n'

        completion = prompt_gpt_for_hypothesis(EXPERIMENT_DIR / '../../context/hypothesis_selection.txt',
                                               complete_prompt, key, config.LLM_PATTERN_SELECTION_TEMPERATURE, config.LLM_PATTERN_SELECTION_MODEL)
        final_selection_text = completion.choices[0].message.content

        # Extract the hypothesis number
        N_extracted = int(
            re.search(r"HYPOTHESIS NUMBER #(\d+)", final_selection_text).group(1))

        selected_hypothesis = open(
            EXPERIMENT_DIR / f"hypothesis/generated_hypotheses/initial_stage/hypothesis_{N_extracted}.txt", "r").read()
        hypothesis_details = load_json(
            EXPERIMENT_DIR / f"hypothesis/generated_hypotheses/initial_stage/hypothesis_details_{N_extracted}.json")

        # Add the motivation behind selecting it
        hypothesis_details['motivation'] = final_selection_text
        with open(EXPERIMENT_DIR / 'hypothesis/selected_hypothesis/hypothesis_details.json', 'w') as file:
            json.dump(hypothesis_details, file, indent=4)

    else:
        selected_hypothesis = open(
            EXPERIMENT_DIR / f"hypothesis/generated_hypotheses/initial_stage/hypothesis_{N}.txt", "r").read()
        hypothesis_details = load_json(
            EXPERIMENT_DIR / f"hypothesis/generated_hypotheses/initial_stage/hypothesis_details_{N}.json")

        # Add the motivation behind selecting it
        hypothesis_details['motivation'] = 'not applicable'
        with open(EXPERIMENT_DIR / 'hypothesis/selected_hypothesis/hypothesis_details.json', 'w') as file:
            json.dump(hypothesis_details, file, indent=4)

    '''
    Generates additional experimental plan variations
    '''

    variants = '####### OPTION 1 #######\n\n' + selected_hypothesis
    with open(EXPERIMENT_DIR / 'hypothesis/generated_hypotheses/second_stage/hypothesis_1.txt', 'w') as file:
        file.write(selected_hypothesis)

    for variants_iter in range(2, (variations+1)):
        variants += f'\n\n####### OPTION {variants_iter} #######\n'
        completion = prompt_gpt_for_hypothesis(EXPERIMENT_DIR / '../../context/hypgen_per_pattern_temp.txt', hypothesis_details['initial_prompt'],
                                               key, config.LLM_VARIATION_TEMPERATURE, config.LLM_VARIATION_MODEL)
        new_variant = completion.choices[0].message.content

        with open(EXPERIMENT_DIR / f'hypothesis/generated_hypotheses/second_stage/hypothesis_{variants_iter}.txt', 'w') as file:
            file.write(new_variant)

        variants += new_variant

    with open(Path(output_folder) / "versions/v_variation.txt", "w") as file:
        file.write(completion.model)

    with open(EXPERIMENT_DIR / 'hypothesis/selected_hypothesis/variants.txt', 'w') as file:
        file.write(variants)

    completion = prompt_gpt_for_hypothesis(EXPERIMENT_DIR / '../../context/refine_plan.txt', variants, key,
                                           config.LLM_EXPERIMENT_SELECTION_TEMPERATURE, config.LLM_EXPERIMENT_SELECTION_MODEL)

    with open(Path(output_folder) / "versions/v_final_selection.txt", "w") as file:
        file.write(completion.model)

    final_selection_text = completion.choices[0].message.content

    with open(EXPERIMENT_DIR / 'hypothesis/selected_hypothesis/reasoning.txt', 'w') as file:
        file.write(final_selection_text)

    N_extracted = int(re.search(r"HYPOTHESIS NUMBER #(\d+)",
                      final_selection_text).group(1))
    finalized_hypothesis = open(
        EXPERIMENT_DIR / f"hypothesis/generated_hypotheses/second_stage/hypothesis_{N_extracted}.txt", "r").read()

    with open(EXPERIMENT_DIR / 'hypothesis/selected_hypothesis/hypothesis.txt', 'w') as file:
        file.write(finalized_hypothesis)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Hypothesis generation")
    parser.add_argument("--output_folder", required=True,
                        type=str, help="Experiment folder")
    parser.add_argument("--N", required=True, type=int,
                        help="Number of hypotheses to select")

    args = parser.parse_args()
    hypothesis_generation(args.output_folder,  args.N)
