#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 17:47:07 2024

@author: danbru
"""
import os
import json
import re

from openai import OpenAI

import config

def prompt_gpt_for_hypothesis(context_path, prompted_statements, key, T ,llm_version='gpt-4o'):
    
    os.environ["OPENAI_API_KEY"] = key    
    client = OpenAI()
   
    with open(context_path, 'r') as file:
        # Read the entire file content
        context = file.read()
   
    completion = client.chat.completions.create(
        model = llm_version,
        temperature = T,
        messages=[
            {"role": "system", "content": context},
            {
                "role": "user",
                "content": prompted_statements
            }
        ]
    )
    return completion

def prompt_gpt_for_experimental_plan(context_path, hypothesis, key):
    
    os.environ["OPENAI_API_KEY"] = key
    client = OpenAI()

    with open(context_path, 'r') as file:
        # Read the entire file content
        exp_context = file.read()
    
    completion = client.chat.completions.create(
        model = config.LLM_AUTOFORMALIZATION_MODEL,
        temperature = config.LLM_AUTOFORMALIZATION_TEMPERATURE,
        messages=[
            {"role": "system", "content": exp_context},
            {
                "role": "user",
                "content": hypothesis
            }
        ]
    )
    
    # Save the model used for the prompt
    return completion

def is_valid_json(filename):    
    try:
        with open(filename, 'r') as json_file:
            json.load(json_file)  # Attempt to load the file
        return True  # If no exception is raised, it's valid
    except json.JSONDecodeError:
        return False  # If a JSONDecodeError is raised, it's not valid
    except Exception:
        return False  # Any other exception also indicates failure

def retry_experimental_plan(context_path, hypothesis, key, retries=3):

    attempts = 0
    success = False
    
    while attempts < retries and not success:
        completion = prompt_gpt_for_experimental_plan(context_path,hypothesis, key)
        experimental_plan = completion.choices[0].message.content
        json_string_cleaned = experimental_plan.replace("\\'", "'")
        
        # Sometimes the output is formatted as markdown, code below is to avoid any errors when compiling into json
        json_string_cleaned = json_string_cleaned.replace("'''", "")
        json_string_cleaned = json_string_cleaned.replace("json", "")
        json_string_cleaned = json_string_cleaned.replace("```", "")
        
        try:
            # Attempt to load the JSON
            data = json.loads(json_string_cleaned)
            success = True  # Mark success to exit the loop
            return data, experimental_plan, completion
            
        except json.JSONDecodeError as e:
            print(json_string_cleaned)
            attempts += 1
            print(f"Failed to load JSON on attempt {attempts}: {e}")
            if attempts < retries:
                print(f"Retrying... ({retries - attempts} attempts left)")
            else:
                print("Maximum retry limit reached. Exiting. Please regenerate hypothesis.")
                exit()

def extract_statements(text, numbers):
    
    # Create a regex pattern to match each numbered statement
    pattern = r"(\d+\.\s[\s\S]*?(?=\d+\.\s|$))"
    
    # Find all statements in the text that match the pattern
    all_statements = re.findall(pattern, text)

    # Filter out statements based on the given list of numbers
    extracted_statements = [statement for statement in all_statements if int(statement.split('.')[0]) in numbers]

    return extracted_statements

