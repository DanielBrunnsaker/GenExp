#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 17:47:07 2024

@author: danbru
"""

def prompt_gpt_for_hypothesis(prompted_statements, temp, key):
    
    from openai import OpenAI
    import os
    
    os.environ["OPENAI_API_KEY"] = key
    print('Generating hypotheses... \n')
    
    client = OpenAI()
    
    # Open a text file in read mode
    with open('../context/hypgen_context.txt', 'r') as file:
        # Read the entire file content
        context = file.read()
    
    #prompted_statements = "1. Cells with higher than normal levels of intracellular proline in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: Cell(A):-phenotype(A,'decreased resistance to chemicals',B),compound_name(B,niclosamide) 2. Cells with higher than normal levels of intracellular proline in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: Cell(A):-phenotype(A,'decreased resistance to chemicals',B),compound_name(B,'5-[[3-(1-phenylethoxy)-4-(2-phenylethoxy)phenyl]methylene]-4-oxo-2-thioxo-3-thiazolidineacetic acid') 3. Cells with higher than normal levels of intracellular proline in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: Cell(A):-phenotype(A,'increased chemical compound accumulation',B),compound_name(B,'lithium(1+)')"
    
    completion = client.chat.completions.create(
        #model="gpt-4o-mini",
        model="gpt-4o",
        temperature = temp,
        messages=[
            {"role": "system", "content": context},
            {
                "role": "user",
                "content": prompted_statements
            }
        ]
    )
    
    message = completion.choices[0].message
    hypothesis = message.content
    
    return hypothesis

def prompt_gpt_for_experimental_plan(hypothesis, key):
    
    from openai import OpenAI
    import os
    
    os.environ["OPENAI_API_KEY"] = key

    client = OpenAI()
    
    print('Designing experimental plan... \n')
    # Redefine a new one.
    # Open a text file in read mode
    with open('../context/expdesign_context.txt', 'r') as file:
        # Read the entire file content
        exp_context = file.read()
    
    completion = client.chat.completions.create(
        #model="gpt-4o-mini",
        model="gpt-4o",
        temperature = 0.0,
        messages=[
            {"role": "system", "content": exp_context},
            {
                "role": "user",
                "content": hypothesis
            }
        ]
    )
    
    message = completion.choices[0].message
    experimental_plan = message.content
    return experimental_plan

def is_valid_json(filename):
    
    import json
    
    try:
        with open(filename, 'r') as json_file:
            json.load(json_file)  # Attempt to load the file
        return True  # If no exception is raised, it's valid
    except json.JSONDecodeError:
        return False  # If a JSONDecodeError is raised, it's not valid
    except Exception:
        return False  # Any other exception also indicates failure

def retry_experimental_plan(hypothesis, key, retries=3):
    
    import json
    
    attempts = 0
    success = False
    
    while attempts < retries and not success:
        experimental_plan = prompt_gpt_for_experimental_plan(hypothesis, key)
        json_string_cleaned = experimental_plan.replace("\\'", "'")
        
        try:
            # Attempt to load the JSON
            data = json.loads(json_string_cleaned)
            
            print("JSON saved successfully!")
            success = True  # Mark success to exit the loop
            
            return data, experimental_plan
            
        except json.JSONDecodeError as e:
            attempts += 1
            print(f"Failed to load JSON on attempt {attempts}: {e}")
            if attempts < retries:
                print(f"Retrying... ({retries - attempts} attempts left)")
            else:
                print("Maximum retry limit reached. Exiting.")
                return None
            
            
def safety_feasibility_prompt(full_prompt, key):
 
     from openai import OpenAI
     import os
     
     os.environ["OPENAI_API_KEY"] = key
     print('Selecting clauses based on safety and feasibility... \n')
     
     client = OpenAI()
     
     # Open a text file in read mode
     with open('../context/feasibility_context.txt', 'r') as file:
         # Read the entire file content
         context = file.read()
     
     #prompted_statements = "1. Cells with higher than normal levels of intracellular proline in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: Cell(A):-phenotype(A,'decreased resistance to chemicals',B),compound_name(B,niclosamide) 2. Cells with higher than normal levels of intracellular proline in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: Cell(A):-phenotype(A,'decreased resistance to chemicals',B),compound_name(B,'5-[[3-(1-phenylethoxy)-4-(2-phenylethoxy)phenyl]methylene]-4-oxo-2-thioxo-3-thiazolidineacetic acid') 3. Cells with higher than normal levels of intracellular proline in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: Cell(A):-phenotype(A,'increased chemical compound accumulation',B),compound_name(B,'lithium(1+)')"
     
     completion = client.chat.completions.create(
         #model="gpt-4o-mini",
         model="gpt-4o",
         temperature = 0.0,
         messages=[
             {"role": "system", "content": context},
             {
                 "role": "user",
                 "content": full_prompt
             }
         ]
     )
     
     message = completion.choices[0].message
     hypothesis = message.content
     
     #print(hypothesis)
 
     return hypothesis

def extract_statements(text, numbers):
    
    import re
    
    # Create a regex pattern to match each numbered statement
    pattern = r"(\d+\.\s[\s\S]*?(?=\d+\.\s|$))"
    
    # Find all statements in the text that match the pattern
    all_statements = re.findall(pattern, text)

    # Filter out statements based on the given list of numbers
    extracted_statements = [statement for statement in all_statements if int(statement.split('.')[0]) in numbers]

    return extracted_statements