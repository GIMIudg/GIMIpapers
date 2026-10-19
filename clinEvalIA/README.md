# ClinEval

ClinEval is a pipeline for the systematic evaluation and refinement of structured clinical instruments using Large Language Models (LLMs).

## Overview

This project implements an automated framework to assess and improve the structural quality of a neurological clinical history form. The pipeline uses an LLM as an evaluator (“LLM-as-judge”) to score each question across four dimensions:

- Clarity  
- Clinical relevance  
- Completeness  
- Granularity  

Based on these evaluations, targeted refinements are applied and re-evaluated iteratively.

## Methodology

1. **Schema Extraction**  
   The clinical form (XLSXForm) is parsed into a structured JSON schema.

2. **LLM-based Evaluation**  
   Each question is evaluated using a rubric-driven prompt with multiple independent runs to ensure stability.

3. **Iterative Refinement**  
   Low-performing sections and questions are identified and improved through controlled modifications.

4. **Re-evaluation**  
   The refined instrument is re-assessed to quantify improvements.

## Key Features

- Structured, reproducible evaluation pipeline  
- Multi-run stability analysis  
- Targeted refinement without altering instrument structure  
- Integration of standardized clinical vocabularies (ICD, LOINC, RxNorm)

## Results

The pipeline demonstrates measurable improvements in structural quality, particularly in granularity and completeness, while maintaining stability across evaluations.

## Limitations

- Dependence on a single LLM for evaluation  
- Partial reliance on human interpretation for applying modifications  
- Limited external validation with clinical experts

## Future Work

- Cross-model evaluation (multiple LLMs)  
- Full automation of the refinement cycle  
- Expanded clinical validation  

## Author

Department of Translational Bioengineering  
University of Guadalajara