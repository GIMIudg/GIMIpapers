
import pandas as pd
import mygene
import csv
import os


# Function to print the header display
def print_header():
    """
    Prints a header display for the XomicsGeneratorTool.

    Returns:
        None

    Displays the tool's name, a brief description of its purpose, and formatting lines for visual clarity.
    """
    print("="*60)
    print("XomicsGeneratorTool")
    print("="*60)
    print("This small tool takes a generic expression matrix (mxn) as an argument")
    print("to return individual files that can be used by the tool")
    print('"XomicsToModel" in MATLAB for the generation of genome-scale models.')
    print("="*60)
    print()


# Function to print the steps
def print_step(step_description):
    """
    Prints a step description with a visual separator.

    Parameters:
        step_description : str
            A description of the current step in the process.

    Returns:
        None

    The description is prefixed with an arrow for clear visibility.
    """
    print(f"-> {step_description}")
    print()


# Function to print the footer display
def print_footer():
    """
    Prints a footer display indicating the successful creation of Xomics files.

    Returns:
        None

    Displays a message indicating that files have been created and are located in the 'XomicsFiles' folder.
    """
    print("="*60)
    print("Xomics files have been successfully created. You can find them in the 'XomicsFiles' folder.")
    print("="*60)


