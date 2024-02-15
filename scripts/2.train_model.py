## Objective: This script uses the xgboost package to predict NBSS parameters as a function of environmental variables

## Functions: All functions required to run this step are included in funcs_predict_NBSS.py

## Modules:
from natsort import natsorted
import pandas as pd
import numpy as np # Arrays
import datetime
from pathlib import Path # Handling of path object
# Main functions to download environmental datasets and predict NBSS:
try:
    from funcs_predict_NBSS import *
except:
    from scripts.funcs_predict_NBSS import *

## Workflow starts here:
path_to_datafile=path_to_git / 'data' /'Model'
# Load dataframe generated on step 0: NBSS parameters per bin per taxon
data=pd.read_csv(path_to_datafile /'Model_input.csv',dtype={'year':int,'month':int})