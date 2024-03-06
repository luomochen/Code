#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------------------------
# This script is used to get the data(enthalpy and pressure) 
# of phase transition calculation.
#--------------------------------------------------------------------------------

import os
import re
import pandas as pd

from Vasp_Read import OutcomeFileRead

files = os.listdir()
data_list = []
file_list = []
converinfo_list = []
# Get the convergence information of all relax.
for file in files:
    if re.search(r"\D_\d+GPa", file):
        file_list.append(file)
        outcome = OutcomeFileRead(file + "/OUTCAR")
        converinfo_list.append(outcome.get_convergence())
# If there have non-convergence calculation, 
# The script will not generate a result file.
if False in converinfo_list:
    for i in range(len(file_list)):
        if converinfo_list[i] == False:
            print(f"{file_list[i]} did not converge!")
else:
    for file in files:
        if re.search(r"\D+_\d+GPa", file):
            outcomes = []
            outcome = OutcomeFileRead(file + "/OUTCAR")
            pressure = outcome.get_pressure()
            pp_type = outcome.get_pp_type()
            enthaply = outcome.get_free_energy()
            outcomes.extend((pressure, pp_type, enthaply[-1]))
            data_list.append(outcomes)
    data_list.sort(key=lambda x:x[0])
    name = ["pressure", "pp", "enthalpy"]
    data = pd.DataFrame(columns=name, data=data_list)
    data.to_csv("outcome.csv", index=False)