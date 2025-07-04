#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------------------------
# This script is used to get the data(enthalpy and pressure) 
# of phase transition calculation.
#--------------------------------------------------------------------------------
import os
import re
import pandas as pd

from convergence import OutcarRead

def enthalpy_pressure_extract():
    """Get the enthalpy and pressure data from all relax.
    """
    data_list = []
    file_list = []
    files = os.listdir()
    for file in files:
        # Match the file name like Number add GPa.
        if re.search(r"\d+(\.\d+)?\s*GPa", file):
            outcomes = []
            file_list.append(file)
            outcar = OutcarRead(file + "/OUTCAR")
            pressure = outcar.get_pressure()
            enthaply = outcar.get_free_energy()
            outcomes.extend((pressure, enthaply[-1]))
            data_list.append(outcomes)
    # Data output, fomart: .csv
    data_list.sort(key=lambda x:x[0])
    name = ["pressure", "enthalpy"]
    data = pd.DataFrame(columns=name, data=data_list)
    data.to_csv("energy.csv", index=False)
    print("energy.scv is generated")
        
def main():
    enthalpy_pressure_extract()      
        
if __name__ == "__main__":
    main()