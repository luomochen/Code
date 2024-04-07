#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------------------------------------------------
# This module is used to read/write vasp file and get output information.
#------------------------------------------------------------------------
import re
from subprocess import getstatusoutput

class StructureFileRead:

    """Read vasp structure file.
    
    This class is used to read vasp structure file like POSCAR, CONTCAR and *.vasp.
    Shell command is used to get structure lattice vector, direct coordinate and so on.
    Warning!!! if you want a higher accuracy of the data, please use string to deliver
    your data.
    
    Attributes:
        filename: a long string indicating the vasp structure file path.  
        float_type: a boolean value that control the delivering data type.
    """ 
    
    def __init__(self, filename, float_type) -> None:
        """Inits class"""
        self.filename = filename
        self.float_type = float_type
        
    def string2float(self, string_2dlist):
        """Transform two dimension string list to float list

        Args:
            string_2dlist (list[list[str]])

        Returns:
            list[list[str|float]]
        """
        if self.float_type == True:
            try:
                string_2dlist = [[float(string) for string in string_list] 
                                 for string_list in string_2dlist]
            except ValueError:
                string_2dlist = [['NaN', 'NaN', 'NaN']]
            return string_2dlist
        else:
            return string_2dlist

    def lattice_vector_matrix(self):
        """Get the lattice vector matrix.

        Returns:
            list[list[float|str]]: The lattice vector matrix which has two kind of type.
        """
        lattice_vector_matrix = []
        for i in range(3):
            input = "sed -n '" + str(i+3) + "p' "+ self.filename
            get_vec = getstatusoutput(input)
            if get_vec[0] == 0:
                lattice_vector = get_vec[1].split()
                lattice_vector_matrix.append(lattice_vector)
            else:
                lattice_vector_matrix.append(['NaN', 'NaN', 'NaN'])
                break
            lattice_vector_matrix = self.string2float(lattice_vector_matrix)
            return lattice_vector_matrix 

    # 获取原子种类。
    def atom_species(self):
        input = "sed -n '6p' " + self.filename
        get_atom_species = getstatusoutput(input)
        if get_atom_species[0] == 0:
            atom_species = get_atom_species[1].split()
            return(atom_species)
        else:
            return ['NaN']
    
    # 获取原子个数。
    def atom_numbers(self):
        input = "sed -n '7p' " + self.filename
        get_atom_numbers = getstatusoutput(input)
        if get_atom_numbers[0] == 0:
            atom_numbers = get_atom_numbers[1].split()
            atom_numbers = [int(number) for number in atom_numbers]
            return(atom_numbers)
        else:
            return ['NaN']

    # 获取原子的分数坐标。
    def atom_direct(self):
        atom_numbers = self.atom_numbers()
        whole_number = 0
        for atom_number in atom_numbers:
            whole_number = whole_number + atom_number
        input = "grep Selective " + self.filename
        get_select = getstatusoutput(input)
        if get_select[0] == 1:
            start_line_number = 9
        else:
            start_line_number = 10
        atom_directs = []
        for i in range(start_line_number, whole_number + start_line_number):
            direct = []
            input = "sed -n '" + str(i) + "p' "+ self.filename
            get_direct = getstatusoutput(input)
            if get_direct[0] == 0:
                direct = get_direct[1].split()
            else:
                direct = ['NaN', 'NaN', 'NaN']
                break
            atom_directs.append(direct)
        atom_directs = self.string2float(atom_directs)
        return atom_directs
    
    # 获取指定原子坐标。
    def assign_atom_direct(self, assign_species):
        atom_species = self.atom_species()
        atom_numbers = self.atom_numbers()
        atom_directs = self.atom_direct()
        if assign_species in atom_species:
            assign_index = atom_species.index(assign_species)
            final_index = 0
            for number in atom_numbers[:assign_index+1]:
                final_index = final_index + number
            start_index = final_index - atom_numbers[assign_index]
            if assign_index == 0:
                assign_atom_direct = self.string2float(atom_directs[:final_index])
            else:
                assign_atom_direct = self.string2float(atom_directs[start_index:final_index])
            return assign_atom_direct
        else:
            return ['NaN', 'NaN', 'NaN']

class OutcomeFileRead:

    def __init__(self, filename) -> None:
        self.filename = filename
    
    def get_pp_type(self):
        input = "gawk '/LEXCH/' " + self.filename
        get_pp_type = getstatusoutput(input)
        if get_pp_type[0] == 0:
            if re.search(r"CA", get_pp_type[1]):
                return ["LDA"]
            if re.search(r"PE", get_pp_type[1]): 
                return ["GGA"]
            else:
                return ['NaN']

    def get_convergence(self):
        input = "gawk '/required accuracy/' " + self.filename
        get_convergence = getstatusoutput(input)
        if get_convergence[0] == 0:
            if re.search(r"required", get_convergence[1]):
                return True
            else:
                return False
        else:
            return ['NaN']
        
    def get_free_energy(self):
        input = "gawk '/TOTEN/' " + self.filename
        get_free_energy = getstatusoutput(input)
        if get_free_energy[0] == 0:
            free_energy_list = re.findall(r"=(.+?)eV", get_free_energy[1])
            try:
                free_energy_list = [float(free_energy.strip()) 
                                    for free_energy in free_energy_list]
            except ValueError:
                free_energy_list = ['NaN']
            return free_energy_list
        else:
            return ['NaN']

    def get_internal_energy(self):
        input = "gawk '/without/' " + self.filename
        get_internal_energy = getstatusoutput(input)
        if get_internal_energy[0] == 0:
            internal_energy_list = re.findall(r"=(.+?)energy", get_internal_energy[1])
            try:
                internal_energy_list = [float(internal_energy.strip()) 
                                        for internal_energy in internal_energy_list]
            except ValueError:
                internal_energy_list = ['NaN']
            return internal_energy_list
        else:
            return ['NaN']
        
    def get_enthalpy(self):
        input = "gawk '/enthalpy/' " + self.filename
        get_enthalpy = getstatusoutput(input)
        if get_enthalpy[0] == 0:
            enthalpy_list = re.findall(r"=(.+?)eV", get_enthalpy[1])
            try:
                enthalpy_list = [float(enthalpy.strip()) for enthalpy in enthalpy_list]
            except ValueError:
                enthalpy_list = ['NaN']
            return enthalpy_list
        else:
            return ['NaN']
        
    def get_pressure(self):
        input = "gawk '/PSTRESS/' " + self.filename
        get_pressure = getstatusoutput(input)
        if get_pressure[0] == 0:
            pressure = re.findall(r"=(.+?)pullay", get_pressure[1])
            try:
                pressure = float(pressure[0].strip()) / 10
            except ValueError:
                pressure = ['NaN']
            return pressure
        else:
            pressure = ['NaN']
            
    def get_frequency(self):
        input = "gawk '/THz/' " + self.filename
        get_frequency = getstatusoutput(input)
        if get_frequency[0] == 0:
            frequency_list = re.findall(r"=(.+?)THz", get_frequency[1])
            try:
                frequency_list = [float(frequency.strip()) for frequency in frequency_list]
            except ValueError:
                frequency_list = ['NaN']
            return frequency_list
        else:
            return ['NaN']