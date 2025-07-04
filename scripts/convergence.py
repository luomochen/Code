import re
from subprocess import getstatusoutput

class OutcarRead:

    def __init__(self, filename) -> None:
        self.filename = filename

    def get_convergence(self):
        """Determine the convergence of calculation.
        """
        converg_info_req = "gawk '/required accuracy/' " + self.filename
        converg_info = getstatusoutput(converg_info_req)
        if converg_info[0] == 0:
            if re.search(r"required", converg_info[1]):
                return True
            else:
                return False
        else:
            return 'NaN'

    def get_free_energy(self):
        if self.get_convergence() == True:
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
        elif self.get_convergence() == False:
            return ['UnConv']
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
                pressure = 'NaN'
            return pressure
        else:
            pressure = 'NaN'

def main():
    energy = OutcarRead("OUTCAR").get_free_energy()
    print(energy)

if __name__ == "__main__":
    main()