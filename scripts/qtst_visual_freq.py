#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-----------------------------------------------
# The script is used to visulaize the imaginary frequency
# aiming to manully adjust atoms' site and find the stable site.
#-----------------------------------------------
import re
import sys
import subprocess
import numpy as np

class CmdRrror(Exception):
    def __init__(self, errorinfo):
        super().__init__(self)  # 初始化父类
        self.errorinfo = errorinfo

    def __str__(self):
        return 'Command Execution Error: ' + self.errorinfo

def execCmd(command):
    pipe = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (content, error) = (pipe.stdout.readlines(), pipe.stderr.read().decode())
    if content:
        for i in range(len(content)):
            content[i] = content[i].decode()
        return content

    if error != "":
        raise CmdRrror(error)

def writeXYZ(file_name, num_structures, elements, num_atoms, coordinates):
    """
    :param num_structures: number of structures in xyz file, int
    :param elements: [element1, element2, ...], str
    :param num_atoms: [number_of_element1, number_of_element2, ...], int
    :param coordinates: [[x1, y1, z1], [x2, y2, z2], ...], float
    :return:
    """
    total_atoms = sum(num_atoms)
    index = 0
    with open(file_name, 'w') as output_file:
        while num_structures != 0:
            output_file.write('%d\n' % total_atoms)
            output_file.write('create form python\n')
            for element_id, element in enumerate(elements):
                for element_index in range(num_atoms[element_id]):
                    output_file.write("%2s    %16.10f    %16.10f    %16.10f\n" %
                                      (element, coordinates[index][0], coordinates[index][1],
                                       coordinates[index][2]))
                    index += 1
            num_structures -= 1

if len(sys.argv) < 3:
    print('')
    print('Usage: %s freq1 freq2 ... frames scale' % sys.argv[0].split('/')[-1])
    print('')
    exit(1)

space = re.compile(r'\s+')

def readFreq(file_name):
    """
    :param file_name: str
    :return: [coordinates, dcoordinates]
    coordinates: 原子坐标(x,y,z)
    dcoordinates: 震动坐标(dx,dy,dz)
    """
    coordinates = []
    dcoordinates = []
    with open(file_name, 'r') as freq_file:
        content = freq_file.readlines()
    index = 2
    while index < len(content):
        line = content[index].strip()
        if line == '':
            break

        line = list(map(float, space.split(line)))
        coordinates.append(np.array([line[0], line[1], line[2]]))
        dcoordinates.append(np.array([line[3], line[4], line[5]]))
        index += 1

    return np.stack(coordinates), np.stack(dcoordinates)


def get_elements_info():
    # 元素符号
    cmd = "grep 'VRHFIN' OUTCAR"
    try:
        content = execCmd(cmd)
    except CmdRrror:
        print(CmdRrror)
        print("")
        exit(1)
    pattern = re.compile(r'=(.*?):')
    elements = []
    for line in content:
        line = line.strip()
        elements.append(pattern.search(line).group(1))

    # 各元素原子数量
    cmd = "grep 'ions per type' OUTCAR"
    try:
        atom_num = execCmd(cmd)
    except CmdRrror:
        print(CmdRrror)
        print("")
        exit(1)

    num_atoms = list(map(int, space.split(atom_num[0].strip().split('=')[-1].strip())))
    return elements, num_atoms


print('')
print('################ This script makes animation of vibration ################')
print('')
elements, num_atoms = get_elements_info()
frames = int(sys.argv[-2])
dframe = 1 / float(frames)
scale = float(sys.argv[-1])
for freq_file in sys.argv[1:-2]:
    print('                            Processing %s' % freq_file)
    freq_pos, vibration = readFreq(freq_file)

    num_structures = 0
    coordinates = []

    # 0→1
    for i in range(0, frames + 1):
        num_structures += 1
        coordinates.append(freq_pos + vibration * dframe * i * scale)

    # 1→0
    for i in range(frames, -1, -1):
        num_structures += 1
        coordinates.append(freq_pos + vibration * dframe * i * scale)

    # 0→-1
    for i in range(-1, -frames - 1, -1):
        num_structures += 1
        coordinates.append(freq_pos + vibration * dframe * i * scale)

    # -1→0
    for i in range(-frames, 0):
        num_structures += 1
        coordinates.append(freq_pos + vibration * dframe * i * scale)

    coordinates = np.concatenate(coordinates)

    writeXYZ('%s.xyz' % freq_file, num_structures, elements, num_atoms, coordinates)

print('')
print('                ------------------ Done ------------------\n')

def main():
    pass

if __name__ == "__main__":
    main()