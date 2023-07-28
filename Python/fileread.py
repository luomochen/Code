import linecache as lc
filename = 'CONTCAR.txt'
with open(filename) as file_object:
    coordinate = lc.getline(filename, 9)
    print(coordinate)