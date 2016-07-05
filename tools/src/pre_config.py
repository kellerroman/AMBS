#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################
import struct
import sys

from configparser import SafeConfigParser

parser = SafeConfigParser()
parser.read('config.cfg')

integer_init = -1
real_init = 0.0

integers = []
reals = []
strings = []

integers.append(["control","iterations",integer_init])
integers.append(["output","solution",integer_init])
integers.append(["output","residual",integer_init])

integers.append(["model","equation",integer_init])
integers.append(["model","turbulence",integer_init])

integers.append(["disc","space_order",integer_init])
integers.append(["disc","riemann_solver",integer_init])
integers.append(["disc","timestep_method",integer_init])
integers.append(["disc","time_order",integer_init])


reals.append(["disc","CFL",real_init])
reals.append(["disc","timestep",real_init])

int_in_file=[]
real_in_file=[]
string_in_file=[]


for var in integers:
    int_in_file.append(0)

for var in reals:
    real_in_file.append(0)

for var in strings:
    string_in_file.append(0)


version = 10000*len(integers)+100*len(reals)+len(strings)
version = 1

for section_name in parser.sections():
    #print ('Section:', section_name)
#     print ('  Options:', parser.options(section_name))
    for name, value in parser.items(section_name):
        #print ('  %s = %s' % (name, value))
        name = name.lower()
        value = value.lower()
        section_name = section_name.lower()
        for i,var in enumerate(integers):
            if (var[0].lower() == section_name):
                if (var[1].lower() == name):
                    var[2] = value
                    int_in_file[i] = 1

        for i,var in enumerate(reals):
            if (var[0].lower() == section_name):
                if (var[1].lower() == name):
                    var[2] = value
                    real_in_file[i] = 1

        for i,var in enumerate(strings):
            if (var[0].lower() == section_name):
                if (var[1].lower() == name):
                    var[2] = value
                    string_in_file[i] = 1
    
for i,var in enumerate(int_in_file):
    if (var == 0):
        print (integers[i][1],"not found in file")
        exit (1)
#
#
#   KOMMANDOZEILEN PARAMETER
#
#
for i in sys.argv[1:]:
    section_name = i.split("/")[0].strip().lower()
    name = i.split("/")[1].split("=")[0].strip().lower()
    value = i.split("/")[1].split("=")[1].strip().lower()
    for i,var in enumerate(integers):
        if (var[0].lower() == section_name):
            if (var[1].lower() == name):
                var[2] = value
                int_in_file[i] = 1

    for i,var in enumerate(reals):
        if (var[0].lower() == section_name):
            if (var[1].lower() == name):
                var[2] = value
                real_in_file[i] = 1

    for i,var in enumerate(strings):
        if (var[0].lower() == section_name):
            if (var[1].lower() == name):
                var[2] = value
                string_in_file[i] = 1
# Write binary data to a file
with open('config.bin', 'wb') as f:
    f.write(version.to_bytes(4, byteorder='little', signed=True))
    for var in integers:
        f.write(int(var[2]).to_bytes(4, byteorder='little', signed=True))
    for var in reals:
        f.write(struct.pack('d',float(var[2])))
