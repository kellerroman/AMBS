#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################
import struct

from configparser import SafeConfigParser

parser = SafeConfigParser()
parser.read('config.cfg')
version = 2
iter        = 1
solution    = 1
residual    = 1
equation    = 1
dimension   = 1
space_disc  = 1
space_order = 1
dt_method    = 2
riemann_solver = 2
CFL         = 0.1
timestep    = 1E-4
boundary_cell_output = 0

for section_name in parser.sections():
    #print ('Section:', section_name)
#     print ('  Options:', parser.options(section_name))
    for name, value in parser.items(section_name):
        #print ('  %s = %s' % (name, value))
        name = name.lower()
        value = value.lower()
        section_name = section_name.lower()
        if (section_name == "control" and name == "iterations"):
            iter = int(value)
        elif (section_name == "output" and name == "solution"):
            solution = int(value)
        elif (section_name == "output" and name == "residual"):
            residual = int(value)
        elif (section_name == "output" and name == "boundary cell"):
            if (value == "false"):
                boundary_cell_output = 0
            elif (value == "true"):
                boundary_cell_output = 1
            else:
                boundary_cell_output = int(value)
        elif (section_name == "disc" and name == "timestep"):
            timestep = float(value)
        elif (section_name == "disc" and name == "cfl"):
            CFL = float(value)
        elif (section_name == "disc" and name == "timestep method"):
            if (value == "dt"):
                dt_method = 1
            elif (value == "cfl"):
                dt_method = 2
            elif (value == "mindt"):
                dt_method = 3
            else:
                print ("Timestep Method nicht erkannt:",value)
                stop
        elif (section_name == "disc" and name == "riemann solver"):
            if (value == "rhll"):
                riemann_solver = 1
            elif (value == "roe"):
                riemann_solver = 2
#             elif (value == "ausm"):
#                 riemann_solver = 3
#             elif (value == "ausmp"):
#                 riemann_solver = 4
#             elif (value == "ausmpw"):
#                 riemann_solver = 5
#             elif (value == "ausmpwp"):
#                 riemann_solver =6
            else:
                print ("Riemann Solver nicht erkannt:",value)
                stop 
        elif (section_name == "disc" and name == "space disc"):
            if (value == "notvd"):
                space_disc = 0
            elif (value == "muscl"):
                space_disc = 1
            else:
                print ("Spatial Discretization mechanism nicht erkannt:",value)
                stop 
        
        elif (section_name == "disc" and name == "space order"):
            space_order = int(value)  
        else:
            print (section_name,":",name,"=",value)
    
# Write binary data to a file
with open('config.bin', 'wb') as f:
    f.write(version                      .to_bytes(4, byteorder='little', signed=True))
    f.write(iter                             .to_bytes(4, byteorder='little', signed=True))
    f.write(solution                    .to_bytes(4, byteorder='little', signed=True))
    f.write(residual                    .to_bytes(4, byteorder='little', signed=True))
    f.write(equation                   .to_bytes(4, byteorder='little', signed=True))
    f.write(dimension                 .to_bytes(4, byteorder='little', signed=True))
    f.write(space_disc                .to_bytes(4, byteorder='little', signed=True))
    f.write(space_order             .to_bytes(4, byteorder='little', signed=True))
    f.write(dt_method                .to_bytes(4, byteorder='little', signed=True))
    f.write(riemann_solver       .to_bytes(4, byteorder='little', signed=True))
    f.write(struct.pack('d',CFL))
    f.write(struct.pack('d',timestep))
    f.write(boundary_cell_output       .to_bytes(4, byteorder='little', signed=True))
