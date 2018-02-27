# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:46:26 2017

@author: MichaelEK

Script to change the paths to packages/modules from the pythonpath and create __init__.py files if needed.
"""

import os, fnmatch, re

############################
## Get pythonpath and file path

#py_path = os.environ['PYTHONPATH'].split(os.pathsep)[0]
py_path = r'D:\Executables\PythonScripts'

print('Python path is ' + py_path)

print('Absolute path is ' + os.path.abspath(__file__))

path1 = os.path.dirname(os.path.abspath(__file__))
path2 = os.path.join(path1, 'core')

print('Path 1 is ' + path1)
print('Path 2 is ' + path2)

############################
## Search for all files recursively

append_path1 = os.path.relpath(path2, py_path)
append_path2 = append_path1.replace('\\', '.') + '.'
append_path3 = append_path1.replace('\\', '.') + ' '

print('Append path 1 is ' + append_path1)
print('Append path 2 is ' + append_path2)

path_base = path1

print('Project path is ' + path_base)

for root, dirs, files in os.walk(path_base):

    ## Create an init file if it doesn't exist
    init1 = os.path.join(root, '__init__.py')
    try:
        fh = open(init1,'r')
    except:
        fh = open(init1,'w')

    print('Root is ' + root)

    ## Iterate through each file and change the 'from core.' paths
    for filename in fnmatch.filter(files, '*.py'):
        # Read in the file
        filename1 = os.path.join(root, filename)
        with open(filename1, 'r') as file :
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace('from core.', 'from ' +  append_path2)
        filedata = filedata.replace('from core ', 'from ' +  append_path3)

        # Write the file out again
        with open(filename1, 'w') as file:
            file.write(filedata)
