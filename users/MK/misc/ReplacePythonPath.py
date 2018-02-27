# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:46:26 2017

@author: MichaelEK

Script to change the paths to packages/modules from the pythonpath and create __init__.py files if needed.
"""

import os, fnmatch

############################
## Get pythonpath and file path

#py_path = os.environ['PYTHONPATH'].split(os.pathsep)[0]
py_path = r'D:\Executables\PythonScripts'

print('Python path is ' + py_path)

path1 = os.path.dirname(os.path.abspath(__file__))
path2 = os.path.join(os.path.split(path1)[0] , 'core')

############################
## Search for all files recursively

append_path1 = os.path.relpath(path2, py_path)
append_path2 = append_path1.replace('\\', '.') + '.'

path_base = os.path.split(path1)[0]

print('Project path is ' + path_base)

for root, dirs, files in os.walk(path_base):

    ## Create an init file if it doesn't exist
    init1 = os.path.join(root, '__init__.py')
    try:
        fh = open(init1,'r')
    except:
        fh = open(init1,'w')

    ## Iterate through each file and change the 'from users.MK.misc.test_project.core.' paths
    for filename in fnmatch.filter(files, '*.py'):
        # Read in the file
        filename1 = os.path.join(root, filename)
        with open(filename1, 'r') as file :
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace('from users.MK.misc.test_project.core.', 'from ' +  append_path2)

        # Write the file out again
        with open(filename1, 'w') as file:
            file.write(filedata)
