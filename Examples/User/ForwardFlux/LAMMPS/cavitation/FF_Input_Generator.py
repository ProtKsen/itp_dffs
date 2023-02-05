#! /usr/bin/env python

# This file is used to take a template input json file for forward flux
# and create a input file for ssages that uses multiple walkers.

import json 
import numpy as np
import os
from random import randint

class LessPrecise(float):
    def __repr__(self):
        return str(self)

if not os.path.exists('ini_files'):
        os.makedirs('ini_files');

if not os.path.exists('profiles'):
        os.makedirs('profiles');

if not os.path.exists('dumpsr'):
        os.makedirs('dumpsr');

if not os.path.exists('logs'):
        os.makedirs('logs');

if not os.path.exists('ptu'):
        os.makedirs('ptu');



def lammps_random_number(lammps_filename, nWalkers):
    """Takes a lammps input file, finds the seed number of langevin, replaces it with a random number, and
	generates new input files
	"""
    randon_numbers = set()
    while (len(randon_numbers) < nWalkers):
        randon_numbers.add(randint(0, 100000))
    randon_numbers = list(randon_numbers)

    in_lammps_file = open(lammps_filename, 'r', encoding='cp1252')

    for i in range(0, nWalkers):
        out_lammps_file = open('ini_files/'+lammps_filename + '-' + str(i), 'w', encoding='cp1252')
        for line in in_lammps_file:
            if line.find('seed_date equal') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 3:
                        line = line.replace(columns[3], str(randon_numbers[i]))
                        line = line + '\n'
            if line.find('log') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 1:
                        line = line.replace(columns[1], 'logs/log'+str(i))
                        line = line + '\n'
            if line.find('all ave/time') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 10:
                        line = line.replace(columns[11], 'ptu/pt_'+str(i)+'.txt')
                        line = line + '\n'
            if line.find('read_restart') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 1:
                        line = line.replace(columns[1], 'restarts/rest_'+str(i)+'_.200000')  
                        line = line + '\n'
            if line.find('read_dump') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 1:
                        line = line.replace(columns[1], 'dumps/dump_end'+str(i)+'.txt')  
                        line = line + '\n'
            if line.find('read_restart') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 1:
                        line = line.replace(columns[1], 'rests_finish/rest_'+str(i))  
                        line = line + '\n'
            if line.find('dump gcg') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 3:
                        line = line.replace(columns[5], 'dumpsr/dump_'+str(i)+'.txt')
                        line = line + '\n'
            if line.find('fix aveRUN') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 3:
                        line = line.replace(columns[12], 'profiles/Profiles_run_'+str(i)+'.txt')
                        line = line + '\n'
            out_lammps_file.write(line)
        in_lammps_file.seek(0)
        out_lammps_file.close()

    in_lammps_file.close()


# User must set these variables
nWalkers = 4
input_filename = "cr_nvt"
nsurf = 51
interfaces = np.linspace(1, 101, nsurf ) 
trials = np.empty(nsurf , dtype=int)
trials.fill(10000)
# Use if you have many equally-spaced interfaces

#Generate the new lammps input files
lammps_random_number(input_filename, nWalkers)

# Open template and load in the json data. 
root = {}
with open('input.json') as f:
    root = json.load(f)


root["walkers"] = nWalkers
root['methods'][0]['interfaces'] = interfaces.tolist()
root['methods'][0]['trials'] = trials.tolist()

input = []
for i in range(0, nWalkers):
    input.append('ini_files/'+input_filename+ '-' + str(i))
root['input'] = input
	
# Convert python dictionary into JSON file
with open('Input_end.json', 'w') as f:
        json.dump(root, f, indent=4, separators=(',', ': '))
