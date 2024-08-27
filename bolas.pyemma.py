import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyemma
from pyemma.util.contexts import settings
import mdtraj as md
import glob
import pandas as pd
from itertools import product
import pickle
import copy
import seaborn as sns
from pandas import DataFrame
import scipy
import operator
from scipy.optimize import minimize
import math
from pymbar import MBAR
from pymbar import FES
from mpl_toolkits.mplot3d import axes3d

#=====================================Load Trajs===============================================

pdb = md.load("./TPS_22500.cp/gn_005_xpr2.noWAT.nc.start.pdb")


files_TPS_2500_cp = glob.glob("./TPS_2500.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_5000_cp = glob.glob("./TPS_5000.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_7500 = glob.glob("./TPS_7500.BOLAS-2/*.noWAT.align.nc")
files_TPS_10000 = glob.glob("./TPS_10000.BOLAS-2/*.noWAT.align.nc")
files_TPS_12500 = glob.glob("./TPS_12500.BOLAS-2/*.noWAT.align.nc")
files_TPS_15000_cp = glob.glob("./TPS_15000.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_17500_cp = glob.glob("./TPS_17500.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_20000_cp = glob.glob("./TPS_20000.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_22500_cp = glob.glob("./TPS_22500.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_25000 = glob.glob("./TPS_25000.BOLAS-2/*.noWAT.align.nc")
files_TPS_27500_cp = glob.glob("./TPS_27500.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_30000 = glob.glob("./TPS_30000.BOLAS-2/*.noWAT.align.nc")
files_TPS_32500 = glob.glob("./TPS_32500.BOLAS-2/*.noWAT.align.nc")
files_TPS_35000_cp = glob.glob("./TPS_35000.cp.BOLAS-2/*.noWAT.align.nc")
files_TPS_37500_cp = glob.glob("./TPS_37500.cp.BOLAS-2/*.noWAT.align.nc")

files_TPS_2500_cp.sort() 
files_TPS_5000_cp.sort() 
files_TPS_7500.sort()
files_TPS_10000.sort()
files_TPS_12500.sort()
files_TPS_15000_cp.sort() 
files_TPS_17500_cp.sort()
files_TPS_20000_cp.sort()
files_TPS_22500_cp.sort() 
files_TPS_25000.sort()
files_TPS_27500_cp.sort()
files_TPS_30000.sort()
files_TPS_35000_cp.sort()
files_TPS_37500_cp.sort()

files = files_TPS_2500_cp + files_TPS_5000_cp + files_TPS_7500 + files_TPS_10000 + files_TPS_12500 + files_TPS_15000_cp + files_TPS_17500_cp + files_TPS_20000_cp + files_TPS_22500_cp + files_TPS_25000 + files_TPS_27500_cp + files_TPS_30000 + files_TPS_35000_cp + files_TPS_37500_cp 


dist_feat=pyemma.coordinates.featurizer(pdb)                            
dist_feat.add_distances(np.array([[105,2416],[66,1956]]))   # 3P [66] - 62P [1956]; 4O4' [105] - 76N1 [2416]             
dist_data=pyemma.coordinates.load(files, features=dist_feat)
dist_data_concatenated = np.concatenate(dist_data)




dist_data_2500=pyemma.coordinates.load(files_TPS_2500_cp, features=dist_feat)
dist_data_2500_concatenated = np.concatenate(dist_data_2500)

dist_data_5000=pyemma.coordinates.load(files_TPS_5000_cp, features=dist_feat)
dist_data_5000_concatenated = np.concatenate(dist_data_5000)

dist_data_7500=pyemma.coordinates.load(files_TPS_7500, features=dist_feat)
dist_data_7500_concatenated = np.concatenate(dist_data_7500)

dist_data_10000=pyemma.coordinates.load(files_TPS_10000, features=dist_feat)
dist_data_10000_concatenated = np.concatenate(dist_data_10000)

dist_data_12500=pyemma.coordinates.load(files_TPS_12500, features=dist_feat)
dist_data_12500_concatenated = np.concatenate(dist_data_12500)

dist_data_15000=pyemma.coordinates.load(files_TPS_15000_cp, features=dist_feat)
dist_data_15000_concatenated = np.concatenate(dist_data_15000)

dist_data_17500=pyemma.coordinates.load(files_TPS_17500_cp, features=dist_feat)
dist_data_17500_concatenated = np.concatenate(dist_data_17500)

dist_data_20000=pyemma.coordinates.load(files_TPS_20000_cp, features=dist_feat)
dist_data_20000_concatenated = np.concatenate(dist_data_20000)

dist_data_22500=pyemma.coordinates.load(files_TPS_22500_cp, features=dist_feat)
dist_data_22500_concatenated = np.concatenate(dist_data_22500)

dist_data_25000=pyemma.coordinates.load(files_TPS_25000, features=dist_feat)
dist_data_25000_concatenated = np.concatenate(dist_data_25000)

dist_data_27500=pyemma.coordinates.load(files_TPS_27500_cp, features=dist_feat)
dist_data_27500_concatenated = np.concatenate(dist_data_27500)

dist_data_30000=pyemma.coordinates.load(files_TPS_30000, features=dist_feat)
dist_data_30000_concatenated = np.concatenate(dist_data_30000)

dist_data_35000=pyemma.coordinates.load(files_TPS_35000_cp, features=dist_feat)
dist_data_35000_concatenated = np.concatenate(dist_data_35000)

dist_data_37500=pyemma.coordinates.load(files_TPS_37500_cp, features=dist_feat)
dist_data_37500_concatenated = np.concatenate(dist_data_37500)


def read_file(filename):
   """Read contents of the specified file.

   Parameters:
   -----------
   filename : str
      The name of the file to be read

   Returns:
   lines : list of str
      The contents of the file, split by line

   """

   infile = open(filename, 'r')
   lines = infile.readlines()
   infile.close()

   return lines
   

race_accept_TPS_2500_cp = read_file("./TPS_2500.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_5000_cp = read_file("./TPS_5000.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_7500 = read_file("./TPS_7500.BOLAS-2/race_accepted.log")
race_accept_TPS_10000 = read_file("./TPS_10000.BOLAS-2/race_accepted.log")
race_accept_TPS_12500 = read_file("./TPS_12500.BOLAS-2/race_accepted.log")
race_accept_TPS_15000_cp = read_file("./TPS_15000.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_17500_cp = read_file("./TPS_17500.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_20000_cp = read_file("./TPS_20000.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_22500_cp = read_file("./TPS_22500.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_25000 = read_file("./TPS_25000.BOLAS-2/race_accepted.log")
race_accept_TPS_27500_cp = read_file("./TPS_27500.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_30000 = read_file("./TPS_30000.BOLAS-2/race_accepted.log")
race_accept_TPS_32500 = read_file("./TPS_32500.BOLAS-2/race_accepted.log")
race_accept_TPS_35000_cp = read_file("./TPS_35000.cp.BOLAS-2/race_accepted.log")
race_accept_TPS_37500_cp = read_file("./TPS_37500.cp.BOLAS-2/race_accepted.log")


accepted_ind_TPS_2500_cp=[]
for line in race_accept_TPS_2500_cp:
    accepted_ind_TPS_2500_cp.append(line.split()[5])
pot_TPS_2500_cp=[]
for item in accepted_ind_TPS_2500_cp:
    if int(item)<10:
        f1=open("./TPS_2500.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_2500_cp.append(pot)
        f1.close()
        f2=open("./TPS_2500.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_2500_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_2500.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_2500_cp.append(pot)
        f1.close()
        f2=open("./TPS_2500.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_2500_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_2500.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_2500_cp.append(pot)
        f1.close()
        f2=open("./TPS_2500.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_2500_cp.append(pot)
        f2.close()
        
        


accepted_ind_TPS_5000_cp=[]
for line in race_accept_TPS_5000_cp:
    accepted_ind_TPS_5000_cp.append(line.split()[5])
pot_TPS_5000_cp=[]
for item in accepted_ind_TPS_5000_cp:
    if int(item)<10:
        f1=open("./TPS_5000.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_5000_cp.append(pot)
        f1.close()
        f2=open("./TPS_5000.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_5000_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_5000.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_5000_cp.append(pot)
        f1.close()
        f2=open("./TPS_5000.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_5000_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_5000.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_5000_cp.append(pot)
        f1.close()
        f2=open("./TPS_5000.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_5000_cp.append(pot)
        f2.close()
        



accepted_ind_TPS_7500=[]
for line in race_accept_TPS_7500:
    accepted_ind_TPS_7500.append(line.split()[5])
pot_TPS_7500=[]
for item in accepted_ind_TPS_7500:
    if int(item)<10:
        f1=open("./TPS_7500.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_7500.append(pot)
        f1.close()
        f2=open("./TPS_7500.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_7500.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_7500.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_7500.append(pot)
        f1.close()
        f2=open("./TPS_7500.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_7500.append(pot)
        f2.close()
    else:
        f1=open("./TPS_7500.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_7500.append(pot)
        f1.close()
        f2=open("./TPS_7500.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_7500.append(pot)
        f2.close()
        



accepted_ind_TPS_10000=[]
for line in race_accept_TPS_10000:
    accepted_ind_TPS_10000.append(line.split()[5])
pot_TPS_10000=[]
for item in accepted_ind_TPS_10000:
    if int(item)<10:
        f1=open("./TPS_10000.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_10000.append(pot)
        f1.close()
        f2=open("./TPS_10000.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_10000.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_10000.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_10000.append(pot)
        f1.close()
        f2=open("./TPS_10000.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_10000.append(pot)
        f2.close()
    else:
        f1=open("./TPS_10000.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_10000.append(pot)
        f1.close()
        f2=open("./TPS_10000.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_10000.append(pot)
        f2.close()
        
        


accepted_ind_TPS_12500=[]
for line in race_accept_TPS_12500:
    accepted_ind_TPS_12500.append(line.split()[5])
pot_TPS_12500=[]
for item in accepted_ind_TPS_12500:
    if int(item)<10:
        f1=open("./TPS_12500.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_12500.append(pot)
        f1.close()
        f2=open("./TPS_12500.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_12500.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_12500.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_12500.append(pot)
        f1.close()
        f2=open("./TPS_12500.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_12500.append(pot)
        f2.close()
    else:
        f1=open("./TPS_12500.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_12500.append(pot)
        f1.close()
        f2=open("./TPS_12500.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_12500.append(pot)
        f2.close()        
        

        
accepted_ind_TPS_15000_cp=[]
for line in race_accept_TPS_15000_cp:
    accepted_ind_TPS_15000_cp.append(line.split()[5])
pot_TPS_15000_cp=[]
for item in accepted_ind_TPS_15000_cp:
    if int(item)<10:
        f1=open("./TPS_15000.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_15000_cp.append(pot)
        f1.close()
        f2=open("./TPS_15000.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_15000_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_15000.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_15000_cp.append(pot)
        f1.close()
        f2=open("./TPS_15000.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_15000_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_15000.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_15000_cp.append(pot)
        f1.close()
        f2=open("./TPS_15000.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_15000_cp.append(pot)
        f2.close()
        


accepted_ind_TPS_17500_cp=[]
for line in race_accept_TPS_17500_cp:
    accepted_ind_TPS_17500_cp.append(line.split()[5])
pot_TPS_17500_cp=[]
for item in accepted_ind_TPS_17500_cp:
    if int(item)<10:
        f1=open("./TPS_17500.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_17500_cp.append(pot)
        f1.close()
        f2=open("./TPS_17500.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_17500_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_17500.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_17500_cp.append(pot)
        f1.close()
        f2=open("./TPS_17500.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_17500_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_17500.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_17500_cp.append(pot)
        f1.close()
        f2=open("./TPS_17500.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_17500_cp.append(pot)
        f2.close()
        




accepted_ind_TPS_20000_cp=[]
for line in race_accept_TPS_20000_cp:
    accepted_ind_TPS_20000_cp.append(line.split()[5])
pot_TPS_20000_cp=[]
for item in accepted_ind_TPS_20000_cp:
    if int(item)<10:
        f1=open("./TPS_20000.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_20000_cp.append(pot)
        f1.close()
        f2=open("./TPS_20000.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_20000_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_20000.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_20000_cp.append(pot)
        f1.close()
        f2=open("./TPS_20000.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_20000_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_20000.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_20000_cp.append(pot)
        f1.close()
        f2=open("./TPS_20000.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_20000_cp.append(pot)
        f2.close()
        




accepted_ind_TPS_22500_cp=[]
for line in race_accept_TPS_22500_cp:
    accepted_ind_TPS_22500_cp.append(line.split()[5])
pot_TPS_22500_cp=[]
for item in accepted_ind_TPS_22500_cp:
    if int(item)<10:
        f1=open("./TPS_22500.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_22500_cp.append(pot)
        f1.close()
        f2=open("./TPS_22500.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_22500_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_22500.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_22500_cp.append(pot)
        f1.close()
        f2=open("./TPS_22500.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_22500_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_22500.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_22500_cp.append(pot)
        f1.close()
        f2=open("./TPS_22500.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_22500_cp.append(pot)
        f2.close()
        




accepted_ind_TPS_25000=[]
for line in race_accept_TPS_25000:
    accepted_ind_TPS_25000.append(line.split()[5])
pot_TPS_25000=[]
for item in accepted_ind_TPS_25000:
    if int(item)<10:
        f1=open("./TPS_25000.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_25000.append(pot)
        f1.close()
        f2=open("./TPS_25000.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_25000.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_25000.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_25000.append(pot)
        f1.close()
        f2=open("./TPS_25000.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_25000.append(pot)
        f2.close()
    else:
        f1=open("./TPS_25000.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_25000.append(pot)
        f1.close()
        f2=open("./TPS_25000.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_25000.append(pot)
        f2.close()    
        
        
        


accepted_ind_TPS_27500_cp=[]
for line in race_accept_TPS_27500_cp:
    accepted_ind_TPS_27500_cp.append(line.split()[5])
pot_TPS_27500_cp=[]
for item in accepted_ind_TPS_27500_cp:
    if int(item)<10:
        f1=open("./TPS_27500.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_27500_cp.append(pot)
        f1.close()
        f2=open("./TPS_27500.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_27500_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_27500.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_27500_cp.append(pot)
        f1.close()
        f2=open("./TPS_27500.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_27500_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_27500.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_27500_cp.append(pot)
        f1.close()
        f2=open("./TPS_27500.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_27500_cp.append(pot)
        f2.close()
        




accepted_ind_TPS_30000=[]
for line in race_accept_TPS_30000:
    accepted_ind_TPS_30000.append(line.split()[5])
pot_TPS_30000=[]
for item in accepted_ind_TPS_30000:
    if int(item)<10:
        f1=open("./TPS_30000.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_30000.append(pot)
        f1.close()
        f2=open("./TPS_30000.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_30000.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_30000.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_30000.append(pot)
        f1.close()
        f2=open("./TPS_30000.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_30000.append(pot)
        f2.close()
    else:
        f1=open("./TPS_30000.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_30000.append(pot)
        f1.close()
        f2=open("./TPS_30000.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_30000.append(pot)
        f2.close()    
        
        


accepted_ind_TPS_32500=[]
for line in race_accept_TPS_32500:
    accepted_ind_TPS_32500.append(line.split()[5])
pot_TPS_32500=[]
for item in accepted_ind_TPS_32500:
    if int(item)<10:
        f1=open("./TPS_32500.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_32500.append(pot)
        f1.close()
        f2=open("./TPS_32500.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_32500.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_32500.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_32500.append(pot)
        f1.close()
        f2=open("./TPS_32500.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_32500.append(pot)
        f2.close()
    else:
        f1=open("./TPS_32500.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_32500.append(pot)
        f1.close()
        f2=open("./TPS_32500.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_32500.append(pot)
        f2.close()    
        
        
        
                


accepted_ind_TPS_35000_cp=[]
for line in race_accept_TPS_35000_cp:
    accepted_ind_TPS_35000_cp.append(line.split()[5])
pot_TPS_35000_cp=[]
for item in accepted_ind_TPS_35000_cp:
    if int(item)<10:
        f1=open("./TPS_35000.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_35000_cp.append(pot)
        f1.close()
        f2=open("./TPS_35000.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_35000_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_35000.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_35000_cp.append(pot)
        f1.close()
        f2=open("./TPS_35000.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_35000_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_35000.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_35000_cp.append(pot)
        f1.close()
        f2=open("./TPS_35000.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_35000_cp.append(pot)
        f2.close()
        



accepted_ind_TPS_37500_cp=[]
for line in race_accept_TPS_37500_cp:
    accepted_ind_TPS_37500_cp.append(line.split()[5])
pot_TPS_37500_cp=[]
for item in accepted_ind_TPS_37500_cp:
    if int(item)<10:
        f1=open("./TPS_37500.cp.BOLAS-2/gn_00"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_37500_cp.append(pot)
        f1.close()
        f2=open("./TPS_37500.cp.BOLAS-2/gn_00"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_37500_cp.append(pot)
        f2.close()
    elif int(item)<100:
        f1=open("./TPS_37500.cp.BOLAS-2/gn_0"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_37500_cp.append(pot)
        f1.close()
        f2=open("./TPS_37500.cp.BOLAS-2/gn_0"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_37500_cp.append(pot)
        f2.close()
    else:
        f1=open("./TPS_37500.cp.BOLAS-2/gn_"+item+"_xpr1.eptot-vs-time.dat","r")
        lines=f1.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_37500_cp.append(pot)
        f1.close()
        f2=open("./TPS_37500.cp.BOLAS-2/gn_"+item+"_xpr2.eptot-vs-time.dat","r")
        lines=f2.readlines()
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                pot = float(tokens[1])
                pot_TPS_37500_cp.append(pot)
        f2.close()
        


        
#race_files = ["./TPS_2500.cp.BOLAS-2/race_history.log","./TPS_5000.cp.BOLAS-2/race_history.log","./TPS_7500.BOLAS-2/race_history.log","./TPS_10000.BOLAS-2/race_history.log","./TPS_12500.BOLAS-2/race_history.log","./TPS_15000.cp.BOLAS-2/race_history.log","./TPS_17500.cp.BOLAS-2/race_history.log","./TPS_20000.cp.BOLAS-2/race_history.log","./TPS_22500.cp.BOLAS-2/race_history.log","./TPS_25000.BOLAS-2/race_history.log","./TPS_27500.cp.BOLAS-2/race_history.log","./TPS_30000.BOLAS-2/race_history.log","./TPS_32500.BOLAS-2/race_history.log","./TPS_35000.cp.BOLAS-2/race_history.log","./TPS_37500.cp.BOLAS-2/race_history.log"]



#scale = []
#accept_traj = []
#traj_chstp=[]
#dlsh=[]


#for race_log in race_files:
#    f=open(race_log,"r")
#    lines = f.readlines()
#    #scale = []
#    #accept_traj = []
#    for line in lines:
#        if "TRAJS ACCEPTED" in line:
#            accept_traj.append(line.split()[5])
#            traj_chstp.append(int(lines[lines.index(line)-6].split()[1]))
#            traj_prev = lines[lines.index(line)-7].split()[1]
#            if traj_prev.split("_")[2]=="xpr2":
#                scale.append([-1,1])
#            else:
#                scale.append([1,-1])
#            if lines[lines.index(line)-13].split()[0]=="SHIFTING":
#                dlsh.append(int(lines[lines.index(line)-12].split()[3]))
#            else:
#                dlsh.append(0)
#    f.close()


#traj_order = []
#traj_order_byf = []
#frame_order = []
#dist_data_update=[]
#traj_length=[]

#coord_data_update=[]
#for i in range(int(len(dist_data)/2)):
#    if int(accept_traj[i]) < 10:
#        traj_name ="gn_00"+accept_traj[i]+"_xpr"
#    elif int(accept_traj[i]) >= 10 and int(accept_traj[i]) < 100:
#        traj_name ="gn_0"+accept_traj[i]+"_xpr"
#    else:
#        traj_name ="gn_"+accept_traj[i]+"_xpr"
#    if scale[i]==[1,-1]:
#        dist_subset=np.append(dist_data[i*2+1][::-1],dist_data[i*2],axis=0)
#        #coord_subset=np.append(coord_data[i*2+1][::-1],coord_data[i*2],axis=0)
#        traj = [traj_name+"2",traj_name+"1"]
#        traj_order = traj_order + traj
#        frame_a = list(range(len(dist_data[i*2+1])))
#        frame_a.reverse()
#        frame_b = list(range(len(dist_data[i*2])))
#        frame_order = frame_order + [fa+1 for fa in frame_a] + [fb+1 for fb in frame_b]
#        traj_order_byf = traj_order_byf + [traj_name+"2"]*len(dist_data[i*2+1]) + [traj_name+"1"]*len(dist_data[i*2])
#    else:
#        dist_subset=np.append(dist_data[i*2][::-1],dist_data[i*2+1],axis=0)
#        #coord_subset=np.append(coord_data[i*2][::-1],coord_data[i*2+1],axis=0)
#        traj = [traj_name+"1",traj_name+"2"]
#        traj_order = traj_order + traj
#        frame_a = list(range(len(dist_data[i*2])))
#        frame_a.reverse()
#        frame_b = list(range(len(dist_data[i*2+1])))
#        frame_order = frame_order + [fa+1 for fa in frame_a] + [fb+1 for fb in frame_b]
#        traj_order_byf = traj_order_byf + [traj_name+"1"]*len(dist_data[i*2]) + [traj_name+"2"]*len(dist_data[i*2+1])
#    dist_data_update.append(dist_subset)
#    traj_length.append([len(dist_data[i*2]),len(dist_data[i*2+1])])
#    #coord_data_update.append(coord_subset)

#dist_data_concatenated = np.concatenate(dist_data_update)





# Parameters
K = 14 # number of umbrellas 
temperature=298.15
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K
T_k = np.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
nbins = 100 # number of bins for 1D PMF
N_k = [len(w) for w in [dist_data_2500_concatenated, dist_data_5000_concatenated, dist_data_7500_concatenated,dist_data_10000_concatenated,dist_data_12500_concatenated,dist_data_15000_concatenated,dist_data_17500_concatenated,dist_data_20000_concatenated,dist_data_22500_concatenated,dist_data_25000_concatenated,dist_data_27500_concatenated,dist_data_30000_concatenated,dist_data_35000_concatenated,dist_data_37500_concatenated]] 
N_max=max(N_k) # maximum number of snapshots/window


#CV1 = dist_data_concatenated[:,0]
#CV2 = dist_data_concatenated[:,1]
#cv1_min=min(CV1)
#cv1_max=max(CV1)
#cv2_min=min(CV2)
#cv2_max=max(CV2)
beta_k = 1.0/(kB*T_k)
nbins=100


pot_data=[pot_TPS_2500_cp,pot_TPS_5000_cp,pot_TPS_7500,pot_TPS_10000,pot_TPS_12500,pot_TPS_15000_cp,pot_TPS_17500_cp,pot_TPS_20000_cp,pot_TPS_22500_cp,pot_TPS_25000,pot_TPS_27500_cp,pot_TPS_30000,pot_TPS_35000_cp,pot_TPS_37500_cp]
dist_data_list_all=[dist_data_2500_concatenated, dist_data_5000_concatenated, dist_data_7500_concatenated,dist_data_10000_concatenated,dist_data_12500_concatenated,dist_data_15000_concatenated,dist_data_17500_concatenated,dist_data_20000_concatenated,dist_data_22500_concatenated,dist_data_25000_concatenated,dist_data_27500_concatenated,dist_data_30000_concatenated,dist_data_35000_concatenated,dist_data_37500_concatenated] 


U_kn = np.zeros([K,N_max], np.float64)
x_kn = np.zeros([K, N_max], np.float64)
y_kn = np.zeros([K, N_max], np.float64)
for i in range(K):
    for t in range(len(pot_data[i])):
        U_kn[i,t]= pot_data[i][t]
        x_kn[i,t] = dist_data_list_all[i][t,0]
        y_kn[i,t] = dist_data_list_all[i][t,1]

u_kln = np.zeros([K,K,N_max], np.float32) # u_kln[k,l,n] is reduced potential energy of trajectory segment n of temperature k evaluated at temperature l
for k in range(K):
   for l in range(K):
      u_kln[k,l,0:N_k[k]] = beta_k[l] * U_kn[k,0:N_k[k]]
      
#Ubias = np.zeros((len(dist_data_concatenated),2))
#a_para=np.ones(14)/sum(np.ones(14))
#weights = np.concatenate([a_para[0]*np.ones(len(dist_data_2500_concatenated)), a_para[1]*np.ones(len(dist_data_5000_concatenated)),a_para[2]*np.ones(len(dist_data_7500_concatenated)),a_para[3]*np.ones(len(dist_data_10000_concatenated)),a_para[4]*np.ones(len(dist_data_12500_concatenated)),a_para[5]*np.ones(len(dist_data_15000_concatenated)),a_para[6]*np.ones(len(dist_data_17500_concatenated)),a_para[7]*np.ones(len(dist_data_20000_concatenated)),a_para[8]*np.ones(len(dist_data_22500_concatenated)),a_para[9]*np.ones(len(dist_data_25000_concatenated)),a_para[10]*np.ones(len(dist_data_27500_concatenated)),a_para[11]*np.ones(len(dist_data_30000_concatenated)),a_para[12]*np.ones(len(dist_data_35000_concatenated)),a_para[13]*np.ones(len(dist_data_37500_concatenated))])

z_bolas_scipy, xHedges_bolas, yHedges_bolas, binnumber_bolas = scipy.stats.binned_statistic_2d(*dist_data_concatenated.T, None, 'count', bins=nbins, expand_binnumbers=True)

set_bins=100

mask_kn = np.zeros([K,N_max], dtype=np.bool_)
for k in range(0,K):
   mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
indices = np.where(mask_kn)

    

dx = xHedges_bolas[1]-xHedges_bolas[0]
dy = yHedges_bolas[1]-yHedges_bolas[0]
# Assign torsion bins
bin_kn = np.zeros([K,N_max], np.int16) # bin_kn[k,n] is the index of which histogram bin sample n from temperature index k belongs to
nbins = 0
bin_counts = list()
bin_centers = list() # bin_centers[i] is a (phi,psi) tuple that gives the center of bin i
for i in range(set_bins):
   for j in range(set_bins):
      # Determine (phi,psi) of bin center.
      x = xHedges_bolas[0] + dx * (i + 0.5)
      y = yHedges_bolas[0] + dy * (j + 0.5)

      # Determine which configurations lie in this bin.
      in_bin = (x-dx/2 <= x_kn[indices]) & (x_kn[indices] < x+dx/2) & (y-dy/2 <= y_kn[indices]) & (y_kn[indices] < y+dy/2)

      # Count number of configurations in this bin.
      bin_count = in_bin.sum()

      # Generate list of indices in bin.
      indices_in_bin = (indices[0][in_bin], indices[1][in_bin])

      if (bin_count > 0):
         # store bin (phi,psi)
         bin_centers.append( (x, y) )
         bin_counts.append( bin_count )

         # assign these conformations to the bin index
         bin_kn[indices_in_bin] = nbins

         # increment number of bins
         nbins += 1




# --- Define MBAR object, handling non-uniformity ---
mbar = MBAR(u_kln, N_k,verbose=True)  # Pass sample counts explicitly for non-uniform windows
fes = FES(U_kn, bin_kn)



print("2D PMF")
print("")
print("%8s %6s %6s %8s %10s %10s" % ('bin', 'DIST 1', 'DIST 2', 'N', 'f', 'df'))

for i in range(nbins):
   print('%8d %6.1f %6.1f %8d %10.3f %10.3f' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i], f_i[i], df_i[i]))

# --- Define PMF grid ---
#cv1_grid = np.linspace(cv1_min, cv1_max, nbins)  # Replace with your ranges and grid points
#cv2_grid = np.linspace(cv2_min, cv2_max, nbins)

# --- Calculate PMF ---
#pmf_grid = mbar.computePMF(cv1_grid, cv2_grid, Ubias)




#//def get_histogram(
#//        xall, yall, nbins=100,
#//        weights=None, avoid_zero_count=False):
#//    """Compute a two-dimensional histogram.
#//
#//    Parameters
#//    ----------
#//    xall : ndarray(T)
#//        Sample x-coordinates.
#//    yall : ndarray(T)
#//        Sample y-coordinates.
#//    nbins : int, optional, default=100
#//        Number of histogram bins used in each dimension.
#//    weights : ndarray(T), optional, default=None
#//        Sample weights; by default all samples have the same weight.
#//    avoid_zero_count : bool, optional, default=True
#//        Avoid zero counts by lifting all histogram elements to the
#//        minimum value before computing the free energy. If False,
#//        zero histogram counts would yield infinity in the free energy.
#//
#//    Returns
#//    -------
#//    x : ndarray(nbins, nbins)
#//        The bins' x-coordinates in meshgrid format.
#//    y : ndarray(nbins, nbins)
#//        The bins' y-coordinates in meshgrid format.
#//    z : ndarray(nbins, nbins)
#//        Histogram counts in meshgrid format.
#//
#//    """
#//    z, xedge, yedge = np.histogram2d(
#//        xall, yall, bins=nbins, weights=weights)
#//    x = 0.5 * (xedge[:-1] + xedge[1:])
#//    y = 0.5 * (yedge[:-1] + yedge[1:])
#//    if avoid_zero_count:
#//        z = np.maximum(z, np.min(z[z.nonzero()]))
#//    return x, y, z.T # transpose to match x/y-directions
#//
#//
#//def _to_density(z):
#//    """Normalize histogram counts.
#//
#//    Parameters
#//    ----------
#//    z : ndarray(T)
#//        Histogram counts.
#//
#//    """
#//    return z / float(z.sum())
#//
#//
#//def _to_free_energy(z, minener_zero=False):
#//    """Compute free energies from histogram counts.
#//
#//    Parameters
#//    ----------
#//    z : ndarray(T)
#//        Histogram counts.
#//    minener_zero : boolean, optional, default=False
#//        Shifts the energy minimum to zero.
#//
#//    Returns
#//    -------
#//    free_energy : ndarray(T)
#//        The free energy values in units of kT.
#//
#//    """
#//    pi = _to_density(z)
#//    free_energy = np.inf * np.ones(shape=z.shape)
#//    nonzero = pi.nonzero()
#//    free_energy[nonzero] = -np.log(pi[nonzero])
#//    if minener_zero:
#//        free_energy[nonzero] -= np.min(free_energy[nonzero])
#//    return free_energy
#//
#//
#//
#//nbins=100
#//minener_zero=True
#//kT=1.0
#//weights=None
#//avoid_zero_count=False
#//minener_zero=True
#//
#//
#//#x, y, z = get_histogram(xall, yall, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//#f = _to_free_energy(z, minener_zero=minener_zero) * kT
#//
#//
#//x_2500, y_2500, z_2500 = get_histogram(*dist_data_2500_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_2500_raw = _to_free_energy(z_2500, minener_zero=minener_zero) * kT
#//
#//x_5000, y_5000, z_5000 = get_histogram(*dist_data_5000_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_5000_raw = _to_free_energy(z_5000, minener_zero=minener_zero) * kT
#//
#//x_7500, y_7500, z_7500 = get_histogram(*dist_data_7500_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_7500_raw = _to_free_energy(z_7500, minener_zero=minener_zero) * kT
#//
#//x_10000, y_10000, z_10000 = get_histogram(*dist_data_10000_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_10000_raw = _to_free_energy(z_10000, minener_zero=minener_zero) * kT
#//
#//x_12500, y_12500, z_12500 = get_histogram(*dist_data_12500_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_12500_raw = _to_free_energy(z_12500, minener_zero=minener_zero) * kT
#//
#//x_15000, y_15000, z_15000 = get_histogram(*dist_data_15000_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_15000_raw = _to_free_energy(z_15000, minener_zero=minener_zero) * kT
#//
#//x_17500, y_17500, z_17500 = get_histogram(*dist_data_17500_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_17500_raw = _to_free_energy(z_17500, minener_zero=minener_zero) * kT
#//
#//x_20000, y_20000, z_20000 = get_histogram(*dist_data_20000_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_20000_raw = _to_free_energy(z_20000, minener_zero=minener_zero) * kT
#//
#//x_22500, y_22500, z_22500 = get_histogram(*dist_data_22500_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_22500_raw = _to_free_energy(z_22500, minener_zero=minener_zero) * kT
#//
#//x_25000, y_25000, z_25000 = get_histogram(*dist_data_25000_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#//f_25000_raw = _to_free_energy(z_25000, minener_zero=minener_zero) * kT
#//
#x_27500, y_27500, z_27500 = get_histogram(*dist_data_27500_concatenated.T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)
#f_27500_raw = _to_free_energy(z_27500, minener_zero=minener_zero) * kT




fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True, sharey=True)
# the * operator used in a function call is used to unpack
# the iterable variable into its single elements. 
pyemma.plots.plot_density(*dist_data_concatenated.T, ax=axes[0])
pyemma.plots.plot_free_energy(*dist_data_concatenated.T, ax=axes[1], legacy=False)
for ax in axes.flat:
    ax.set_xlabel('DIST 1')
    ax.set_aspect('equal')
axes[0].set_ylabel('DIST 2')
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.png")



data_index=list(range(len(dist_data_concatenated)))
data_sele=[i for i in data_index if i not in [2989755,2989756]]

nclusters=75
cluster_list=list(range(nclusters))

fig, ax = plt.subplots(figsize=(7, 6))
# the * operator used in a function call is used to unpack
# the iterable variable into its single elements. 
#pyemma.plots.plot_density(*dist_data_concatenated.T, ax=axes[0])
pyemma.plots.plot_free_energy(*dist_data_concatenated[data_sele,:].T, ax=ax, legacy=False)
ax.scatter(dist_data_center[:,0],dist_data_center[:,1], s=10, c='white', marker="o")
#for i, txt in enumerate(cluster_list):
#    ax.annotate(txt, (dist_data_center[i,0], dist_data_center[i,1]))
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.1113.label_center.png",dpi=1200)
#fig.savefig("DIST_FEL_density.BOLAS.0928.png",dpi=1200)


fig, ax = plt.subplots(figsize=(7, 6))
# the * operator used in a function call is used to unpack
# the iterable variable into its single elements. 
#pyemma.plots.plot_density(*dist_data_concatenated.T, ax=axes[0])
pyemma.plots.plot_free_energy(*dist_data_concatenated.T, ax=ax, legacy=False)
#ax.scatter(dist_data_center[:,0],dist_data_center[:,1], s=10, c='white', marker="o")
#for i, txt in enumerate(cluster_list):
#    ax.annotate(txt, (dist_data_center[i,0], dist_data_center[i,1]))
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.011824.png",dpi=1200)



fig, ax = plt.subplots(figsize=(7, 6))
# the * operator used in a function call is used to unpack
# the iterable variable into its single elements. 
#pyemma.plots.plot_density(*dist_data_concatenated.T, ax=axes[0])
pyemma.plots.plot_free_energy(*dist_data_2500_concatenated.T, ax=ax, legacy=False)
#ax.scatter(dist_data_center[:,0],dist_data_center[:,1], s=10, c='white', marker="o")
#for i, txt in enumerate(cluster_list):
#    ax.annotate(txt, (dist_data_center[i,0], dist_data_center[i,1]))
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.2500D.011824.png",dpi=1200)





def get_histogram(
        xall, yall, nbins=100,
        weights=None, avoid_zero_count=False):
    """Copied from pyemma
    """
    z, xedge, yedge = np.histogram2d(
        xall, yall, bins=nbins, weights=weights)
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    if avoid_zero_count:
        z = np.maximum(z, np.min(z[z.nonzero()]))
    return x, y, z.T # transpose to match x/y-directions

def _to_density(z):
    """Copied from pyemma
    """
    return z / float(z.sum())


def _to_free_energy(z, minener_zero=False):
    """Copied from pyemma
    """
    pi = _to_density(z)
    free_energy = np.inf * np.ones(shape=z.shape)
    nonzero = pi.nonzero()
    free_energy[nonzero] = -np.log(pi[nonzero])
    if minener_zero:
        free_energy[nonzero] -= np.min(free_energy[nonzero])
    return free_energy


def _prune_kwargs(kwargs):
    """Copied from pyemma
    """
    allowed_keys = [
        'corner_mask', 'alpha', 'locator', 'extend', 'xunits',
        'yunits', 'antialiased', 'nchunk', 'hatches', 'zorder']
    ignored = [key for key in kwargs.keys() if key not in allowed_keys]
    for key in ignored:
        _warn(
            '{}={} is not an allowed optional parameter and will'
            ' be ignored'.format(key, kwargs[key]))
        kwargs.pop(key, None)
    return kwargs

def plot_map(
        x, y, z, ax=None, cmap=None,
        ncontours=100, vmin=None, vmax=None, levels=None,
        cbar=True, cax=None, cbar_label=None,
        cbar_orientation='vertical', norm=None,
        **kwargs):
    """Copied from pyemma
    """
    import matplotlib.pyplot as _plt
    if ax is None:
        fig, ax = _plt.subplots()
    else:
        fig = ax.get_figure()
    mappable = ax.contourf(
        x, y, z, ncontours, norm=norm,
        vmin=vmin, vmax=vmax, cmap=cmap,
        levels=levels, **_prune_kwargs(kwargs))
    misc = dict(mappable=mappable)
    if cbar_orientation not in ('horizontal', 'vertical'):
        raise ValueError(
            'cbar_orientation must be "horizontal" or "vertical"')
    if cbar:
        if cax is None:
            cbar_ = fig.colorbar(
                mappable, ax=ax, orientation=cbar_orientation)
        else:
            cbar_ = fig.colorbar(
                mappable, cax=cax, orientation=cbar_orientation)
        if cbar_label is not None:
            cbar_.set_label(cbar_label)
        misc.update(cbar=cbar_)
    return fig, ax, misc


def plot_free_energy(
        xall, yall, weights=None, ax=None, nbins=100, ncontours=100,
        offset=-1, avoid_zero_count=False, minener_zero=True, kT=1.0,
        vmin=None, vmax=None, cmap='nipy_spectral', cbar=True,
        cbar_label='free energy / kT', cax=None, levels=None,
        legacy=True, ncountours=None, cbar_orientation='vertical',
        **kwargs):
    """Copied from pyemma
    """
    if legacy:
        _warn(
            'Legacy mode is deprecated is will be removed in the'
            ' next major release. Until then use legacy=False',
            DeprecationWarning)
        cmap = _get_cmap(cmap)
        if offset != -1:
            _warn(
                'Parameter offset is deprecated and will be ignored',
                DeprecationWarning)
        if ncountours is not None:
            _warn(
                'Parameter ncountours is deprecated;'
                ' use ncontours instead',
                DeprecationWarning)
            ncontours = ncountours
        if vmin is None:
            vmin = 0.0
    else:
        if offset != -1:
            raise ValueError(
                'Parameter offset is not allowed outside legacy mode')
        if ncountours is not None:
            raise ValueError(
                'Parameter ncountours is not allowed outside'
                ' legacy mode; use ncontours instead')
    x, y, z = get_histogram(
        xall, yall, nbins=nbins, weights=weights,
        avoid_zero_count=avoid_zero_count)
    f = _to_free_energy(z, minener_zero=minener_zero) * kT
    fig, ax, misc = plot_map(
        x, y, f, ax=ax, cmap=cmap,
        ncontours=ncontours, vmin=vmin, vmax=vmax, levels=levels,
        cbar=cbar, cax=cax, cbar_label=cbar_label,
        cbar_orientation=cbar_orientation, norm=None,
        **kwargs)
    if legacy:
        return fig, ax
    return fig, ax, misc

def get_histogram_grid(
        xedge, yedge, nbins=100,
        weights=None, avoid_zero_count=False):
    """Copied from pyemma; Only x, y part
    """
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    return x, y # transpose to match x/y-directions

def plot_free_energy_simplify(
        xedge, yedge, P_u,weights=None, ax=None, nbins=100, ncontours=100,
        offset=-1, avoid_zero_count=False, minener_zero=True, kT=1.0,
        vmin=None, vmax=None, cmap='nipy_spectral', cbar=True,
        cbar_label='free energy / kT', cax=None, levels=None,
        legacy=True, ncountours=None, cbar_orientation='vertical',
        **kwargs):
    """Copied from pyemma; Only last part
    """
    x, y=get_histogram_grid(xedge,yedge)
    free_energy = np.inf * np.ones(shape=P_u.shape)
    nonzero = P_u.nonzero()
    free_energy[nonzero] = -np.log(P_u[nonzero])
    if minener_zero:
        free_energy[nonzero] -= np.min(free_energy[nonzero])
    fig, ax, misc = plot_map(
        x, y, free_energy, ax=ax, cmap=cmap,
        ncontours=ncontours, vmin=vmin, vmax=vmax, levels=levels,
        cbar=cbar, cax=cax, cbar_label=cbar_label,
        cbar_orientation=cbar_orientation, norm=None,
        **kwargs)
    if legacy:
        return fig, ax
    return fig, ax, misc



fig, ax = plt.subplots(figsize=(7, 6))
plot_free_energy(*dist_data_concatenated.T, ax=ax, legacy=False)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.011824.compare_DIST_FEL_density.png")


fig, ax = plt.subplots(figsize=(7, 6))
plot_free_energy(*dist_data_2500_concatenated.T, ax=ax, legacy=False)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.011824.test_2500.compare_DIST_FEL_density.png")



fig, ax = plt.subplots(figsize=(7, 6))
plot_free_energy_simplify(xHedges_bolas,yHedges_bolas,P_u, ax=ax, legacy=False)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.011824.P_u_test.compare_DIST_FEL_density.png")


pdb_center = md.load("Cluster1.center.pdb")
files_center = glob.glob("Cluster.center.all.pyemma.nc")
dist_center_feat=pyemma.coordinates.featurizer(pdb_center)                            
dist_center_feat.add_distances(np.array([[105,2416],[66,1956]]))   # 3P [66] - 62P [1956]; 4O4' [105] - 76N1 [2416]             
dist_data_center=pyemma.coordinates.load(files_center, features=dist_center_feat)



coord_feat = pyemma.coordinates.featurizer(pdb) 
P_ind = [33,67,97,130,160,190,224,547,580,613,646,680,711,742,773,807,837,868,898,928,961,992,1025,1056,1087,1211,2052,2085,2115,2148,2181,2212,2246,2280,2314,2345,2375,2405,2435] 
P_index = list(map(lambda x:x-1,P_ind)) # P in resi 1-7+18-36+65-77
coord_feat.add_selection(P_index) #3*39
coord_feat.add_residue_COM(list(range(1,8))+list(range(18,37))+list(range(65,78))) #
coord_feat.add_distances(coord_feat.pairs(P_ind,excluded_neighbors=2),periodic=True) #975
coord_data = pyemma.coordinates.load(files, features=coord_feat)



race_files = ["./TPS_2500.cp/race_history.log","./TPS_7500/race_history.log","./TPS_10000/race_history.log","./TPS_12500/race_history.log","./TPS_15000.cp/race_history.log","./TPS_17500.cp/race_history.log","./TPS_20000.cp/race_history.log","./TPS_22500.cp/race_history.log","./TPS_25000/race_history.log","./TPS_27500.cp/race_history.log","./TPS_30000/race_history.log","./TPS_35000.cp/race_history.log", "./TPS_37500.cp/race_history.log","./TPS_2500.cp-2/race_history.log", "./TPS_5000.cp-2/race_history.log", "./TPS_7500.cp-2/race_history.log", "./TPS_10000.cp-2/race_history.log", "./TPS_12500.cp-2/race_history.log", "./TPS_15000.cp-2/race_history.log","./TPS_17500.cp-2/race_history.log", "./TPS_20000.cp-2/race_history.log"]


x, y, z = get_histogram(*dist_data_concatenated[data_sele,:].T, nbins=nbins, weights=weights,avoid_zero_count=avoid_zero_count)

def plot_map(
        x, y, z, ax=None, cmap=None,
        ncontours=100, vmin=None, vmax=None, levels=None,
        cbar=True, cax=None, cbar_label=None,
        cbar_orientation='vertical', norm=None):
    """Plot a two-dimensional map from data on a grid.

    Parameters
    ----------
    x : ndarray(T)
        Binned x-coordinates.
    y : ndarray(T)
        Binned y-coordinates.
    z : ndarray(T)
        Binned z-coordinates.
    ax : matplotlib.Axes object, optional, default=None
        The ax to plot to; if ax=None, a new ax (and fig) is created.
    cmap : matplotlib colormap, optional, default=None
        The color map to use.
    ncontours : int, optional, default=100
        Number of contour levels.
    vmin : float, optional, default=None
        Lowest z-value to be plotted.
    vmax : float, optional, default=None
        Highest z-value to be plotted.
    levels : iterable of float, optional, default=None
        Contour levels to plot.
    cbar : boolean, optional, default=True
        Plot a color bar.
    cax : matplotlib.Axes object, optional, default=None
        Plot the colorbar into a custom axes object instead of
        stealing space from ax.
    cbar_label : str, optional, default=None
        Colorbar label string; use None to suppress it.
    cbar_orientation : str, optional, default='vertical'
        Colorbar orientation; choose 'vertical' or 'horizontal'.
    norm : matplotlib norm, optional, default=None
        Use a norm when coloring the contour plot.

    Optional parameters for contourf (**kwargs)
    -------------------------------------------
    corner_mask : boolean, optional
        Enable/disable corner masking, which only has an effect if
        z is a masked array. If False, any quad touching a masked
        point is masked out. If True, only the triangular corners
        of quads nearest those points are always masked out, other
        triangular corners comprising three unmasked points are
        contoured as usual.
        Defaults to rcParams['contour.corner_mask'], which
        defaults to True.
    alpha : float
        The alpha blending value.
    locator : [ None | ticker.Locator subclass ]
        If locator is None, the default MaxNLocator is used. The
        locator is used to determine the contour levels if they are
        not given explicitly via the levels argument.
    extend : [ ‘neither’ | ‘both’ | ‘min’ | ‘max’ ]
        Unless this is ‘neither’, contour levels are automatically
        added to one or both ends of the range so that all data are
        included. These added ranges are then mapped to the special
        colormap values which default to the ends of the
        colormap range, but can be set via
        matplotlib.colors.Colormap.set_under() and
        matplotlib.colors.Colormap.set_over() methods.
    xunits, yunits : [ None | registered units ]
        Override axis units by specifying an instance of a
        matplotlib.units.ConversionInterface.
    antialiased : boolean, optional
        Enable antialiasing, overriding the defaults. For filled
        contours, the default is True. For line contours, it is
        taken from rcParams[‘lines.antialiased’].
    nchunk : [ 0 | integer ]
        If 0, no subdivision of the domain. Specify a positive
        integer to divide the domain into subdomains of nchunk by
        nchunk quads. Chunking reduces the maximum length of polygons
        generated by the contouring algorithm which reduces the
        rendering workload passed on to the backend and also requires
        slightly less RAM. It can however introduce rendering
        artifacts at chunk boundaries depending on the backend, the
        antialiased flag and value of alpha.
    hatches :
        A list of cross hatch patterns to use on the filled areas.
        If None, no hatching will be added to the contour. Hatching
        is supported in the PostScript, PDF, SVG and Agg backends
        only.
    zorder : float
        Set the zorder for the artist. Artists with lower zorder
        values are drawn first.

    Returns
    -------
    fig : matplotlib.Figure object
        The figure in which the used ax resides.
    ax : matplotlib.Axes object
        The ax in which the map was plotted.
    misc : dict
        Contains a matplotlib.contour.QuadContourSet 'mappable' and,
        if requested, a matplotlib.Colorbar object 'cbar'.

    """
    import matplotlib.pyplot as _plt
    if ax is None:
        fig, ax = _plt.subplots()
    else:
        fig = ax.get_figure()
    mappable = ax.contourf(
        x, y, z, ncontours, norm=norm,
        vmin=vmin, vmax=vmax, cmap=cmap)
    misc = dict(mappable=mappable)
    if cbar_orientation not in ('horizontal', 'vertical'):
        raise ValueError(
            'cbar_orientation must be "horizontal" or "vertical"')
    if cbar:
        if cax is None:
            cbar_ = fig.colorbar(
                mappable, ax=ax, orientation=cbar_orientation)
        else:
            cbar_ = fig.colorbar(
                mappable, cax=cax, orientation=cbar_orientation)
        if cbar_label is not None:
            cbar_.set_label(cbar_label)
        misc.update(cbar=cbar_)
    return fig, ax, misc

    
#fig, ax, misc = plot_map(x, y, f, ax=ax, cmap=cmap, ncontours=ncontours, vmin=vmin, vmax=vmax, levels=levels, cbar=cbar, cax=cax, cbar_label=cbar_label, cbar_orientation=cbar_orientation, norm=None)





def get_histogram(
        xall, yall, nbins=100,
        weights=None, avoid_zero_count=False):
    """Compute a two-dimensional histogram.

    Parameters
    ----------
    xall : ndarray(T)
        Sample x-coordinates.
    yall : ndarray(T)
        Sample y-coordinates.
    nbins : int, optional, default=100
        Number of histogram bins used in each dimension.
    weights : ndarray(T), optional, default=None
        Sample weights; by default all samples have the same weight.
    avoid_zero_count : bool, optional, default=True
        Avoid zero counts by lifting all histogram elements to the
        minimum value before computing the free energy. If False,
        zero histogram counts would yield infinity in the free energy.

    Returns
    -------
    x : ndarray(nbins, nbins)
        The bins' x-coordinates in meshgrid format.
    y : ndarray(nbins, nbins)
        The bins' y-coordinates in meshgrid format.
    z : ndarray(nbins, nbins)
        Histogram counts in meshgrid format.

    """
    z, xedge, yedge = np.histogram2d(
        xall, yall, bins=nbins, weights=weights)
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    if avoid_zero_count:
        z = np.maximum(z, np.min(z[z.nonzero()]))
    return x, y, z.T # transpose to match x/y-directions


def _to_density(z):
    """Normalize histogram counts.

    Parameters
    ----------
    z : ndarray(T)
        Histogram counts.

    """
    return z / float(z.sum())


def _to_free_energy(z, minener_zero=False):
    """Compute free energies from histogram counts.

    Parameters
    ----------
    z : ndarray(T)
        Histogram counts.
    minener_zero : boolean, optional, default=False
        Shifts the energy minimum to zero.

    Returns
    -------
    free_energy : ndarray(T)
        The free energy values in units of kT.

    """
    pi = _to_density(z)
    free_energy = np.inf * np.ones(shape=z.shape)
    nonzero = pi.nonzero()
    free_energy[nonzero] = -np.log(pi[nonzero])
    if minener_zero:
        free_energy[nonzero] -= np.min(free_energy[nonzero])
    return free_energy



nbins=100
minener_zero=True
kT=1.0


z_bolas, xedge_bolas, yedge_bolas = np.histogram2d(*dist_data_concatenated.T, bins=nbins)
z_bolas_scipy, xHedges_bolas, yHedges_bolas, binnumber_bolas = scipy.stats.binned_statistic_2d(*dist_data_concatenated.T, None, 'count', bins=nbins, expand_binnumbers=True)
f_bolas = _to_free_energy(z_bolas, minener_zero=minener_zero) * kT

f=open("dist_energy_bolas.update.pickle",'wb')
pickle.dump(z_bolas_scipy,f)
pickle.dump(xHedges_bolas,f)
pickle.dump(yHedges_bolas,f)
pickle.dump(binnumber_bolas,f)
pickle.dump(f_bolas,f)
f.close()



histogram_2500, _, _ = np.histogram2d(dist_data_2500_concatenated.T[0,:],dist_data_2500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_5000, _, _ = np.histogram2d(dist_data_5000_concatenated.T[0,:],dist_data_5000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_7500, _, _ = np.histogram2d(dist_data_7500_concatenated.T[0,:],dist_data_7500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_10000, _, _ = np.histogram2d(dist_data_10000_concatenated.T[0,:],dist_data_10000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_12500, _, _ = np.histogram2d(dist_data_15000_concatenated.T[0,:],dist_data_15000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_15000, _, _ = np.histogram2d(dist_data_12500_concatenated.T[0,:],dist_data_12500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_17500, _, _ = np.histogram2d(dist_data_17500_concatenated.T[0,:],dist_data_17500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_20000, _, _ = np.histogram2d(dist_data_20000_concatenated.T[0,:],dist_data_20000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_22500, _, _ = np.histogram2d(dist_data_22500_concatenated.T[0,:],dist_data_22500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_25000, _, _ = np.histogram2d(dist_data_25000_concatenated.T[0,:],dist_data_25000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_27500, _, _ = np.histogram2d(dist_data_27500_concatenated.T[0,:],dist_data_27500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_30000, _, _ = np.histogram2d(dist_data_30000_concatenated.T[0,:],dist_data_30000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_35000, _, _ = np.histogram2d(dist_data_35000_concatenated.T[0,:],dist_data_35000_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])
histogram_37500, _, _ = np.histogram2d(dist_data_37500_concatenated.T[0,:],dist_data_37500_concatenated.T[1,:], bins=[xedge_bolas, yedge_bolas])



f_bolas_2500 = _to_free_energy(histogram_2500, minener_zero=minener_zero) * kT
f_bolas_5000 = _to_free_energy(histogram_5000, minener_zero=minener_zero) * kT
f_bolas_7500 = _to_free_energy(histogram_7500, minener_zero=minener_zero) * kT
f_bolas_10000 = _to_free_energy(histogram_10000, minener_zero=minener_zero) * kT
f_bolas_12500 = _to_free_energy(histogram_12500, minener_zero=minener_zero) * kT
f_bolas_15000 = _to_free_energy(histogram_15000, minener_zero=minener_zero) * kT
f_bolas_17500 = _to_free_energy(histogram_17500, minener_zero=minener_zero) * kT
f_bolas_20000 = _to_free_energy(histogram_20000, minener_zero=minener_zero) * kT
f_bolas_22500 = _to_free_energy(histogram_22500, minener_zero=minener_zero) * kT
f_bolas_25000 = _to_free_energy(histogram_25000, minener_zero=minener_zero) * kT
f_bolas_27500 = _to_free_energy(histogram_27500, minener_zero=minener_zero) * kT
f_bolas_30000 = _to_free_energy(histogram_30000, minener_zero=minener_zero) * kT
f_bolas_35000 = _to_free_energy(histogram_35000, minener_zero=minener_zero) * kT
f_bolas_37500 = _to_free_energy(histogram_37500, minener_zero=minener_zero) * kT




def wham_reweighting(histograms, beta, dh_x,dh_y,max_iterations=1000, convergence_threshold=1e-6):
    num_windows = len(histograms)
    num_bins = histograms[0].shape[0]

    delta_F = np.ones(num_windows)
    a_i=np.ones(num_windows)
    p_i=a_i/a_i.sum()
    P_u_i=list()
    for i in range(num_windows):
        P_u_i.append(histograms[i]/histograms[i].sum())
    P_u=np.sum(P_u_i,axis=0)
    convergence = False

    for iteration in range(max_iterations):
        old_delta_F = delta_F.copy()

        # Reweighting step
        for i in range(num_windows):
            delta_F[i]=np.log(np.trapz(np.trapz(P_u,dx=dh_x,axis=0),dx=dh_y))/(-beta)
            a_i[i] = histograms[i].sum() * np.exp(beta * delta_F[i])       

        old_P_u=np.sum(P_u_i,axis=0)
        p_i=a_i/a_i.sum()
                
        for i in range(num_windows):                     
            P_u_i[i]=p_i[i]*P_u_i[i]

        new_P_u=np.sum(P_u_i,axis=0)
            

        # Check convergence
        if np.linalg.norm(new_P_u - old_P_u) < convergence_threshold:
            convergence = True
            break

    if convergence:
        print(f"Converged after {iteration + 1} iterations.")
    else:
        print("Did not converge within the specified iterations.")

    return a_i, P_u, delta_F


histograms=[histogram_2500,histogram_5000,histogram_7500,histogram_10000,histogram_12500,histogram_15000,histogram_17500,histogram_20000,histogram_22500,histogram_25000,histogram_27500,histogram_30000,histogram_35000,histogram_37500]
dh_x=xHedges_bolas[1]-xHedges_bolas[0] 
dh_y=yHedges_bolas[1]-yHedges_bolas[0] 

temperature=298.15
kB = 0.001985875 # Boltzmann constant in kcal/(mol⋅K)
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))



a_i, P_u, delta_F =  wham_reweighting(histograms, beta, dh_x, dh_y,max_iterations=1000, convergence_threshold=1e-20)






f=open("dist_energy_bolas.WHAM.pickle",'wb')
pickle.dump(a_i,f)
pickle.dump(P_u,f)
pickle.dump(delta_F,f)
pickle.dump(beta,f)
f.close()




def plot_free_energy_den(
        xedge, yedge, P_u, weights=None, ax=None, nbins=100, ncontours=100,
        offset=-1, avoid_zero_count=False, minener_zero=True, kT=1.0,
        vmin=None, vmax=None, cmap='nipy_spectral', cbar=True,
        cbar_label='free energy / kT', cax=None, levels=None,
        legacy=False,ncountours=None, cbar_orientation='vertical',norm=None,
        **kwargs):
    if legacy:
        _warn(
            'Legacy mode is deprecated is will be removed in the'
            ' next major release. Until then use legacy=False',
            DeprecationWarning)
        cmap = _get_cmap(cmap)
        if offset != -1:
            _warn(
                'Parameter offset is deprecated and will be ignored',
                DeprecationWarning)
        if ncountours is not None:
            _warn(
                'Parameter ncountours is deprecated;'
                ' use ncontours instead',
                DeprecationWarning)
            ncontours = ncountours
        if vmin is None:
            vmin = 0.0
    else:
        if offset != -1:
            raise ValueError(
                'Parameter offset is not allowed outside legacy mode')
        if ncountours is not None:
            raise ValueError(
                'Parameter ncountours is not allowed outside'
                ' legacy mode; use ncontours instead')
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    free_energy = np.inf * np.ones(shape=P_u.shape)
    nonzero = P_u.nonzero()
    free_energy[nonzero] = -np.log(P_u[nonzero])
    if minener_zero:
        free_energy[nonzero] -= np.min(free_energy[nonzero])
    fig, ax, misc = plot_map(
        x, y, free_energy.T, ax=ax, cmap=cmap,
        ncontours=ncontours, vmin=vmin, vmax=vmax, levels=levels,
        cbar=cbar, cax=cax, cbar_label=cbar_label,
        cbar_orientation=cbar_orientation, norm=None,
        **kwargs)
    if legacy:
        return fig, ax
    return fig, ax, misc



free_energy = np.inf * np.ones(shape=P_u.shape)
nonzero = P_u.nonzero()
free_energy[nonzero] = -np.log(P_u[nonzero])



def freeE_small(arr_2d, A):
    return [(i, j) for i, row in enumerate(arr_2d) for j, value in enumerate(row) if value < A]


free_energy_small=freeE_small(free_energy,0.5)






fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_free_energy_den(xHedges_bolas, yHedges_bolas, P_u, ax=ax)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS_WHAM.012224.png",dpi=1200)


nclusters=75
cluster_list=list(range(nclusters))
fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_free_energy_den(xHedges_bolas, yHedges_bolas, P_u, ax=ax)
ax.scatter(dist_data_center[:,0],dist_data_center[:,1], s=9, c='white', marker="o")
for i, txt in enumerate(cluster_list):
    ax.annotate(txt, (dist_data_center[i,0], dist_data_center[i,1]),fontsize=2)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS_WHAM.cluster_center.012224.png",dpi=1200)


x = 0.5 * (xHedges_bolas[:-1] + xHedges_bolas[1:])
y = 0.5 * (yHedges_bolas[:-1] + yHedges_bolas[1:])

X=np.ones(shape=(100,100))
Y=np.ones(shape=(100,100))

for i in range(0,100):
    X[i]=X[i]*x[i]
    for j in range(0,100):
        Y[j][i]=Y[j][i]*y[i]


fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")
surf=ax.plot_surface(X, Y, free_energy, cmap='nipy_spectral',edgecolor='none',lw=0.5, rstride=1, cstride=1, alpha=0.5)
ax.contour(X, Y, free_energy, 10, lw=1, colors="k", linestyles="solid")
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=16)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
ax.set_zlabel('Free energy   kcal/mol')
fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the PDF
#ax.view_init(60, 35)
fig.savefig("DIST_FEL_density.BOLAS_WHAM.3D.040124.png",dpi=1200)


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
X, Y = np.mgrid[-1:1:30j, -1:1:30j]
Z = np.sin(np.pi*X)*np.sin(np.pi*Y)
ax.plot_surface(X, Y, Z, cmap="autumn_r", lw=0.5, rstride=1, cstride=1, alpha=0.5)
ax.contour(X, Y, Z, 10, lw=3, cmap="autumn_r", linestyles="solid", offset=-1)
ax.contour(X, Y, Z, 10, lw=3, colors="k", linestyles="solid")
fig.savefig("test3D.040124.png",dpi=1200)



fig = plt.figure(figsize=(7, 6))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(xx, yy, f, rstride=1, cstride=1, cmap='coolwarm', edgecolor='none')



nclusters=75
#cluster_list=list(range(nclusters))
fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_free_energy_den(xHedges_bolas, yHedges_bolas, P_u, ax=ax)
#ax.scatter(dist_data_center[:,0],dist_data_center[:,1], s=9, c='white', marker="o")
#for i, txt in enumerate(cluster_list):
#    ax.annotate(txt, (dist_data_center[i,0], dist_data_center[i,1]),fontsize=2)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1)) 
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS_WHAM.cluster_center.012225.png",dpi=1200)



def plot_free_energy_f_bolas(
        xedge, yedge, f_bolas, weights=None, ax=None, nbins=100, ncontours=100,
        offset=-1, avoid_zero_count=False, minener_zero=True, kT=1.0,
        vmin=None, vmax=None, cmap='nipy_spectral', cbar=True,
        cbar_label='free energy / kT', cax=None, levels=None,
        legacy=True, cbar_orientation='vertical',norm=None,
        **kwargs):
    #f = _to_free_energy(z, minener_zero=minener_zero) * kT
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    mappable = ax.contourf(
        x, y, f_bolas.T, ncontours, norm=norm,
        vmin=vmin, vmax=vmax, cmap=cmap)
    misc = dict(mappable=mappable)
    if cbar_orientation not in ('horizontal', 'vertical'):
        raise ValueError(
            'cbar_orientation must be "horizontal" or "vertical"')
    if cbar:
        if cax is None:
            cbar_ = fig.colorbar(
                mappable, ax=ax, orientation=cbar_orientation)
        else:
            cbar_ = fig.colorbar(
                mappable, cax=cax, orientation=cbar_orientation)
        if cbar_label is not None:
            cbar_.set_label(cbar_label)
        misc.update(cbar=cbar_)
    if legacy:
        return fig, ax
    return fig, ax, misc





fig, ax = plt.subplots(figsize=(7, 6))
fig, ax = plot_free_energy_f_bolas(xHedges_bolas, yHedges_bolas, f_bolas_2500, ax=ax)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS_WHAM.test_2500.012224.png",dpi=1200)




x_window=np.array([[9,13],[10,16],[12.5,17.5],[16,20],[18,24],[19,26],[22.5,27.5],[25,30],[28,31],[29,33],[31.5,33.5],[33,35],[35,37],[36,39]]) 
y_window=np.array([[28,30],[26.5,28.5],[26,28],[25,27],[24,26],[22,25],[20,24],[20,24],[19,21],[17,20],[17,18],[16,17],[14.5,16],[13.5,15]]) 



def calculate_residuals(constant, X0,X1):
    # Your calculation logic here, using the constant and data
    # Return the residuals
    residuals = np.linalg.norm((X0+constant) - X1)
    return residuals


weights= np.ones(len(f_bolas_list))
p_weights=weights/weights.sum()







#initial_guess = 0.0
#result = minimize(calculate_residuals, initial_guess, args=(your_data,))
#optimized_constant = result.x[0]


constant_list=[0]
X_list=[]
Y_list=[]
f_bolas_list=[f_bolas_2500,f_bolas_5000,f_bolas_7500,f_bolas_10000,f_bolas_12500,f_bolas_15000,f_bolas_17500,f_bolas_20000,f_bolas_22500,f_bolas_25000,f_bolas_27500,f_bolas_30000,f_bolas_35000,f_bolas_37500]
for i in range(len(x_window)-1):
    if x_window[i+1][0]<x_window[i][1] and y_window[i][0]<y_window[i+1][1]:
        X,Y = np.mgrid[x_window[i+1][0]:x_window[i][1]:(x_window[i][1]-x_window[i+1][0])/10, y_window[i][0]:y_window[i+1][1]:(y_window[i+1][1]-y_window[i][0])/10]
        X_list.append(X)
        Y_list.append(Y)
        X_bin=np.digitize(X/10,xHedges_bolas)
        Y_bin=np.digitize(Y/10,yHedges_bolas)
        f_bolas_data=f_bolas_list[i]+constant_list[i]
        f_bolas_dataN=f_bolas_list[i+1]
        f_ene_data=np.empty((10,10))
        f_ene_dataN=np.empty((10,10))
        for j in range(len(X)):
            for k in range(len(Y)):
                f_ene_data[j,k]=f_bolas_data[X_bin[j,k],Y_bin[j,k]]
                f_ene_dataN[j,k]=f_bolas_dataN[X_bin[j,k],Y_bin[j,k]]
        initial_guess = 0.0
        result = minimize(calculate_residuals, initial_guess, args=(f_ene_data,f_ene_dataN))
        optimized_constant = result.x[0]
        constant_list.append(optimized_constant)
        print(i)
        print(optimized_constant)
        print(f_ene_data-f_ene_dataN)
        print("\n\n")
    elif x_window[i+1][0]<x_window[i][1] and y_window[i][0]==y_window[i+1][1]:
        X= np.linspace(x_window[i+1][0],x_window[i][1],10)
        X_list.append(X)
        Y=y_window[i][0]
        Y_list.append(y_window[i][0])
        X_bin=np.digitize(X/10,xHedges_bolas)
        Y_bin=np.digitize(Y/10,yHedges_bolas)
        f_bolas_data=f_bolas_list[i]+constant_list[i]
        f_bolas_dataN=f_bolas_list[i+1]
        f_ene_data=np.empty((10,1))
        f_ene_dataN=np.empty((10,1))
        for j in range(len(X)):
            f_ene_data[j]=f_bolas_data[X_bin[j],Y_bin]
            f_ene_dataN[j]=f_bolas_dataN[X_bin[j],Y_bin]
        initial_guess = 0.0
        result = minimize(calculate_residuals, initial_guess, args=(f_ene_data,f_ene_dataN))
        optimized_constant = result.x[0]
        constant_list.append(optimized_constant)
        print(i)
        print(optimized_constant)
        print(f_ene_data-f_ene_dataN)
        print("\n\n")
    elif x_window[i+1][0]==x_window[i][1] and y_window[i][0]<y_window[i+1][1]:
        Y = np.linspace(y_window[i][0],y_window[i+1][1],10)
        X=x_window[i+1][0]
        X_list.append(x_window[i+1][0])
        Y_list.append(Y)
        X_bin=np.digitize(X/10,xHedges_bolas)
        Y_bin=np.digitize(Y/10,yHedges_bolas)
        f_bolas_data=f_bolas_list[i]+constant_list[i]
        f_bolas_dataN=f_bolas_list[i+1]
        f_ene_data=np.empty((1,10))
        f_ene_dataN=np.empty((1,10))
        for k in range(len(Y)):
            f_ene_data[0][k]=f_bolas_data[X_bin,Y_bin[k]]
            f_ene_dataN[0][k]=f_bolas_dataN[X_bin,Y_bin[k]]
        initial_guess = 0.0
        result = minimize(calculate_residuals, initial_guess, args=(f_ene_data,f_ene_dataN))
        optimized_constant = result.x[0]
        constant_list.append(optimized_constant)
        print(i)
        print(optimized_constant)
        print(f_ene_data-f_ene_dataN)
        print("\n\n")
    elif x_window[i+1][0]==x_window[i][1] and y_window[i][0]==y_window[i+1][1]:
        X_list.append(x_window[i+1][0])
        Y_list.append(y_window[i][0])
        X=x_window[i+1][0]
        Y=y_window[i][0]
        X_bin=np.digitize(X/10,xHedges_bolas)
        Y_bin=np.digitize(Y/10,yHedges_bolas)
        f_bolas_data=f_bolas_list[i]
        f_bolas_dataN=f_bolas_list[i+1]
        f_ene_data=f_bolas_data[X_bin,Y_bin]+constant_list[i]
        f_ene_dataN=f_bolas_dataN[X_bin,Y_bin]
        initial_guess = 0.0
        result = minimize(calculate_residuals, initial_guess, args=(f_ene_data,f_ene_dataN))
        optimized_constant = result.x[0]
        constant_list.append(optimized_constant)
        print(i)
        print(optimized_constant)
        print(f_ene_data-f_ene_dataN)
        print("\n\n")
    else:
        constant_list.append(' ') 



f_bolas_merge=np.array(np.ones((100,100)) * np.inf)

for i in range(len(constant_list)):
    f_bolas_data=f_bolas_list[i]
    for m in range(len(f_bolas_2500)):
        for n in range(len(f_bolas_2500)):
            if not math.isinf(f_bolas_data[m][n]):
                f_bolas_merge[m][n]=f_bolas_data[m][n]+constant_list[i]
            
       

def plot_map(
        x, y, z, ax=None, cmap=None,
        ncontours=100, vmin=None, vmax=None, levels=None,
        cbar=True, cax=None, cbar_label=None,
        cbar_orientation='vertical', norm=None):
    """Plot a two-dimensional map from data on a grid.

    Parameters
    ----------
    x : ndarray(T)
        Binned x-coordinates.
    y : ndarray(T)
        Binned y-coordinates.
    z : ndarray(T)
        Binned z-coordinates.
    ax : matplotlib.Axes object, optional, default=None
        The ax to plot to; if ax=None, a new ax (and fig) is created.
    cmap : matplotlib colormap, optional, default=None
        The color map to use.
    ncontours : int, optional, default=100
        Number of contour levels.
    vmin : float, optional, default=None
        Lowest z-value to be plotted.
    vmax : float, optional, default=None
        Highest z-value to be plotted.
    levels : iterable of float, optional, default=None
        Contour levels to plot.
    cbar : boolean, optional, default=True
        Plot a color bar.
    cax : matplotlib.Axes object, optional, default=None
        Plot the colorbar into a custom axes object instead of
        stealing space from ax.
    cbar_label : str, optional, default=None
        Colorbar label string; use None to suppress it.
    cbar_orientation : str, optional, default='vertical'
        Colorbar orientation; choose 'vertical' or 'horizontal'.
    norm : matplotlib norm, optional, default=None
        Use a norm when coloring the contour plot.

    Optional parameters for contourf (**kwargs)
    -------------------------------------------
    corner_mask : boolean, optional
        Enable/disable corner masking, which only has an effect if
        z is a masked array. If False, any quad touching a masked
        point is masked out. If True, only the triangular corners
        of quads nearest those points are always masked out, other
        triangular corners comprising three unmasked points are
        contoured as usual.
        Defaults to rcParams['contour.corner_mask'], which
        defaults to True.
    alpha : float
        The alpha blending value.
    locator : [ None | ticker.Locator subclass ]
        If locator is None, the default MaxNLocator is used. The
        locator is used to determine the contour levels if they are
        not given explicitly via the levels argument.
    extend : [ ‘neither’ | ‘both’ | ‘min’ | ‘max’ ]
        Unless this is ‘neither’, contour levels are automatically
        added to one or both ends of the range so that all data are
        included. These added ranges are then mapped to the special
        colormap values which default to the ends of the
        colormap range, but can be set via
        matplotlib.colors.Colormap.set_under() and
        matplotlib.colors.Colormap.set_over() methods.
    xunits, yunits : [ None | registered units ]
        Override axis units by specifying an instance of a
        matplotlib.units.ConversionInterface.
    antialiased : boolean, optional
        Enable antialiasing, overriding the defaults. For filled
        contours, the default is True. For line contours, it is
        taken from rcParams[‘lines.antialiased’].
    nchunk : [ 0 | integer ]
        If 0, no subdivision of the domain. Specify a positive
        integer to divide the domain into subdomains of nchunk by
        nchunk quads. Chunking reduces the maximum length of polygons
        generated by the contouring algorithm which reduces the
        rendering workload passed on to the backend and also requires
        slightly less RAM. It can however introduce rendering
        artifacts at chunk boundaries depending on the backend, the
        antialiased flag and value of alpha.
    hatches :
        A list of cross hatch patterns to use on the filled areas.
        If None, no hatching will be added to the contour. Hatching
        is supported in the PostScript, PDF, SVG and Agg backends
        only.
    zorder : float
        Set the zorder for the artist. Artists with lower zorder
        values are drawn first.

    Returns
    -------
    fig : matplotlib.Figure object
        The figure in which the used ax resides.
    ax : matplotlib.Axes object
        The ax in which the map was plotted.
    misc : dict
        Contains a matplotlib.contour.QuadContourSet 'mappable' and,
        if requested, a matplotlib.Colorbar object 'cbar'.

    """
    import matplotlib.pyplot as _plt
    if ax is None:
        fig, ax = _plt.subplots()
    else:
        fig = ax.get_figure()
    mappable = ax.contourf(
        x, y, z, ncontours, norm=norm,
        vmin=vmin, vmax=vmax, cmap=cmap)
    misc = dict(mappable=mappable)
    if cbar_orientation not in ('horizontal', 'vertical'):
        raise ValueError(
            'cbar_orientation must be "horizontal" or "vertical"')
    if cbar:
        if cax is None:
            cbar_ = fig.colorbar(
                mappable, ax=ax, orientation=cbar_orientation)
        else:
            cbar_ = fig.colorbar(
                mappable, cax=cax, orientation=cbar_orientation)
        if cbar_label is not None:
            cbar_.set_label(cbar_label)
        misc.update(cbar=cbar_)
    return fig, ax, misc


def plot_free_energy_edit(
        x, y, f, weights=None, ax=None, nbins=100, ncontours=100,
        offset=-1, avoid_zero_count=False, minener_zero=True, kT=1.0,
        vmin=None, vmax=None, cmap='nipy_spectral', cbar=True,
        cbar_label='free energy / kT', cax=None, levels=None,
        legacy=True, ncountours=None, cbar_orientation='vertical',
        **kwargs):
    """Plot a two-dimensional free energy map using a histogram of
    scattered data.

    Parameters
    ----------
    xall : ndarray(T)
        Sample x-coordinates.
    yall : ndarray(T)
        Sample y-coordinates.
    weights : ndarray(T), optional, default=None
        Sample weights; by default all samples have the same weight.
    ax : matplotlib.Axes object, optional, default=None
        The ax to plot to; if ax=None, a new ax (and fig) is created.
        Number of contour levels.
    nbins : int, optional, default=100
        Number of histogram bins used in each dimension.
    ncontours : int, optional, default=100
        Number of contour levels.
    offset : float, optional, default=-1
        Deprecated and ineffective; raises a ValueError
        outside legacy mode.
    avoid_zero_count : bool, optional, default=False
        Avoid zero counts by lifting all histogram elements to the
        minimum value before computing the free energy. If False,
        zero histogram counts would yield infinity in the free energy.
    minener_zero : boolean, optional, default=True
        Shifts the energy minimum to zero.
    kT : float, optional, default=1.0
        The value of kT in the desired energy unit. By default,
        energies are computed in kT (setting 1.0). If you want to
        measure the energy in kJ/mol at 298 K, use kT=2.479 and
        change the cbar_label accordingly.
    vmin : float, optional, default=None
        Lowest free energy value to be plotted.
        (default=0.0 in legacy mode)
    vmax : float, optional, default=None
        Highest free energy value to be plotted.
    cmap : matplotlib colormap, optional, default='nipy_spectral'
        The color map to use.
    cbar : boolean, optional, default=True
        Plot a color bar.
    cbar_label : str, optional, default='free energy / kT'
        Colorbar label string; use None to suppress it.
    cax : matplotlib.Axes object, optional, default=None
        Plot the colorbar into a custom axes object instead of
        stealing space from ax.
    levels : iterable of float, optional, default=None
        Contour levels to plot.
    legacy : boolean, optional, default=True
        Switch to use the function in legacy mode (deprecated).
    ncountours : int, optional, default=None
        Legacy parameter (typo) for number of contour levels.
    cbar_orientation : str, optional, default='vertical'
        Colorbar orientation; choose 'vertical' or 'horizontal'.

    Optional parameters for contourf (**kwargs)
    -------------------------------------------
    corner_mask : boolean, optional
        Enable/disable corner masking, which only has an effect if
        z is a masked array. If False, any quad touching a masked
        point is masked out. If True, only the triangular corners
        of quads nearest those points are always masked out, other
        triangular corners comprising three unmasked points are
        contoured as usual.
        Defaults to rcParams['contour.corner_mask'], which
        defaults to True.
    alpha : float
        The alpha blending value.
    locator : [ None | ticker.Locator subclass ]
        If locator is None, the default MaxNLocator is used. The
        locator is used to determine the contour levels if they are
        not given explicitly via the levels argument.
    extend : [ ‘neither’ | ‘both’ | ‘min’ | ‘max’ ]
        Unless this is ‘neither’, contour levels are automatically
        added to one or both ends of the range so that all data are
        included. These added ranges are then mapped to the special
        colormap values which default to the ends of the
        colormap range, but can be set via
        matplotlib.colors.Colormap.set_under() and
        matplotlib.colors.Colormap.set_over() methods.
    xunits, yunits : [ None | registered units ]
        Override axis units by specifying an instance of a
        matplotlib.units.ConversionInterface.
    antialiased : boolean, optional
        Enable antialiasing, overriding the defaults. For filled
        contours, the default is True. For line contours, it is
        taken from rcParams[‘lines.antialiased’].
    nchunk : [ 0 | integer ]
        If 0, no subdivision of the domain. Specify a positive
        integer to divide the domain into subdomains of nchunk by
        nchunk quads. Chunking reduces the maximum length of polygons
        generated by the contouring algorithm which reduces the
        rendering workload passed on to the backend and also requires
        slightly less RAM. It can however introduce rendering
        artifacts at chunk boundaries depending on the backend, the
        antialiased flag and value of alpha.
    hatches :
        A list of cross hatch patterns to use on the filled areas.
        If None, no hatching will be added to the contour. Hatching
        is supported in the PostScript, PDF, SVG and Agg backends
        only.
    zorder : float
        Set the zorder for the artist. Artists with lower zorder
        values are drawn first.

    Returns
    -------
    fig : matplotlib.Figure object
        The figure in which the used ax resides.
    ax : matplotlib.Axes object
        The ax in which the map was plotted.
    misc : dict
        Contains a matplotlib.contour.QuadContourSet 'mappable' and,
        if requested, a matplotlib.Colorbar object 'cbar'.

    """
    if legacy:
        warnings.warn(
            'Legacy mode is deprecated is will be removed in the'
            ' next major release. Until then use legacy=False',
            DeprecationWarning)
        cmap = plt.get_cmap(cmap)
        if offset != -1:
            warnings.warn(
                'Parameter offset is deprecated and will be ignored',
                DeprecationWarning)
        if ncountours is not None:
            warnings.warn(
                'Parameter ncountours is deprecated;'
                ' use ncontours instead',
                DeprecationWarning)
            ncontours = ncountours
        if vmin is None:
            vmin = 0.0
    else:
        if offset != -1:
            raise ValueError(
                'Parameter offset is not allowed outside legacy mode')
        if ncountours is not None:
            raise ValueError(
                'Parameter ncountours is not allowed outside'
                ' legacy mode; use ncontours instead')
    #x, y, z = get_histogram(
    #    xall, yall, nbins=nbins, weights=weights,
    #    avoid_zero_count=avoid_zero_count)
    #f = _to_free_energy(z, minener_zero=minener_zero) * kT
    fig, ax, misc = plot_map(
        x, y, f, ax=ax, cmap=cmap,
        ncontours=ncontours, vmin=vmin, vmax=vmax, levels=levels,
        cbar=cbar, cax=cax, cbar_label=cbar_label,
        cbar_orientation=cbar_orientation, norm=None,
        **kwargs)
    if legacy:
        return fig, ax
    return fig, ax, misc




def plot_map_edit(
        x, y, z, weights=None, ax=None, nbins=100, ncontours=100,
        offset=-1, avoid_zero_count=False, minener_zero=True, kT=1.0,
        vmin=None, vmax=None, cmap='nipy_spectral', cbar=True,
        cbar_label='free energy / kT', cax=None, levels=None,norm=None,
        legacy=True, ncountours=None, cbar_orientation='vertical',
        **kwargs      ):

    if legacy:
        warnings.warn(
            'Legacy mode is deprecated is will be removed in the'
            ' next major release. Until then use legacy=False',
            DeprecationWarning)
        cmap = plt.get_cmap(cmap)
        if offset != -1:
            warnings.warn(
                'Parameter offset is deprecated and will be ignored',
                DeprecationWarning)
        if ncountours is not None:
            warnings.warn(
                'Parameter ncountours is deprecated;'
                ' use ncontours instead',
                DeprecationWarning)
            ncontours = ncountours
        if vmin is None:
            vmin = 0.0
    else:
        if offset != -1:
            raise ValueError(
                'Parameter offset is not allowed outside legacy mode')
        if ncountours is not None:
            raise ValueError(
                'Parameter ncountours is not allowed outside'
                ' legacy mode; use ncontours instead')    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    mappable = ax.contourf(
        x, y, z, ncontours, norm=norm,
        vmin=vmin, vmax=vmax, cmap=cmap)
    misc = dict(mappable=mappable)
    if cbar_orientation not in ('horizontal', 'vertical'):
        raise ValueError(
            'cbar_orientation must be "horizontal" or "vertical"')
    if cbar:
        if cax is None:
            cbar_ = fig.colorbar(
                mappable, ax=ax, orientation=cbar_orientation)
        else:
            cbar_ = fig.colorbar(
                mappable, cax=cax, orientation=cbar_orientation)
        if cbar_label is not None:
            cbar_.set_label(cbar_label)
        misc.update(cbar=cbar_)
    return fig, ax, misc




x = 0.5 * (xedge_bolas[:-1] + xedge_bolas[1:])
y = 0.5 * (yedge_bolas[:-1] + yedge_bolas[1:])


fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_map(x, y, f_bolas_merge,ax=ax)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.merge_constant.011824.png",dpi=1200)


fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_map_edit(x, y, f_bolas_merge)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.merge_constant.EDIT.011824.png",dpi=1200)




flist=[2500,5000,7500,10000,12500,15000,17500,20000,22500,25000,27500,30000,35000,37500]
i=0
for data in f_bolas_list:
    fig, ax = plt.subplots(figsize=(7, 6))
    fig, ax, misc = plot_map(x, y, data+constant_list[i],ax=ax)
    ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
    ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 4)
    ax.set_title('BOLAS', fontweight='bold',fontsize=24)
    ax.tick_params(labelsize=20)
    ax.set_xticks(np.arange(0,6,1))
    ax.set_yticks(np.arange(0,4,1))
    fig.tight_layout()
    fig.savefig("DIST_FEL_density.BOLAS.f_bolas_"+str(flist[i])+".011824.png",dpi=1200)
    i+=1



fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_map(x, y, f_bolas_2500,ax=ax)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.f_bolas_2500D.011824.png",dpi=1200)




fig, ax = plt.subplots(figsize=(7, 6))
fig, ax, misc = plot_map_edit(x, y, f_bolas_2500)
ax.set_xlabel('DIST 1  /nm',fontsize=20, fontweight='bold')
ax.set_ylabel('DIST 2  /nm',fontsize=20, fontweight='bold')
ax.set_xlim(0, 6)
ax.set_ylim(0, 4)
ax.set_title('BOLAS', fontweight='bold',fontsize=24)
ax.tick_params(labelsize=20)
ax.set_xticks(np.arange(0,6,1))
ax.set_yticks(np.arange(0,4,1))
fig.tight_layout()
fig.savefig("DIST_FEL_density.BOLAS.f_bolas_2500D.EDIT.011824.png",dpi=1200)





len(dist_data_2500_concatenated)                                                                                         
#Out[20]: 450000

len(dist_data_5000_concatenated)                                                                                         
#Out[21]: 332000

len(dist_data_7500_concatenated)                                                                                         
#Out[22]: 464000

len(dist_data_10000_concatenated)                                                                                        
#Out[23]: 267000

len(dist_data_12500_concatenated)                                                                                        
#Out[24]: 380000

len(dist_data_15000_concatenated)                                                                                        
#Out[25]: 189000

len(dist_data_17500_concatenated)                                                                                        
#Out[26]: 52000

len(dist_data_20000_concatenated)                                                                                        
#Out[27]: 433036

len(dist_data_22500_concatenated)                                                                                        
#Out[28]: 420000

len(dist_data_25000_concatenated)                                                                                        
#Out[29]: 2000

len(dist_data_27500_concatenated)                                                                                        
#Out[30]: 19000

len(dist_data_30000_concatenated)                                                                                        
#Out[31]: 73000

len(dist_data_35000_concatenated)                                                                                        
#Out[32]: 387000

len(dist_data_37500_concatenated)                                                                                        
#Out[33]: 309000

[2989755,2989756]
In [62]: start_25000                                                                                                                                                                  
Out[62]: 2987036

In [63]: start_27500                                                                                                                                                                  
Out[63]: 2989757


#start_7500=len(dist_data_5000_concatenated)+len(dist_data_2500_concatenated)
#start_10000=start_7500+len(dist_data_7500_concatenated)
#start_12500=start_10000+len(dist_data_10000_concatenated)
#start_15000=start_12500+len(dist_data_12500_concatenated)
#start_17500=start_15000+len(dist_data_15000_concatenated)
#start_20000=start_17500+len(dist_data_17500_concatenated)
#start_22500=start_20000+len(dist_data_20000_concatenated)
#start_25000=start_22500+len(dist_data_22500_concatenated)
#start_27500=start_25000+len(dist_data_25000_concatenated)
#start_30000=start_27500+len(dist_data_27500_concatenated)
#start_35000=start_30000+len(dist_data_30000_concatenated)
#start_37500=start_35000+len(dist_data_35000_concatenated)


#binnumber_2500 = binnumber_bolas.T[0:len(dist_data_2500_concatenated),]
#binnumber_5000 = binnumber_bolas.T[len(dist_data_2500_concatenated):start_7500,]
#binnumber_7500 = binnumber_bolas.T[start_7500:start_10000,]
#binnumber_10000 = binnumber_bolas.T[start_10000:start_12500,]
#binnumber_12500 = binnumber_bolas.T[start_12500:start_15000,]
#binnumber_15000 = binnumber_bolas.T[start_15000:start_17500,]
#binnumber_17500 = binnumber_bolas.T[start_17500:start_20000,]
#binnumber_20000 = binnumber_bolas.T[start_20000:start_22500,]
#binnumber_22500 = binnumber_bolas.T[start_22500:start_25000,]
#binnumber_25000 = binnumber_bolas.T[start_25000:(start_27500-2),]
#binnumber_27500 = binnumber_bolas.T[(start_27500-2):(start_30000-2),]
#binnumber_30000 = binnumber_bolas.T[(start_30000-2):(start_35000-2),]
#binnumber_35000 = binnumber_bolas.T[(start_35000-2):(start_37500-2),]
#binnumber_37500 = binnumber_bolas.T[(start_37500-2):,]


edges=[xHedges_bolas,yHedges_bolas]

def histogramdd_wedge(sample, edges, irange=None, density=None, weights=None):
    """
    Compute the multidimensional histogram of some data.

    Parameters
    ----------
    sample : (N, D) array, or (N, D) array_like
        The data to be histogrammed.

        Note the unusual interpretation of sample when an array_like:

        * When an array, each row is a coordinate in a D-dimensional space -
          such as ``histogramdd(np.array([p1, p2, p3]))``.
        * When an array_like, each element is the list of values for single
          coordinate - such as ``histogramdd((X, Y, Z))``.

        The first form should be preferred.

    bins : sequence or int, optional
        The bin specification:

        * A sequence of arrays describing the monotonically increasing bin
          edges along each dimension.
        * The number of bins for each dimension (nx, ny, ... =bins)
        * The number of bins for all dimensions (nx=ny=...=bins).

    range : sequence, optional
        A sequence of length D, each an optional (lower, upper) tuple giving
        the outer bin edges to be used if the edges are not given explicitly in
        `bins`.
        An entry of None in the sequence results in the minimum and maximum
        values being used for the corresponding dimension.
        The default, None, is equivalent to passing a tuple of D None values.
    density : bool, optional
        If False, the default, returns the number of samples in each bin.
        If True, returns the probability *density* function at the bin,
        ``bin_count / sample_count / bin_volume``.
    weights : (N,) array_like, optional
        An array of values `w_i` weighing each sample `(x_i, y_i, z_i, ...)`.
        Weights are normalized to 1 if density is True. If density is False,
        the values of the returned histogram are equal to the sum of the
        weights belonging to the samples falling into each bin.

    Returns
    -------
    H : ndarray
        The multidimensional histogram of sample x. See density and weights
        for the different possible semantics.
    edges : list
        A list of D arrays describing the bin edges for each dimension.

    See Also
    --------
    histogram: 1-D histogram
    histogram2d: 2-D histogram

    Examples
    --------
    >>> r = np.random.randn(100,3)
    >>> H, edges = np.histogramdd(r, bins = (5, 8, 4))
    >>> H.shape, edges[0].size, edges[1].size, edges[2].size
    ((5, 8, 4), 6, 9, 5)

    """

    try:
        # Sample is an ND-array.
        N, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = np.atleast_2d(sample).T
        N, D = sample.shape

    nbin = np.empty(D, np.intp)

    # normalize the range argument
    if irange is None:
        irange = (None,) * D
    elif len(irange) != D:
        raise ValueError('range argument must have one entry per dimension')

    # Create edge arrays
    for i in range(D):
        nbin[i] = len(edges[i]) + 1  # includes an outlier on each end
    #    dedges[i] = np.diff(edges[i])

    # Compute the bin number each sample falls into.
    Ncount = tuple(
        # avoid np.digitize to work around gh-11022
        np.searchsorted(edges[i], sample[:, i], side='right')
        for i in range(D)
    )

    # Using digitize, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right edge to be
    # counted in the last bin, and not as an outlier.
    for i in range(D):
        # Find which points are on the rightmost edge.
        on_edge = (sample[:, i] == edges[i][-1])
        # Shift these points one bin to the left.
        Ncount[i][on_edge] -= 1

    # Compute the sample indices in the flattened histogram matrix.
    # This raises an error if the array is too large.
    xy = np.ravel_multi_index(Ncount, nbin)

    # Compute the number of repetitions in xy and assign it to the
    # flattened histmat.
    hist = np.bincount(xy, weights, minlength=nbin.prod())

    # Shape into a proper matrix
    hist = hist.reshape(nbin)

    # This preserves the (bad) behavior observed in gh-7845, for now.
    hist = hist.astype(float, casting='safe')

    # Remove outliers (indices 0 and -1 for each dimension).
    core = D*(slice(1, -1),)
    hist = hist[core]

    if density:
        # calculate the probability density function
        s = hist.sum()
        for i in range(D):
            shape = np.ones(D, int)
            shape[i] = nbin[i] - 2
            hist = hist / dedges[i].reshape(shape)
        hist /= s

    if (hist.shape != nbin - 2).any():
        raise RuntimeError(
            "Internal Shape Error")
    return hist


z_2500 = histogramdd_wedge(dist_data_2500_concatenated,edges)
f_2500 = _to_free_energy(z_2500, minener_zero=minener_zero) * kT

z_5000 = histogramdd_wedge(dist_data_5000_concatenated,edges)
f_5000 = _to_free_energy(z_5000, minener_zero=minener_zero) * kT

z_7500 = histogramdd_wedge(dist_data_7500_concatenated,edges)
f_7500 = _to_free_energy(z_7500, minener_zero=minener_zero) * kT

z_10000 = histogramdd_wedge(dist_data_10000_concatenated,edges)
f_10000 = _to_free_energy(z_10000, minener_zero=minener_zero) * kT

z_12500 = histogramdd_wedge(dist_data_12500_concatenated,edges)
f_12500 = _to_free_energy(z_12500, minener_zero=minener_zero) * kT

z_15000 = histogramdd_wedge(dist_data_15000_concatenated,edges)
f_15000 = _to_free_energy(z_15000, minener_zero=minener_zero) * kT

z_17500 = histogramdd_wedge(dist_data_17500_concatenated,edges)
f_17500 = _to_free_energy(z_17500, minener_zero=minener_zero) * kT

z_20000 = histogramdd_wedge(dist_data_20000_concatenated,edges)
f_20000 = _to_free_energy(z_20000, minener_zero=minener_zero) * kT

z_22500 = histogramdd_wedge(dist_data_22500_concatenated,edges)
f_22500 = _to_free_energy(z_22500, minener_zero=minener_zero) * kT

z_25000 = histogramdd_wedge(dist_data_25000_concatenated,edges)
f_25000 = _to_free_energy(z_25000, minener_zero=minener_zero) * kT

z_27500 = histogramdd_wedge(dist_data_27500_concatenated,edges)
f_27500 = _to_free_energy(z_27500, minener_zero=minener_zero) * kT

z_30000 = histogramdd_wedge(dist_data_30000_concatenated,edges)
f_30000 = _to_free_energy(z_30000, minener_zero=minener_zero) * kT

z_35000 = histogramdd_wedge(dist_data_35000_concatenated,edges)
f_35000 = _to_free_energy(z_35000, minener_zero=minener_zero) * kT

z_37500 = histogramdd_wedge(dist_data_37500_concatenated,edges)
f_37500 = _to_free_energy(z_37500, minener_zero=minener_zero) * kT















