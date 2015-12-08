import numpy as np
import sys
sys.path.append("../interface")
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
from subprocess import call
rain=np.genfromtxt("rain.dat")*1000*1000

#set hyperparameters
hyperparam=[0.0373,2,0,44.78,0.112,0.01,0.01]

#specify correlation lenght and the amount of design data sets
cor_len=[0.1]*8
design=e.design("design_outputs.dat","design_pars.dat",64)
design.pick_first(16)

#initialize emulator
emu=e.emu(design,rain,hyperparam,cor_len,e_ini=0,art="kalm",gamma=1.7,variance='false')

#condition
emu.condition()

# emulate, this can be the be inserted into a loop or whatever
emu.emulate(design.test_pars[7])

# plot results
emu.plot(design,["emu","swmm"])

