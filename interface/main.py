import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
from subprocess import call
import sys
rain=np.genfromtxt("../example/rain.dat")*1000*1000
hyperparam=[0.0373,2,0,44.78,0.112,0.01,0.01]
cor_len=[5]*8
lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,0,0])
upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
design=e.design("../example/design_outputs.dat","../example/design_pars.dat",64)
names=["Impervious area","Width","Slope","stor. imp.","$n_{imp}$","stor. per.",
       "% imp. area w/o dep. stor.","$n_{con}$","$\sigma^2_e$","$\sigma^2_b$"]
design.pick_first(64)
emu=e.emu(design,rain,hyperparam,cor_len,e_ini=3000,art="kalm",gamma=1.0,variance='false')
emu.condition()
errscale=0.01
eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=errscale)
eli.names=names

emu.emulate(design.test_pars[5])
emu.plot(design,["emu","swmm"])

