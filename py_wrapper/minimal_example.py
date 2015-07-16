import numpy as np
import emulib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import emcee
import triangle
import common_functions as cf
import run_swwm_for_calib_output as swmm
from scipy import stats
from scipy import optimize
import datetime
from multiprocessing import Pool
#
# set the id of the measurement we are using, 4 and for 9 for calibration or 11 for validation
event_id="4"
#load design data and other stuff
design_data_name=("min_ex_design_data.dat")
design_data_all=np.genfromtxt(event_id.join(design_data_name))
design_pars_name=("min_ex_design_pars.dat")
design_pars_all=np.genfromtxt(event_id.join(design_pars_name))
rain_name=("min_ex_rain_input.dat")
pars_physical=np.genfromtxt("min_ex_physical_pars.dat")
# set parameters of the emulator
n_total=16
nt=16
t=rain.shape[0]
no_pars=4
# set parameters of the emulator, to be fiddled with
m=1
n_u=16
test_set=3
dim_obs=1
cor_len=2
v_ini=1000
e_ini=3000 #bug, this needs to be inferred from the design data
input_dim=1
lambda_dim=1
# set the percentile of included data judged by the relative gradients 
hyperparam=np.zeros(lambda_dim+2*input_dim+dim_obs-1)
# set pars of the linmod for each event, accuracy isn't terribly important here
    hyperparam[0]=0.0000931
    hyperparam[1]=2000000
    hyperparam[2]=10
# set layouts, sample design data (useful in case n_total>>n_u),
# load validation data
obs_layout=np.arange(0,t*dim_obs,dim_obs)
indices=cf.pick_indicies(n_total,n_u)
design_data=design_data_all[indices.astype(int),0:t*dim_obs]
design_pars=design_pars_all[indices.astype(int),:]
test_data=design_data_all[n_total:n_total+nt,0:t]
test_pars=design_pars_all[n_total:n_total+nt,:]
param=test_pars[test_set]

# run the emulator, Kalman filter version
# conditionoing should not take too long and has to be done only once
conditioned=emulib.condition_kalman(m,dim_obs,n_u,t,no_pars,cor_len,input_dim,
                                lambda_dim,hyperparam,design_data,design_pars,rain,
                                    pars_physical,v_ini,e_ini)

# evaluation of the emulator for a particular parameter vector
# should be done in no time, this function can evaluated many times
result=emulib.evaluate_kalman(conditioned[0],conditioned[1],conditioned[2],
                          conditioned[3],m,dim_obs,n_u,t,no_pars,
                          cor_len,input_dim,lambda_dim,hyperparam,
                          param,design_data,design_pars,rain,pars_physical,v_ini,e_ini)
mean_kalm=result[0,0,:,0]


# run the emulator, non-Kalman filter version
# conditioning will take very long and will take up a lot of memory
z_prime=emulib.condition_nonkalman(m,dim_obs,n_u,t,no_pars,cor_len,input_dim,
                          lambda_dim,hyperparam,design_data_xt,design_pars_xt,rain,
                          pars_physical,v_ini,e_ini)

# evaluation of the emulator for a particular parameter vector
# should be done in no time, this function can evaluated many times
mean_nonkalm=emulib.evaluate_nonkalman(z_prime,m,dim_obs,n_u,t,no_pars,
                          cor_len,input_dim,lambda_dim,hyperparam,
                          param,design_data_xt,design_pars_xt,rain,pars_physical,v_ini,e_ini)


