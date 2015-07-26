import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
design=e.design("design_data_full.dat","design_pars_full.dat",32)
measurement_path="measurement.dat"
design.pick_first(32)
hyperparam=[0.0000231,2000000,0]
pars_physical=[44.78,0.112,0.01,0.01,1000]
cor_len=[1.3]*10
rain=np.genfromtxt("rain_4_emu.dat")
names=["Impervious area","Width","Slope","$n_{imp}$","storage imp.","storage per.","% of imp. area w/o dep. sto.","$n_{con}$","Tue","Zue","$\sigma^2_e$","$\sigma^2_b$"]
hyperparam=[0.0000231,0.000231,2000000,0]
emu=e.emu(design,rain,pars_physical,hyperparam,cor_len,art="nonkalm",m=2)
emu.condition()
emu.emulate(design.test_pars[4])
emu.plot(design,["emu","swmm"])

# hyperparam=[0.0000231,0.000231,2000000,0]
# emu=e.emu(design,rain,pars_physical,hyperparam,cor_len,m=2)
# emu.condition()
# emu.emulate(design.test_pars[4])
# emu.plot(design,["emu","swmm"])

# import scipy.optimize as opt
# emu=e.emu(design,rain,pars_physical,hyperparam,cor_len)
# mybounds=e.MyBounds()
# ret=opt.basinhopping(emu.objective_rmse,[2.4]*10,minimizer_kwargs={"method":"L-BFGS-B", "args":(design,16)},accept_test=mybounds,disp=True)
# cor_len=ret.x
# emu=e.emu(design,rain,pars_physical,hyperparam,cor_len)

# emu.improve(design,16)


# # # calibration pars
# # lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1,0,0])
# # upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
lower_par=np.array([0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1,0,0])
upper_par=np.array([1.1,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1,1])
# # start_value=[0.7]*12
# # bounds_first_estimate=c.MyBounds(xmin=lower_par,xmax=upper_par)
# # cal=c.likelihood(measurement_path,emu,lower_par,upper_par)
# # ret=cal.first_estimate(start_value,bounds_first_estimate)
# # print(ret)


# # import matplotlib.pyplot as plt
# # t=np.arange(design.data.shape[1])
# # yf0=np.fft.fft(design.data[0])
# # yf1=np.fft.fft(design.data[1])
# # xf=np.linspace(0,1/2,t.shape[0]/2)
# # # freq=np.fft.fftfreq(t.shape[-1])
# # # freq=np.fft.fftshift(freq)
# # plt.clf()
# # plt.figure(figsize=(8,8))
# # plt.plot(xf,2/t.shape[0]*np.abs(yf0[0:t.shape[0]/2]),xf,2/t.shape[0]*np.abs(yf1[0:t.shape[0]/2]))
# # plt.grid()
# # plt.savefig('fft.pdf',dpi=500)
# # plt.close()


# emu.condition()

elikelihood=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.1)
# MCMC, (package emcee)
import emcee
elikelihood.sampler_pars(24,2002)
sampler = emcee.EnsembleSampler(elikelihood.walkers, elikelihood.ndim, elikelihood.lnprob)
                                # ,threads=1)
sampler.run_mcmc(elikelihood.pos, elikelihood.length)
elikelihood.print_info_write_chain(sampler)
sampler.pool.close()

# # param_w_errormod=np.genfromtxt("../data/pars_calibrated_swmm.dat")
# # emu.emulate(elikelihood.max_posterior[0:10])
# elikelihood.better_loglikelihood(elikelihood.max_posterior)

# swmm=s.swmm()
# slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par)
# slikelihood.better_loglikelihood(elikelihood.max_posterior)


