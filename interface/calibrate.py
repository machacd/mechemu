import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
design=e.design("design_data_full.dat","design_pars_full.dat",128)
measurement_path="measurement.dat"
design.pick_first(32)
pars_physical=[44.78,0.112,0.01,0.01,1000]
cor_len=[2.0]*10
rain=np.genfromtxt("rain_4_emu.dat")
names=["Impervious area","Width","Slope","$n_{imp}$","storage imp.","storage per.","% of imp. area w/o dep. sto.","$n_{con}$","Tue","Zue","$\sigma^2_e$","$\sigma^2_b$"]
# hyperparam=[0.0000231,0.000231,2000000,0]
hyperparam=[0.00000831,2000000,0]

lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1,0,0])
upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
emu=e.emu(design,rain,pars_physical,hyperparam,cor_len,e_ini=3000,art="kalm")
emu.condition()
# hist of loglikelihood values for the emulator
eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.01)

# lliks_e=eli.log_liks_ddata(design.test_pars[0:128])
# hist,bins=np.histogram(lliks_e,50)
# width = 0.91 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# import matplotlib.pyplot as plt
# plt.clf()
# plt.bar(center, hist, align='center', width=width)
# plt.savefig("hist_e.pdf")


# swmm=s.swmm()
# slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=0.01)
# lliks_s=slikelihood.log_liks_ddata(design.test_pars[0:128])

# difference=np.abs(np.array(lliks_s)-np.array(lliks_e))
# hist,bins=np.histogram(difference,80)
# width = 0.91 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# import matplotlib.pyplot as plt
# plt.clf()
# plt.bar(center, hist, align='center', width=width)
# plt.savefig("hist_diff.pdf")




# swmm=s.swmm()
# eli.improve_emulator_for_lnlik(swmm,16)

import emcee
eli.sampler_pars(24,2005)
sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                eli.lnprob,threads=8)
sampler.run_mcmc(eli.pos, eli.length)
eli.print_info_write_chain(sampler)
sampler.pool.close()

# eli.chainz(sampler)


# cf.compare_two_posteriors(lower_par,upper_par,names,"samples_swmm_2000_0.1_b.dat",
#                           "samples_emulator_2000_1.dat")

cf.compare_two_posteriors(lower_par,upper_par,names,"samples_emulator_2001_0.01.dat",
                          "samples_swmm_2000_0.01.dat")

cf.compare_two_posteriors(lower_par,upper_par,names,"samples_emulator_2000_0.01.dat",
                          "samples_emulator_2005_0.01.dat")

cf.plot_triangle(lower_par,upper_par,names,"samples_swmm_2000_0.01.dat")
# %timeit emu.emulate(design.test_pars[19])
# emu.plot(design,["emu","swmm"])

# hyperparam=[0.0000231,0.000231,2000000,0]
# emu=e.emu(design,rain,pars_physical,hyperparam,cor_len)
# emu.condition()
# emu.emulate(design.test_pars[4])
# emu.plot(design,["emu","swmm"])

# import scipy.optimize as opt
# emu=e.emu(design,rain,pars_physical,hyperparam,cor_len)
# mybounds=e.MyBounds()
# ret=opt.basinhopping(emu.objective_rmse,[2.4]*10,minimizer_kwargs={"method":"L-BFGS-B", "args":(design,16)},accept_test=mybounds,disp=True)
# cor_len=ret.x
# emu=e.emu(design,rain,pars_physical,hyperparam,cor_len)

# emu.improve(design,15)

# initial=(lower_par+upper_par])/2

# middle=(lower_par[0:10]+upper_par[0:10])/2
# swmm=s.swmm()
# swmm.run(middle)
# design.extend(swmm.result,middle)


# import scipy.optimize as opt
# ret=opt.differential_evolution(elikelihood.lnprob,bounds=list(zip(lower_par,upper_par)),args=[False],disp=True, popsize=2,maxiter=5)
# bounds=c.MyBounds()
# ret=opt.basinhopping(elikelihood.lnprob,initial,niter=100,T=10,stepsize=0.2,minimizer_kwargs={"method":"L-BFGS-B", "args":(False)},accept_test=bounds,disp=True,interval=1)

# # # calibration pars
# lower_par=np.array([0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1,0,0])
# upper_par=np.array([1.1,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1,1])
# start_value=[0.7]*12
# bounds_first_estimate=c.MyBounds(xmin=lower_par,xmax=upper_par)
# cal=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.1)
# ret=cal.first_estimate(start_value,bounds_first_estimate)
# print(ret)


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

########VALIDATION PART

emu.condition()
elikelihood=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.01)
emu.condition()
param=np.genfromtxt("max_posterior_swmm_2000_0.01.dat")
emu.emulate(param[0:10])
eli.better_loglikelihood(param)

# %timeit emu.emulate(design.test_pars[19])
emu.plot(design,["emu","swmm","measurement"],likelihood=eli,swmm=swmm)


swmm=s.swmm()

slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=0.01)
slikelihood.better_loglikelihood(param)
