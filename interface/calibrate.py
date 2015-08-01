import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
design=e.design("design_data_full.dat","design_pars_full.dat",32)
measurement_path="measurement.dat"
design.pick_first(32)
pars_physical=[44.78,0.112,0.01,0.01,1000]
cor_len=[1.5]*10
rain=np.genfromtxt("rain_4_emu.dat")
names=["Impervious area","Width","Slope","$n_{imp}$","storage imp.","storage per.","% of imp. area w/o dep. sto.","$n_{con}$","Tue","Zue","$\sigma^2_e$","$\sigma^2_b$"]
# hyperparam=[0.0000231,0.000231,2000000,0]
hyperparam=[0.0000231,2000000,0]


lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1,0,0])
upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
emu=e.emu(design,rain,pars_physical,hyperparam,cor_len)

emu.condition()
# eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.01)

# swmm=s.swmm()
# eli.improve_emulator_for_lnlik(swmm,16)

# import emcee
# eli.sampler_pars(24,2000)
# sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
#                                 eli.lnprob,threads=8)
# sampler.run_mcmc(eli.pos, eli.length)
# eli.print_info_write_chain(sampler)
# sampler.pool.close()

# eli.chainz(sampler)


# cf.compare_two_posteriors(lower_par,upper_par,names,"samples_emulator_500_0.01.dat",
                          # "samples_emulator_501_0.01.dat")


%timeit emu.emulate(design.test_pars[5])
emu.plot(design,["emu","swmm"])

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



# emu.condition()
# elikelihood=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.1)
# # MCMC, (package emcee)

# import emcee
# elikelihood.sampler_pars(24,2002)
# sampler = emcee.EnsembleSampler(elikelihood.walkers, elikelihood.ndim,
#                                 elikelihood.lnprob,threads=8)
# sampler.run_mcmc(elikelihood.pos, elikelihood.length)
# elikelihood.print_info_write_chain(sampler)
# sampler.pool.close()

emu.condition()
param=np.genfromtxt("max_posterior_emulator_500_0.01.dat")
emu.emulate(param[0:10])
eli.better_loglikelihood(param)

swmm=s.swmm()

slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=0.01)
slikelihood.better_loglikelihood(param)


