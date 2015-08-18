import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
measurement_path="measurement.dat"
pars_physical=[44.78,0.112,0.01,0.01,1000]
rain=np.genfromtxt("rain_4_emu.dat")
# hyperparam=[0.0000231,0.000231,2000000,0]
hyperparam=[0.00000131,2000000,0]
cor_len=[10]*8
lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0])
upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
design=e.design("design_data_full.dat","design_pars_full.dat",128)
design.pick_first(16)
emu=e.emu(design,rain,pars_physical,hyperparam,cor_len,e_ini=3000,art="kalm",gamma=1)
emu.condition()
eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.5)

cf.plot_posteriors(eli,["samples_emulator_1000_0.005.dat","samples_emulator_1001_0.005.dat",
                        "samples_emulator_1002_0.005.dat","samples_emulator_2000_0.005.dat"])

# cf.plot_posteriors(eli,["samples_emulator_100_0.005.dat","samples_emulator_101_0.005.dat",
#                         "samples_emulator_102_0.005.dat"])
# %timeit emu.emulate(design.test_pars[10])
# emu.plot(design,["emu","swmm"],swmm=swmm,likelihood=eli)


swmm=s.swmm()
# eli.add_candidates(swmm)
# eli.improve_emulator_for_lnlik(swmm,48,design)

no_chains=24
base_length=1000
new_len=base_length
top_length=base_length*no_chains
lower_length=int(np.floor(base_length*no_chains/3*2))
from subprocess import call
import emcee

for k in np.arange(3):
    eli.sampler_pars(no_chains,new_len)
    sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                    eli.lnprob,threads=8)
    sampler.run_mcmc(eli.pos, eli.length)
    eli.print_info_write_chain(sampler)
    sampler.pool.close()
    filename=cf.create_file_name(["samples","emulator",str(new_len),"0.005",".dat"])
    samples=np.genfromtxt(filename)[np.random.randint(lower_length,top_length,size=16)][:,0:-2]
    np.savetxt("candidates.dat",samples,fmt="%10.5f")
    eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.005)
    eli.add_candidates(swmm)
    # call(["mv","candidates.dat","candidates_1.dat"])
    new_len+=1

eli.sampler_pars(no_chains,1000)
sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                eli.lnprob,threads=8)
sampler.run_mcmc(eli.pos, eli.length)
eli.print_info_write_chain(sampler)
sampler.pool.close()





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

# emu.condition()
# elikelihood=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=0.01)
# emu.condition()

# param=np.genfromtxt("max_posterior_emulator_2000_0.01.dat")
# emu.emulate(param[0:10])

swmm=s.swmm()

slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=0.005)

param=np.hstack((design.test_pars[10],[0.2,0.3]))
eli.lnprob(param)
slikelihood.lnprob(param)

# emu.emulate(design.test_pars[17])
# emu.plot(design,["emu","swmm","measurement"],likelihood=eli,swmm=swmm)



# slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=0.01)

samples=np.genfromtxt("samples_swmm_200_0.005.dat")[np.random.randint(2400,4800,size=200)]


# samples=samples[:,0:7]


lliks_s=slikelihood.log_liks_ddata(samples)

# samples=np.genfromtxt("lik_comp_samples.dat")
# lliks_s=np.genfromtxt("lik_swmm.dat")

# # hist of loglikelihood values for the emulator

swmm=s.swmm()
eli.add_candidates(swmm)

lliks_e=eli.log_liks_ddata(samples)
hist,bins=np.histogram(lliks_e,50)
width = 0.91 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
import matplotlib.pyplot as plt
plt.clf()
plt.bar(center, hist, align='center', width=width)
plt.savefig("hist_e.pdf")
difference=np.abs(np.array(lliks_s)-np.array(lliks_e))
hist,bins=np.histogram(difference,20)
width = 0.91 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
import matplotlib.pyplot as plt
plt.clf()
plt.bar(center, hist, align='center', width=width)
plt.savefig("hist_diff.pdf")





cf.compare_two_posteriors(eli,"samples_swmm_2000_0.005.dat",
                          "samples_swmm_2000_0.005.dat")

cf.compare_two_posteriors(eli,"samples_swmm_2001_0.005.dat",
                          "samples_swmm_2001_0.005.dat")
# cf.compare_two_posteriors(lower_par,upper_par,names,"samples_emulator_2001_0.01.dat",
#                           "samples_swmm_2000_0.01.dat")

cf.compare_two_posteriors(eli,"samples_emulator_600_0.005.dat",
                          "samples_swmm_200_0.005.dat")

# cf.compare_two_posteriors(lower_par,upper_par,names,"samples_emulator_500_0.005.dat",
#                           "samples_emulator_5000_0.005.dat")

# cf.plot_triangle(lower_par,upper_par,names,"samples_swmm_2000_0.01.dat")


