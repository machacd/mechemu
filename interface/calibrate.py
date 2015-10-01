import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
from mpi4py import MPI
from subprocess import call
import emcee
import sys
from emcee.utils import MPIPool
measurement_path="measurement.dat"
pars_physical=[44.78,0.112,0.01,0.01,1000]
rain=np.genfromtxt("rain_4_emu.dat")
# hyperparam=[0.0000231,0.000231,2000000,0]
hyperparam=[0.00000973,2050407,0]
# cor_len=[1]*2
# lower_par=np.array([0.5,0.5,0,0])
# upper_par=np.array([1.1,1.5,1,1])
# design=e.design("design_data_twopar.dat","design_pars_twopar.dat",32)
# names=["Impervious area","Width","$\sigma^2_e$","$\sigma^2_b$"]
# design.pick_first(32)
# cor_len=[2]*4
# lower_par=np.array([0.5,0.5,0.5,0.5,0,0])
# upper_par=np.array([1.1,1.5,1.5,1.5,1,1])
# design=e.design("design_data_four.dat","design_pars_four.dat",48)
# names=["Impervious area","Width","Slope","stor. imp.","$\sigma^2_e$","$\sigma^2_b$"]
# design.pick_first(48)
cor_len=[5]*8
lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,0,0])
upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
design=e.design("design_data_full.dat","design_pars_full.dat",64)
names=["Impervious area","Width","Slope","stor. imp.","$n_{imp}$","stor. per.",
       "% imp. area w/o dep. stor.","$n_{con}$","$\sigma^2_e$","$\sigma^2_b$"]
design.pick_first(64)
emu=e.emu(design,rain,pars_physical,hyperparam,cor_len,e_ini=3000,art="kalm",gamma=1.5)
emu.condition()
errscale=0.01
eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=errscale)
eli.names=names
# swmm=s.swmm()
# emu.emulate(design.test_pars[2])

pool = MPIPool()
if not pool.is_master():
        pool.wait()
        sys.exit(0)
call(["rm","candidates_all.dat"])
swmm=s.swmm()
no_chains=32
base_length=540
new_len=base_length
top_length=base_length*no_chains
lower_length=int(np.floor(base_length*no_chains/3*2))
# eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=errscale)
# eli.add_candidates(swmm)
for k in np.arange(0):
    eli.sampler_pars(no_chains,new_len)
    sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                    eli.lnprob,pool=pool)
    sampler.run_mcmc(eli.pos, eli.length)
    eli.print_info_write_chain(sampler)
    filename=cf.create_file_name(["samples","emulator",str(new_len),str(errscale),".dat"])
    impsamples=np.genfromtxt(filename)[new_len+np.arange(16)*new_len-1][:,0:-2]
    means=np.mean(impsamples,0)
    np.savetxt("candidates.dat",(impsamples-means)*1.9+means,fmt="%10.5f")
    eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=errscale)
    eli.add_candidates(swmm)
    call(["./append_candidates.sh"])
    # cf.send_notification("iterace "+str(k)+" dokoncena na siam13 pro 4 param")
    new_len+=1
eli.sampler_pars(no_chains,544)
sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                eli.lnprob,pool=pool)
sampler.run_mcmc(eli.pos, eli.length)
eli.print_info_write_chain(sampler)
cf.send_notification("iterative run pro 8 param dokoncen na si15")
pool.close()

sys.exit(0)

# emu.plot(design,["emu","swmm","measurement"],swmm=swmm,likelihood=eli)
#                          "samples_swmm_508_0.01.dat"])

# cf.plot_posteriors(eli,["samples_emulator_480_0.01.dat","samples_swmm_500_0.01_two.dat"])


sampleNames=[cf.create_file_name(["samples","emulator",str(i),
                                  str(errscale),".dat"]) for i in np.arange(570,575)]
sampleNames.append("samples_swmm_1008_0.01.dat")
cf.plot_posteriors_iter(eli,sampleNames)

cf.plot_posteriors(eli,["samples_emulator_902_0.01.dat","samples_emulator_903_0.01.dat",
                        "samples_swmm_1009_0.01.dat"])



################################################################################
# production plots
cf.plot_posteriors(eli,["samples_emulator_1040_0.01.dat","samples_emulator_534_0.01.dat",
                        "samples_swmm_1014_0.01.dat"])

cf.plot_posteriors(eli,["samples_emulator_1008_0.01.dat","samples_emulator_1019_0.01.dat",
                        "samples_swmm_1008_0.01.dat"])

cf.plot_posteriors(eli,["samples_emulator_1002_0.01.dat",
                        "samples_emulator_1023_0.01.dat","samples_swmm_1002_0.01.dat"])

cf.plot_posteriors(eli,["samples_emulator_808_0.01.dat","samples_emulator_820_0.01.dat",
                        "samples_swmm_1008_0.01.dat"])

swmm=s.swmm()
pars_emu=eli.maxpo("samples_emulator_1019_0.01.dat")
pars_swmm=eli.maxpo("samples_emulator_1008_0.01.dat")
emu.plot_for_paper(design,pars_swmm,pars_emu,swmm,eli)


swmm=s.swmm()
samples=np.genfromtxt("samples_emulator_820_0.01.dat")
no_of_samples=320
no_of_threads=8
pars=np.zeros((no_of_samples,samples.shape[1]-2))
for i in np.arange(no_of_samples):
    for j in np.arange(samples.shape[1]-2):
        pars[i,j]=samples[np.random.randint(samples.shape[0]),j]
def swmm_multi(index):
    swmm.run(pars[index,:])
    return swmm.result
samples_swmm_output=np.zeros((no_of_samples,swmm.t))
from multiprocessing import Pool
for i in np.arange(0,no_of_samples,no_of_threads):
    p=Pool(no_of_threads)
    output_parallel_session=p.map(swmm_multi,np.arange(i,i+no_of_threads).tolist())
    samples_swmm_output[i:i+no_of_threads,:]=np.array(output_parallel_session)
    p.close()

samples_swmm_output_both=[]
samples_swmm_output_both.append(samples_swmm_output)

samples_swmm_output_both.append(samples_swmm_output)

emu.plot_for_paper_cumulative(samples_swmm_output_both,eli)



################################################################################


# cf.plot_posteriors(eli,["samples_emulator_1040_0.01.dat","samples_emulator_800_0.01.dat",
#                         "samples_swmm_1014_0.01.dat"])

# cf.plot_posteriors(eli,["samples_emulator_380_0.01.dat"
#                         ,"samples_emulator_383_0.01.dat", "samples_emulator_389_0.01.dat",
#                         "samples_swmm_508_0.01.dat",])


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

# emu.emulate(design.test_pars[17])
# emu.plot(design,["emu","swmm","measurement"],likelihood=eli,swmm=swmm)



# generate lliks

# no_of_sets=24
# no_of_threads=8
# new_len=1004
# # samples_idx=new_len+np.random.randint(no_chains,size=no_of_sets)*new_len-1-np.random.randint(np.floor(new_len*0.3),size=no_of_sets)
# samples_idx=new_len+np.arange(24)*new_len-1
# samples=np.genfromtxt("samples_swmm_1004_0.01.dat")[samples_idx]
# # samples=[np.hstack((design.test_pars[i],0.5,0.5)) for i in np.arange(64)]

# swmm=s.swmm()
# slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=errscale)
# def swmm_multi(index):
#     return slikelihood.lnprob(samples[index]) 
# def emu_multi(index):
#     return eli.lnprob(samples[index]) 
# lliks_s=[0]*no_of_sets
# from multiprocessing import Pool
# for i in np.arange(0,no_of_sets,no_of_threads):
#     p=Pool(no_of_threads)
#     output_parallel_session=p.map(swmm_multi,np.arange(i,i+no_of_threads).tolist())
#     lliks_s[i:i+no_of_threads]=np.array(output_parallel_session)
#     p.close()
# lliks_range=np.max(lliks_s)-np.min(lliks_s)

# lliks_e=[0]*no_of_sets
# for i in np.arange(0,no_of_sets,no_of_threads):
#     p=Pool(no_of_threads)
#     output_parallel_session=p.map(emu_multi,np.arange(i,i+no_of_threads).tolist())
#     lliks_e[i:i+no_of_threads]=np.array(output_parallel_session)
#     p.close()
# # # hist of loglikelihood values for the emulator
# hist,bins=np.histogram(lliks_s,50)
# width = 0.91 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# import matplotlib.pyplot as plt
# plt.clf()
# plt.bar(center, hist, align='center', width=width)
# plt.savefig("hist_e.pdf")
# difference=(np.array(lliks_s)-np.array(lliks_e))/lliks_range*100
# print(np.mean(abs(difference)))
# hist,bins=np.histogram(difference,40)
# width = 0.91 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# import matplotlib.pyplot as plt
# plt.clf()
# plt.bar(center, hist, align='center', width=width)
# plt.savefig("hist_diff.pdf")



# # # do some plotting

# swmm=s.swmm()
# maxpo=np.genfromtxt("./max_posterior_swmm_504_0.01.dat")
# emu.emulate(maxpo[0:4])
# emu.plot(design,["emu","swmm","measurement","bias"],swmm=swmm,likelihood=eli,full_pars=maxpo)

# param=np.hstack((design.test_pars[10],[0.4,0.6]))

# eli.lnprob(maxpo)

# slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=errscale)
# slikelihood.lnprob(maxpo)


# eli.lnprob(maxpo)

# slikelihood.lnprob(maxpo)
# param=np.genfromtxt("../data/pars_calibrated.dat")


# model=swmm.run_swmm(param[0:no_pars],"4","own")
# observations=measurement

# def bias_mean(param,model,observation,no_pars):
#     sigb=cf.inverse_box_cox(param[-1],0.35)
#     sige=cf.inverse_box_cox(param[-2],0.35)
#     covariance=cov_mat_b(sigb,10)\
#         +cov_mat_e(sige)
#     result=np.dot(cov_mat_b(sigb,10),
#                np.linalg.solve(covariance,observations-model))
#     return result
# bmean=bias_mean(param,model,observations,no_pars)
 
# cf.plot_triangle(lower_par,upper_par,names,"samples_swmm_1008_0.01.dat")

