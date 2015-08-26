import numpy as np
import classes_emu as e
import classes_cal as c
import classes_swmm as s
import common_functions as cf
measurement_path="measurement.dat"
pars_physical=[44.78,0.112,0.01,0.01,1000]
rain=np.genfromtxt("rain_4_emu.dat")
# hyperparam=[0.0000231,0.000231,2000000,0]
hyperparam=[0.0000091,3000000,0]
cor_len=[10]*8
lower_par=np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,0,0])
upper_par=np.array([1.1,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1])
design=e.design("design_data_full.dat","design_pars_full.dat",128)
# cor_len=[10]*2
# lower_par=np.array([0.5,0.5,0,0])
# upper_par=np.array([1.1,1.5,5,5])
# design=e.design("design_data_twopar.dat","design_pars_twopar.dat",128)
design.pick_first(32)
emu=e.emu(design,rain,pars_physical,hyperparam,cor_len,e_ini=3000,art="kalm",gamma=1)
emu.condition()
errscale=0.01
# eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=errscale)



swmm=s.swmm()
emu.emulate(design.test_pars[1])
emu.plot(design,["emu","swmm","measurement"],swmm=swmm,likelihood=eli)



cf.plot_posteriors(eli,[cf.create_file_name(["samples","emulator",str(i),
                                             str(errscale),".dat"]) for i in np.arange(453,457)])

cf.plot_posteriors(eli,["samples_emulator_554_0.01.dat","samples_swmm_501_0.01.dat"])


swmm=s.swmm()
# eli.add_candidates(swmm)
# eli.improve_emulator_for_lnlik(swmm,48,design)
no_chains=24
base_length=453
new_len=base_length
top_length=base_length*no_chains
lower_length=int(np.floor(base_length*no_chains/3*2))
from subprocess import call
call(["rm","candidates_all.dat"])
import emcee
for k in np.arange(4):
    eli.sampler_pars(no_chains,new_len)
    sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                    eli.lnprob,threads=8)
    sampler.run_mcmc(eli.pos, eli.length)
    eli.print_info_write_chain(sampler)
    sampler.pool.close()
    filename=cf.create_file_name(["samples","emulator",str(new_len),str(errscale),".dat"])
    impsamples=np.genfromtxt(filename)[new_len+np.arange(4)*new_len-1][:,0:-2]
    np.savetxt("candidates.dat",impsamples,fmt="%10.5f")
    eli=c.likelihood(measurement_path,emu,lower_par,upper_par,errscale=errscale)
    eli.add_candidates(swmm)
    call(["./append_candidates.sh"])
    new_len+=1
    cf.send_notification("iterace "+str(k)+" dokoncena na giantocto")
new_len-=1
eli.sampler_pars(no_chains,480)
sampler = emcee.EnsembleSampler(eli.walkers, eli.ndim,
                                eli.lnprob,threads=8)
sampler.run_mcmc(eli.pos, eli.length)
eli.print_info_write_chain(sampler)
sampler.pool.close()
cf.send_notification("gintocto run pro 96 + zlepseni dokoncen, s dobrym modelem")

cf.plot_posteriors(eli,["samples_emulator_480_0.01.dat","samples_swmm_500_0.01_two.dat"])

# cf.plot_posteriors(eli,["samples_emulator_300_0.01.dat","samples_swmm_501_0.01.dat"])





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


param=np.hstack((design.test_pars[10],[0.4,0.6]))
eli.lnprob(param)

slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=errscale)
slikelihood.lnprob(param)

# emu.emulate(design.test_pars[17])
# emu.plot(design,["emu","swmm","measurement"],likelihood=eli,swmm=swmm)



# generate lliks

no_of_sets=128
no_of_threads=8
samples_idx=new_len+np.random.randint(no_chains,size=no_of_sets)*new_len-1-np.random.randint(np.floor(new_len*0.3),size=no_of_sets)
samples=np.genfromtxt("samples_swmm_501_0.01.dat")[samples_idx]

# samples=[np.hstack((design.test_pars[i],0.45,0.6)) for i in np.arange(128)]

swmm=s.swmm()
slikelihood=c.likelihood(measurement_path,swmm,lower_par,upper_par,errscale=errscale)
def swmm_multi(index):
    return slikelihood.lnprob(samples[index]) 
def emu_multi(index):
    return eli.lnprob(samples[index]) 
lliks_s=[0]*no_of_sets
from multiprocessing import Pool
for i in np.arange(0,no_of_sets,no_of_threads):
    p=Pool(no_of_threads)
    output_parallel_session=p.map(swmm_multi,np.arange(i,i+no_of_threads).tolist())
    lliks_s[i:i+no_of_threads]=np.array(output_parallel_session)
    p.close()

lliks_e=[0]*no_of_sets
for i in np.arange(0,no_of_sets,no_of_threads):
    p=Pool(no_of_threads)
    output_parallel_session=p.map(emu_multi,np.arange(i,i+no_of_threads).tolist())
    lliks_e[i:i+no_of_threads]=np.array(output_parallel_session)
    p.close()
# # hist of loglikelihood values for the emulator
hist,bins=np.histogram(lliks_s,50)
width = 0.91 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
import matplotlib.pyplot as plt
plt.clf()
plt.bar(center, hist, align='center', width=width)
plt.savefig("hist_e.pdf")
difference=np.array(lliks_s)-np.array(lliks_e)
print(np.mean(abs(difference)))
hist,bins=np.histogram(difference,20)
width = 0.91 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
import matplotlib.pyplot as plt
plt.clf()
plt.bar(center, hist, align='center', width=width)
plt.savefig("hist_diff.pdf")



# # do some plotting

# maxpo=np.genfromtxt("./max_posterior_emulator_1000_0.01.dat")
# emu.emulate(maxpo[0:8])
# emu.plot(design,["emu","swmm","measurement"],swmm=swmm,likelihood=eli)

# eli.lnprob(maxpo)

# slikelihood.lnprob(maxpo)
