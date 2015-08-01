import numpy as np
import common_functions as cf
from scipy import stats
class likelihood(object):
    def __init__(self,measurement_path,result_producing_thing,lower_bounds,
                 upper_bounds,errscale=0.2):
        self.measurement=np.genfromtxt(measurement_path)
        self.upper_bounds=upper_bounds
        self.lower_bounds=lower_bounds
        self.errscale=errscale
        self.result_producing_thing=result_producing_thing
        if result_producing_thing.typ=="emulator":
            self.t=result_producing_thing.dd.shape[1]
        elif result_producing_thing.typ=="swmm":
            self.t=result_producing_thing.t
        
        cov_mat_b_base=np.zeros((self.t,self.t))
        exponentials=np.zeros(self.t)
        tau=10
        for i in np.arange(self.t):
            exponentials[i]=np.exp(-1/tau*i)
        for i in np.arange(self.t):
            for j in np.arange(0,i+1):
                cov_mat_b_base[i,j]=exponentials[i-j]
        self.cov_mat_b_base=cf.symmetrize(cov_mat_b_base)
        sds=np.zeros(self.t)+1
        self.cov_mat_e_base=np.diagflat(sds)

    def better_loglikelihood(self,param_e):
        if self.result_producing_thing.typ=="emulator":
            self.result_producing_thing.emulate(param_e[0:10])
        if self.result_producing_thing.typ=="swmm":
            self.result_producing_thing.run(param_e[0:10])
        data=stats.boxcox(abs(100+self.measurement),0.35)
        mean=stats.boxcox(abs(100+self.result_producing_thing.result),0.35)
        covariance=param_e[11]*self.cov_mat_b_base+\
            self.cov_mat_e_base*param_e[10]
        lik=-0.5*np.linalg.slogdet(covariance)[1]-\
            0.5*np.dot(mean-data,np.linalg.solve(covariance,mean-data))
        print(lik)
        return lik

    def logprior(self,param_e):
        prior=0
        down=self.lower_bounds
        up=self.upper_bounds
        if np.all(down < param_e) and np.all(param_e < up):
            for i in np.arange(down.shape[0]):
                if (i<10):
                    mu=(up[i]+down[i])/2
                    sc=(up[i]-down[i])/10
                    prior+=np.log(stats.truncnorm.pdf(param_e[i],(down[i]-mu)/sc,
                                                      (up[i]-mu)/sc,
                                                      loc=mu,scale=sc))
                else:
                    prior+=np.log(stats.gamma.pdf(param_e[i],1,scale=self.errscale))
            return prior
        return -np.inf

    def lnprob(self,param_e,maximize=True):
        lp=self.logprior(param_e)
        if not np.isfinite(lp):
            if maximize:    
                return -np.inf
            else:
                return np.inf
        ll=self.better_loglikelihood(param_e)
        if np.isnan(ll):
            if maximize:    
                return -np.inf
            else:
                return np.inf
        if maximize:    
            return lp + ll
        else:
            return -lp -ll

    def sampler_pars(self,nwalkers,length):
        self.walkers=nwalkers
        self.length=length
        self.ndim = self.upper_bounds.shape[0]
        self.pos = [self.lower_bounds+(self.upper_bounds-self.lower_bounds)*\
                    np.random.rand(self.ndim) for i in range(nwalkers)]
        self.extent=[(1,1)]*self.ndim
        for i in np.arange(self.ndim):
            self.extent[i]=(self.lower_bounds[i],self.upper_bounds[i])

    def improve_emulator_for_lnlik(self,swmm,improvement_steps):
        import scipy.optimize as opt
        # from multiprocessing import Pool
        j=0
        # threads=8
        while j<improvement_steps:
            ret=opt.differential_evolution(self.lnprob,
                                           bounds=list(zip(self.lower_bounds,
                                                           self.upper_bounds))
                                           ,args=[False]
                                           ,disp=True, popsize=10,maxiter=10,
                                           polish=False)
            swmm.run(ret.x[0:10])
            self.result_producing_thing.dd=np.vstack((self.result_producing_thing.dd,swmm.result))
            self.result_producing_thing.dp=np.vstack((self.result_producing_thing.dp,
                                                ret.x[0:10]))
            self.result_producing_thing.condition()
            with open("candidates.dat", 'ab') as file:
                np.savetxt(file,ret.x,fmt='%10.5f', newline=' ')
            with open("candidates.dat", 'a') as file:
                file.writelines("\n")
            j+=1


    def print_info_write_chain(self,sampler):
        import datetime
        print("Mean acceptance fraction: {0:.3f}"
                              .format(np.mean(sampler.acceptance_fraction)))
        print(["Finished at: ",datetime.datetime.now()])
        samples = sampler.chain[:,:, :].reshape((-1,self.ndim))
        filename="samples_"
        filename+=self.result_producing_thing.typ
        filename+="_"
        filename+=str(self.length)
        filename+="_"
        filename+=str(self.errscale)
        filename+=".dat"
        np.savetxt(filename,samples)
        self.max_posterior = np.zeros(self.ndim)
        for i in np.arange(self.ndim):
            self.max_posterior[i]=stats.mode(samples[:,i])[0]
            print(stats.mode(samples[:,i])[1])
        filename="max_posterior_"
        filename+=self.result_producing_thing.typ
        filename+="_"
        filename+=str(self.length)
        filename+="_"
        filename+=str(self.errscale)
        filename+=".dat"
        np.savetxt(filename,self.max_posterior,fmt='%10.5f')

    def chainz(self,sampler):
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        chain_no=np.random.randint(0,self.walkers)
        plt.clf()
        f,axes=plt.subplots(2,6,figsize=(24,8))
        row=0
        col=0
        for i in np.arange(self.ndim):
            axes[row,col].plot(sampler.chain[chain_no,:,i])
            row+=1
            if (row==2):
                row=0
                col+=1
        f.tight_layout()
        filename="chainz_"
        filename+=self.result_producing_thing.typ
        filename+="_"
        filename+=str(self.length)
        filename+="_"
        filename+=str(self.errscale)
        filename+=".pdf"
        f.savefig(filename)

