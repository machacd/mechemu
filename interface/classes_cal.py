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
            self.ini_dd=result_producing_thing.dp.shape[0]
            self.added_counter=0
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
        self.names=["Impervious area","Width","Slope","$n_{imp}$","storage imp.","storage per.","% of imp. area w/o dep. sto.","$n_{con}$","Tue","Zue","$\sigma^2_e$","$\sigma^2_b$"]

    def better_loglikelihood(self,param_e):
        if self.result_producing_thing.typ=="emulator":
            self.result_producing_thing.emulate(param_e[0:-2])
        if self.result_producing_thing.typ=="swmm":
            self.result_producing_thing.run(param_e[0:-2])
        data=stats.boxcox((self.measurement>0)*self.measurement+0.01,0.35)
        mean=stats.boxcox((self.result_producing_thing.result>0)*self.result_producing_thing.result+0.01,0.35)
        covariance=param_e[-1]*self.cov_mat_b_base+\
            self.cov_mat_e_base*param_e[-2]
        lik=-0.5*np.linalg.slogdet(covariance)[1]-\
            0.5*np.dot(mean-data,np.linalg.solve(covariance,mean-data))-\
            0.5*self.t*np.log(2*np.pi)
        return lik

    def logprior(self,param_e):
        prior=0
        down=self.lower_bounds
        up=self.upper_bounds
        if np.all(down < param_e) and np.all(param_e < up):
            for i in np.arange(self.lower_bounds.shape[0]):
                prior+=np.log(self.prior_dist(param_e[i],i))
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
            print(lp,ll)
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

    def improve_emulator_for_lnlik(self,swmm,improvement_steps,design):
        import scipy.optimize as opt
        # from multiprocessing import Pool
        j=0
        # threads=8
        while j<improvement_steps:
            ret=opt.differential_evolution(self.lnprob,
                                           bounds=list(zip(self.lower_bounds,
                                                           self.upper_bounds))
                                           ,args=[False]
                                           ,disp=True, popsize=20,maxiter=20,
                                           polish=False)
            pars=ret.x[0:8]
            if cf.closest_distance(self.result_producing_thing.dp,pars)<0.5:
                pars=design.pars_all[self.ini_dd+self.added_counter]
                swmm.result=design.data_all[self.ini_dd+self.added_counter]
                self.added_counter+=1
            else:
                swmm.run(pars)
                with open("candidates.dat", 'ab') as file:
                    np.savetxt(file,pars,fmt='%10.5f', newline=' ')
                with open("candidates.dat", 'a') as file:
                    file.writelines("\n")
            self.result_producing_thing.dd=np.vstack((self.result_producing_thing.dd,swmm.result))
            self.result_producing_thing.dp=np.vstack((self.result_producing_thing.dp,
                                                pars))
            self.result_producing_thing.condition()
            j+=1

    def prior_dist(self,par,i):
        down=self.lower_bounds[i]
        up=self.upper_bounds[i]
        if (i<self.lower_bounds.shape[0]-2):
            mu=(up+down)/2
            sc=(up-down)/7
            prior=stats.truncnorm.pdf(par,(down-mu)/sc,
                                           (up-mu)/sc,
                                           loc=mu,scale=sc)
        else:
            prior=stats.gamma.pdf(par,1,scale=self.errscale)
        return prior


#service functions

    def print_info_write_chain(self,sampler):
        import datetime
        print("Mean acceptance fraction: {0:.3f}"
                              .format(np.mean(sampler.acceptance_fraction)))
        print(["Finished at: ",datetime.datetime.now()])
        samples = sampler.chain[:,:, :].reshape((-1,self.ndim))
        filename=cf.create_file_name(["samples",self.result_producing_thing.typ,
                                     str(self.length),str(self.errscale),".dat"])
        np.savetxt(filename,samples)
        self.max_posterior = np.zeros(self.ndim)
        for i in np.arange(self.ndim):
            self.max_posterior[i]=stats.mode(samples[:,i])[0]
            print(stats.mode(samples[:,i])[1])
        filename=cf.create_file_name(["max_posterior",self.result_producing_thing.typ,
                                     str(self.length),str(self.errscale),".dat"])
        np.savetxt(filename,self.max_posterior,fmt='%10.5f')

#extra experiments

    def add_candidates(self,swmm):
        candidates=np.genfromtxt("candidates.dat")
        candidates=cf.remove_duplicate_rows(candidates)
        for i in np.arange(candidates.shape[0]):
            swmm.run(candidates[i])
            # add at the end
            self.result_producing_thing.dd=np.vstack((self.result_producing_thing.dd,swmm.result))
            self.result_producing_thing.dp=np.vstack((self.result_producing_thing.dp,
                                                      candidates[i]))
            # remove something from the beginning, so the total number is constant
            self.result_producing_thing.dd=np.delete(self.result_producing_thing.dd,0,0)
            self.result_producing_thing.dp=np.delete(self.result_producing_thing.dp,0,0)
        self.result_producing_thing.condition()

