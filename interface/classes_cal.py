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
    
    def cov_mat_b(self,sig_b,tau):
        result=np.zeros((self.t,self.t))
        exponentials=np.zeros(self.t)
        for i in np.arange(self.t):
            exponentials[i]=np.exp(-1/tau*i)
        for i in np.arange(self.t):
            for j in np.arange(0,i+1):
                result[i,j]=sig_b*exponentials[i-j]
        return cf.symmetrize(result)

    def cov_mat_e(self,sig_e):
        sds=np.zeros(self.t)+sig_e
        result=np.diagflat(sds)
        return result

    def better_loglikelihood(self,param_e):
        if self.result_producing_thing.typ=="emulator":
            self.result_producing_thing.emulate(param_e[0:10])
        if self.result_producing_thing.typ=="swmm":
            self.result_producing_thing.run(param_e[0:10])
        data=stats.boxcox(abs(100+self.measurement),0.35)
        mean=stats.boxcox(abs(100+self.result_producing_thing.result),0.35)
        sig2=param_e[11]
        covariance=self.cov_mat_b(sig2,10)+\
            self.cov_mat_e(param_e[10])
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



    def first_estimate(self,start_value,mybounds):
        import scipy.optimize as opt
        ret=opt.basinhopping(self.lnprob,start_value,niter=10000,T=100,stepsize=0.5,minimizer_kwargs={"method":"L-BFGS-B", "args":(False)},accept_test=mybounds)
        return ret

# these bounds are ony for the first_estimate
class MyBounds(object):
    def __init__(self, xmax=[0.5]*12, xmin=[1.5]*12):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin
