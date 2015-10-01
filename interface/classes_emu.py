import numpy as np
import common_functions as cf
import emulib
class design(object):
    def __init__(self,path_d,path_p,n_total):
        self.data_all=np.genfromtxt(path_d)
        self.pars_all=np.genfromtxt(path_p)
        self.n_total=n_total
        self.data=self.data_all[0,:]
        self.pars=self.pars_all[0,:]
        self.test_data=self.data_all[self.n_total:self.
                                            data_all.shape[0],:]
        self.test_pars=self.pars_all[self.n_total:
                                            self.pars_all.shape[0],:]
        self.n_test=self.test_pars.shape[0]
    
    def pick_random(self,n):
        if n>self.n_total:
            raise RuntimeError('You cannot pick more dd than is the tot #')
        indices=np.zeros(n)
        for i in np.arange(n):
            indices[i]=np.random.randint(0,self.n_total)
            j=0
            while j<i:
                if indices[j]!=indices[i]:
                    j+=1
                else:
                    indices[i]=np.random.randint(0,self.n_total)
                    j=0
        self.data=self.data_all[indices.astype(int),:]
        self.pars=self.pars_all[indices.astype(int),:]
        self.n=n

    def pick_first(self,n):
        if n>self.n_total:
            raise RuntimeError('You cannot pick more dd than is the tot #')
        self.data=self.data_all[0:n,:]
        self.pars=self.pars_all[0:n,:]
        self.n=n

    def extend(self,data,pars):
        self.data=np.vstack((self.data,data))
        self.pars=np.vstack((self.pars,pars))
        self.n+=1
    
    def remove_test_set(self,index):
        if index>self.n_test:
            raise RuntimeError('Index is too large')
        self.test_pars=np.delete(self.test_pars,index,0)
        self.test_data=np.delete(self.test_data,index,0)
        self.n_test-=1



class emu(object):
    def __init__(self,design,inp,other_pars,hyperparameters,cor_len,
                 m=1,d_obs=1,e_ini=1000,v_ini=1000,
                 art="kalm",input_dim=1,gamma=2):
        self.m=m
        self.d_obs=d_obs
        self.dp=design.pars
        self.dd=design.data
        self.inp=inp
        self.other_pars=other_pars
        self.e_ini=e_ini
        self.v_ini=v_ini
        self.art=art
        self.cor_len=cor_len
        self.gamma=gamma
        self.lambda_dim=m
        self.input_dim=input_dim
        self.hyperparameters=hyperparameters
        self.typ="emulator"

    def condition(self):
        if self.art=="kalm":
            self.conditioned=emulib.condition_kalman(self.m,
                                                self.d_obs,
                                                self.dd.shape[0],
                                                self.dd.shape[1],
                                                self.dp.shape[1],
                                                self.cor_len,
                                                self.gamma,
                                                self.input_dim,
                                                self.lambda_dim,
                                                self.hyperparameters,
                                                self.dd,
                                                self.dp,
                                                self.inp,
                                                self.other_pars,
                                                self.v_ini,
                                                self.e_ini)
        elif self.art=="nonkalm":
            self.conditioned=emulib.condition_nonkalman(self.m,
                                                self.d_obs,
                                                self.dd.shape[0],
                                                self.dd.shape[1],
                                                self.dp.shape[1],
                                                self.cor_len,
                                                self.gamma,
                                                self.input_dim,
                                                self.lambda_dim,
                                                self.hyperparameters,
                                                self.dd,
                                                self.dp,
                                                self.inp,
                                                self.other_pars,
                                                self.v_ini,
                                                self.e_ini)

    def emulate(self,pars):
        self.pars=pars
        if self.art=="kalm":
            emulated=emulib.evaluate_kalman(self.conditioned[0],
                                               self.conditioned[1],
                                               self.conditioned[2],
                                               self.conditioned[3],
                                               self.m,
                                               self.d_obs,
                                               self.dd.shape[0],
                                               self.dd.shape[1],
                                               self.dp.shape[1],
                                               self.cor_len,
                                               self.gamma,
                                               self.input_dim,
                                               self.lambda_dim,
                                               self.hyperparameters,
                                               self.pars,
                                               self.dd,
                                               self.dp,
                                               self.inp,
                                               self.other_pars,
                                               self.v_ini,
                                               self.e_ini)
            self.result=emulated[0,0,:,0]
        elif self.art=="nonkalm":
            self.result=emulib.evaluate_nonkalman(self.conditioned,
                                               self.m,
                                               self.d_obs,
                                               self.dd.shape[0],
                                               self.dd.shape[1],
                                               self.dp.shape[1],
                                               self.cor_len,
                                               self.gamma,
                                               self.input_dim,
                                               self.lambda_dim,
                                               self.hyperparameters,
                                               self.pars,
                                               self.dd,
                                               self.dp,
                                               self.inp,
                                               self.other_pars,
                                               self.v_ini,
                                               self.e_ini)
            
    def improve(self,design,improvement_steps):
        # from multiprocessing import Pool
        j=0
        # threads=8
        while j<improvement_steps:
            # need to add a condition to be already initiated
            totals=np.zeros((design.n_test,3))
            self.condition()
            # multicore processing, does not work, try PATHOS instead of MULTIP
            # def emulated(index):
            #     self.emulate(design.test_pars[index])
            #     return(cf.rmse(self.result,design.test_data[index]))
            # magazine=np.zeros(threads)
            # mag_idx=0
            # i=0
            # while i<design.n_test:
            #     if mag_idx<threads:
            #         magazine[mag_idx]=i
            #         mag_idx+=1
            #     else:
            #         p=Pool(threads)
            #         output_parallel_session=p.map(emulated,magazine.tolist())
            #         totals[magazine.tolist()]=output_parallel_session
            #         mag_idx=0
            #         p.close()
            #     i+=1
            i=0
            while i<design.n_test:
                self.emulate(design.test_pars[i])
                totals[i,0]=cf.rmse(self.result,design.test_data[i])
                totals[i,1]=cf.sharpness(self.result,design.test_data[i])
                i+=1
            max_0=max(totals[:,0])
            max_1=max(totals[:,1])
            print(max_0)
            print(max_1)
            totals[:,0]=totals[:,0]/max_0
            totals[:,1]=totals[:,1]/max_1
            totals[:,2]=totals[:,1]+totals[:,0]
            worst=np.argmax(totals[:,2])
            current0=np.mean(totals[:,0])
            current1=np.mean(totals[:,1])
            current2=np.mean(totals[:,2])
            print(max(totals[:,2]))
            print(current0)
            print(current1)
            print(current2)
            print("--------")
            j+=1
            design.extend(design.test_data[worst],design.test_pars[worst])
            design.remove_test_set(worst)
            self.dp=design.pars
            self.dd=design.data
 
    def objective_rmse(self,new_cor_len,design,n_test_sample):
        totals=np.zeros(n_test_sample)
        self.cor_len=new_cor_len
        self.condition()
        i=0
        while i<n_test_sample:
            self.emulate(design.test_pars[i])
            totals[i]=cf.rmse(self.result,design.test_data[i])
            i+=1
        current=np.mean(totals)
        # print(current)
        # print(self.cor_len)
        return(current)

    def objective_hyperpars(self,pars):
        self.cor_len=[pars[0]]*self.dp.shape[1]
        self.hyperparameters=np.hstack((pars[1:3],0))
        t=self.dd.shape[1]
        self.condition()
        suma=0
        for i in np.arange(t):
            covariance=self.conditioned[3][:,:,i,0]
            covariance_inv=self.conditioned[3][:,:,i,1]
            sds=np.zeros(16)+1
            covariance=np.diagflat(sds)
            covariance_inv=covariance
            mean=self.conditioned[1][:,i]
            swmm=self.dd[:,i]
            suma+=0.5*np.linalg.slogdet(covariance)[1]+\
            0.5*np.dot(mean-swmm,np.dot(covariance_inv,mean-swmm))+\
            0.5*t*np.log(2*np.pi)
        return suma

    def estimate_hyperpars(self):
        import scipy.optimize as opt
        lower=np.array([0,0,0])
        upper=np.array([0.001,0.01,100000000])
        ret=opt.differential_evolution(self.objective_hyperpars,
                                       bounds=list(zip(lower,
                                                       upper))
                                       ,disp=True, popsize=20,maxiter=20,
                                       polish=False)
        return ret




    def plot(self,design,what,likelihood=0,swmm=0,full_pars=0):
        #«what» can contain:
        #kalm, nkalm, meas, test, swmm
        import matplotlib as mpl
        mpl.use('Agg')
        mpl.rcParams.update({'font.size': 20})
        import matplotlib.pyplot as plt
        plt.clf()
        plt.figure(figsize=(16,16))
        t=self.dd.shape[1]/self.d_obs
        time=np.arange(0,t)
        for i in np.arange(0,self.d_obs):
            plt.subplot(self.d_obs,1,i+1)
            obs_layout=np.arange(i,t*self.d_obs,self.d_obs)
            # simulator=design.test_data[test_set,:]
            indices=np.arange(1,self.dd.shape[0])
            for i in np.nditer(indices):
                plt.plot(time,design.data[i,:],'0.8',linewidth=0.3),
            plt.plot(time,design.data[0,:],'0.8',linewidth=0.3,label="design data"),
            if "swmm" in what:
                found=0
                for i in range(design.n_test):
                    if (self.pars==design.test_pars[i,:]).all():
                        plt.plot(time,design.test_data[i],color="red",label="SWMM",lw=2)
                        found=1
                if found==0:
                    swmm.run(self.pars)
                    plt.plot(time,swmm.result,color="red",label="SWMM",lw=2)
            if "bias" in what:
                plt.plot(time,swmm.result+likelihood.bias_mean(full_pars,swmm.result),
                         color="cyan",label="SWMM with bias",lw=2)
            if "emu" in what:
                plt.plot(time,self.result,color="green",label="emulator",lw=2)
            if "measurement" in what:
                plt.plot(time,likelihood.measurement,color="blue",label="measurement",lw=2)
            # plt.plot(time,measurement[obs_layout],color="black",label="measurement")
            # plt.plot(time,test[obs_layout],color="magenta",label="swmm, calibrated")
            plt.ylabel('Q [l$\cdot$s$^{-1}$]')
            plt.xlabel('t [min]')
            plt.ylim([0,400])
            plt.grid()
        plt.legend(fancybox=True,loc=1)
        plt.savefig('result.pdf',dpi=500)
        plt.close()

    def plot_for_paper(self,design,pars1,pars2,swmm,likelihood):
        import matplotlib as mpl
        mpl.use('Agg')
        mpl.rcParams.update({'font.size': 20})
        import matplotlib.pyplot as plt
        plt.clf()
        plt.figure(figsize=(16,16))
        t=self.dd.shape[1]/self.d_obs
        time=np.arange(0,t)
        maxima=np.zeros(t)
        minima=np.zeros(t)
        for i in time:
            maxima[i]=np.max(design.data[:,i])
            minima[i]=np.min(design.data[:,i])
        for i in np.arange(0,self.d_obs):
            plt.subplot(self.d_obs,1,i+1)
            obs_layout=np.arange(i,t*self.d_obs,self.d_obs)
            swmm.run(pars1[0:-2])
            plt.plot(time,maxima,'0.6',linewidth=0.2),
            plt.plot(time,minima,'0.6',linewidth=0.2),
            plt.fill_between(time,maxima,minima,alpha=0.2,label="Calibration bounds",color="grey")
            plt.plot(time,swmm.result,color="red",label="SWMM, w. SWMM",lw=2)
            plt.plot(time,swmm.result+likelihood.bias_mean(pars1,swmm.result),
                     'r-.',label=" + bias",lw=3)
            swmm.run(pars2[0:-2])
            plt.plot(time,swmm.result,color="green",label="SWMM, w. emulator",lw=2)
            plt.plot(time,swmm.result+likelihood.bias_mean(pars2,swmm.result),
                     'g-.',label=" + bias",lw=3)
            plt.plot(time,likelihood.measurement,'b:',label="measurement",lw=2)
            plt.ylabel('Q [l$\cdot$s$^{-1}$]')
            plt.xlabel('t [min]')
            plt.ylim([0,400])
            plt.grid()
        plt.legend(fancybox=True,loc=1)
        plt.savefig('result.pdf',dpi=500)
        plt.close()

    def plot_for_paper_cumulative(self,samples,likelihood):
        import matplotlib as mpl
        mpl.use('Agg')
        mpl.rcParams.update({'font.size': 28})
        import matplotlib.pyplot as plt
        plt.clf()
        plt.figure(figsize=(24,16))
        t=samples[0].shape[1]
        time=np.arange(0,t)
        obs_layout=np.arange(t)
        labels=["SWMM posterior","emu posterior"]
        colors=["red","green"]
        for i in np.arange(2):
            plt.subplot(1,2,i+1)
            for j in np.arange(samples[i].shape[0]):
                plt.plot(time,samples[i][j,:],alpha=0.05,color=colors[i])
            plt.plot(time,[-1]*t,alpha=1,color=colors[i],label=labels[i])
            plt.plot(time,likelihood.measurement,alpha=0.50,color='black',label="measurement",lw=2)
            plt.ylabel('Q [l$\cdot$s$^{-1}$]')
            plt.xlabel('t [min]')
            plt.ylim([0,400])
            plt.grid()
            plt.legend(fancybox=True,loc=1)
        plt.tight_layout()
        plt.savefig('result_cumulative.pdf',dpi=500)
        plt.close()




class MyBounds(object):
    def __init__(self, xmax=[3.5]*10, xmin=[0.5]*10):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin



            







