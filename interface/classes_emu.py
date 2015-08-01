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
                 art="kalm",input_dim=1):
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
        self.lambda_dim=m
        self.input_dim=input_dim
        self.hyperparameters=hyperparameters
        self.typ="emulator"
        self.degree=1
        self.distances=np.zeros(self.dd.shape[0])+1
    
    def create_distance_matrix(self,degree):
        n=self.dd.shape[0]
        self.distances=np.zeros(n)+1
        self.degree=degree
        for i in np.arange(n):
            minima=np.zeros((2,degree))+np.inf
            for j in np.arange(degree):
                for k in np.arange(n):
                    if not k in minima[0,:] and k!=i:
                        dist=np.linalg.norm(self.dp[i,:]-self.dp[k,:])
                        if dist<minima[1,j]:
                            minima[1,j]=dist
                            minima[0,j]=k
            self.distances[i]=np.mean(minima[1,:])


    def add_distance(self,pars):
        if np.all(self.distances==1):
            return 1
        else:
            n=self.dd.shape[0]
            minima=np.zeros((2,self.degree))+np.inf
            for j in np.arange(self.degree):
                for k in np.arange(n):
                    if not k in minima[0,:]:
                        dist=np.linalg.norm(pars-self.dp[j,:])
                        if dist<minima[1,j]:
                            minima[1,j]=dist
                            minima[0,j]=k
            return np.mean(minima[1,:])



    def condition(self):
        if self.art=="kalm":
            self.conditioned=emulib.condition_kalman(self.m,
                                                self.d_obs,
                                                self.dd.shape[0],
                                                self.dd.shape[1],
                                                self.dp.shape[1],
                                                self.cor_len,
                                                self.input_dim,
                                                self.lambda_dim,
                                                self.hyperparameters,
                                                self.dd,
                                                self.dp,
                                                self.inp,
                                                self.other_pars,
                                                self.v_ini,
                                                self.e_ini,
                                                self.distances)
        elif self.art=="nonkalm":
            self.conditioned=emulib.condition_nonkalman(self.m,
                                                self.d_obs,
                                                self.dd.shape[0],
                                                self.dd.shape[1],
                                                self.dp.shape[1],
                                                self.cor_len,
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
                                               self.input_dim,
                                               self.lambda_dim,
                                               self.hyperparameters,
                                               self.pars,
                                               self.dd,
                                               self.dp,
                                               self.inp,
                                               self.other_pars,
                                               self.v_ini,
                                               self.e_ini,
                                               self.distances,
                                               self.add_distance(self.pars))
            # print(self.add_distance(self.pars))
            self.result=emulated[0,0,:,0]
        elif self.art=="nonkalm":
            self.result=emulib.evaluate_nonkalman(self.conditioned,
                                               self.m,
                                               self.d_obs,
                                               self.dd.shape[0],
                                               self.dd.shape[1],
                                               self.dp.shape[1],
                                               self.cor_len,
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


    def plot(self,design,what):
        #«what» can contain:
        #kalm, nkalm, meas, test, swmm
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        plt.clf()
        plt.figure(figsize=(8,8))
        t=self.dd.shape[1]/self.d_obs
        time=np.arange(0,t)
        for i in np.arange(0,self.d_obs):
            plt.subplot(self.d_obs,1,i+1)
            obs_layout=np.arange(i,t*self.d_obs,self.d_obs)
            # simulator=design.test_data[test_set,:]
            indices=np.arange(1,self.dd.shape[0])
            for i in np.nditer(indices):
                plt.plot(time,design.data[i,:],'0.8',linewidth=0.1),
            plt.plot(time,design.data[0,:],'0.8',linewidth=0.1,label="design data"),
            if "swmm" in what:
                for i in range(design.n_test):
                    if (self.pars==design.test_pars[i,:]).all():
                        plt.plot(time,design.test_data[i],color="red",label="SWMM")
            # plt.plot(time,mean_nonkalm[obs_layout],color="blue",label="emulator nonkalm")
            if "emu" in what:
                plt.plot(time,self.result,color="green",label="emulator")
            # plt.plot(time,measurement[obs_layout],color="black",label="measurement")
            # plt.plot(time,test[obs_layout],color="magenta",label="swmm, calibrated")
            plt.ylabel('Q [l$\cdot$s$^{-1}$]')
            plt.xlabel('t [min]')
            plt.ylim([0,1000])
            plt.grid()
        plt.legend(fancybox=True,loc=1)
        plt.savefig('result.pdf',dpi=500)
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



            







