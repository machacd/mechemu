import numpy as np
def nash_sut(observed,modelled):
    time=observed.shape[0]
    average=sum(observed)/time
    denom=sum((observed-average)**2)
    numer=sum((observed-modelled)**2)
    return 1-numer/denom

def rmse(a,b):
    return np.sqrt(sum((a-b)**2)/a.shape[0])

def sharpness(a,b):
    af=np.fft.fft(a)[np.floor(a.shape[0])/6:np.floor(a.shape[0])/2]
    bf=np.fft.fft(b)[np.floor(b.shape[0])/6:np.floor(a.shape[0])/2]
    return(sum(np.abs(af-bf)))

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

def compare_two_posteriors(lower_par,upper_par,names,path1="samples_emu.dat",
                           path2="samples_swmm.dat"):
    from scipy import stats
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    down=lower_par
    up=upper_par
    emu_samples=np.genfromtxt(path1)
    swmm_samples=np.genfromtxt(path2)
    plt.clf()
    f,axes=plt.subplots(4,3,figsize=(16,16))
    row=0
    col=0
    for i in np.arange(lower_par.shape[0]):
        handle_hist=axes[row,col].hist([emu_samples[:,i]],
                                       bins=100,
                                       alpha=1,
                                       range=[down[i],up[i]],
                                       color='crimson',
                                       histtype='stepfilled',
                                       normed=True)
        handle_hist2=axes[row,col].hist([swmm_samples[:,i]],
                                        bins=100,
                                        alpha=0.73,
                                        range=[down[i],up[i]],
                                        color='black',
                                        histtype='stepfilled',
                                        normed=True)
        axes[row,col].set_title(names[i])
        if (i<10):
            mu=(up[i]+down[i])/2
            sc=(up[i]-down[i])/10
            x=np.arange(down[i],up[i],0.01)
            prior=stats.truncnorm.pdf(x,(down[i]-mu)/sc,
                                              (up[i]-mu)/sc,
                                              loc=mu,scale=sc)
            axes[row,col].plot(x,prior)
        else:
            x=np.arange(down[i],up[i],0.01)
            prior=stats.gamma.pdf(x,1,scale=0.01)
            axes[row,col].plot(x,prior)
        row+=1
        if (row==4):
            row=0
            col+=1
    f.tight_layout()
    f.legend((handle_hist[2]),("First",),fancybox=True,loc=4)
    f.legend((handle_hist2[2]),("Second",),fancybox=True,loc=4)
    f.savefig("posterior.pdf",dpi=500)
    plt.close()

def create_file_name(strings):
    strings=np.array(strings)
    filename=""
    for i in np.arange(strings.shape[0]):
        filename+=strings[i]
        if i!=strings.shape[0]-2 and i!=strings.shape[0]-1:
            filename+="_"
    return filename

def plot_all_chains(input_file,total_chains,no_pars):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    samples=np.genfromtxt(input_file)
    chain_length=samples.shape[0]/total_chains
    plt.clf()
    f,axes=plt.subplots(2,6,figsize=(24,8))
    for i in np.arange(total_chains):
        chain=samples[i*chain_length:(i+1)*chain_length,:]
        row=0
        col=0
        for j in np.arange(no_pars):
            axes[row,col].plot(chain[:,j],'b',alpha=0.54)
            row+=1
            if (row==2):
                row=0
                col+=1
    f.tight_layout()
    filename=create_file_name(["chainz",input_file,".pdf"])
    f.savefig(filename)

def plot_triangle(lower_par,upper_par,names,filename):
    import triangle
    inp_file=np.genfromtxt(filename)
    extent=[(1,1)]*lower_par.shape[0]
    for i in np.arange(lower_par.shape[0]):
        extent[i]=(lower_par[i],upper_par[i])
    fig = triangle.corner(inp_file,labels=names,exents=extent)
    filename=create_file_name(["triangle_",filename,".png"])
    fig.savefig(filename)




