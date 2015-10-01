import numpy as np
def nash_sut(observed,modelled):
    time=observed.shape[0]
    average=sum(observed)/time
    denom=sum((observed-average)**2)
    numer=sum((observed-modelled)**2)
    return 1-numer/denom

def inverse_box_cox(transformed,lambd):
    return(transformed*lambd+1)**(1/lambd)

def rmse(a,b):
    return np.sqrt(sum((a-b)**2)/a.shape[0])

def sharpness(a,b):
    af=np.fft.fft(a)[np.floor(a.shape[0])/6:np.floor(a.shape[0])/2]
    bf=np.fft.fft(b)[np.floor(b.shape[0])/6:np.floor(a.shape[0])/2]
    return(sum(np.abs(af-bf)))

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

def plot_posteriors(eli,file_names):
    from scipy import stats
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams.update({'font.size': 20})
    import matplotlib.pyplot as plt
    no_files=np.array(file_names).shape[0]
    samples=list()
    for i in np.arange(no_files):
        samples.append(np.genfromtxt(file_names[i]))
    down=eli.lower_bounds
    up=eli.upper_bounds
    no_of_pars=up.shape[0]
    plt.clf()
    f,axes=plt.subplots(2,int(np.ceil(no_of_pars/2)),figsize=(24,12))
    row=0
    col=0
    colors='bgrcmykwbgrcmykw'
    linetypes=['--','-','-','-','-','-','-']
    for i in np.arange(up.shape[0]):
        # patches=[]
        lines=[]
        x=np.arange(down[i],up[i],0.01)
        pri=np.zeros(x.shape[0])
        for j in np.arange(x.shape[0]):
            pri[j]=eli.prior_dist(x[j],i)
        leg1,=axes[row,col].plot(x,pri,'k.-')
        D1,p1=stats.ks_2samp(samples[no_files-1][:,i],samples[no_files-3][:,i])
        D2,p2=stats.ks_2samp(samples[no_files-1][:,i],samples[no_files-2][:,i])
        print(D1)
        print(D2)
        print("_________")
        for j in np.arange(no_files):
            kde=stats.gaussian_kde(samples[j][:,i])
            kde.covariance_factor = lambda : .3
            kde._compute_covariance()
            line,=axes[row,col].plot(x,kde.evaluate(x),linestyle=linetypes[j],color=colors[j],lw=2)
            lines.append(line)
        axes[row,col].set_title(eli.names[i])
        if i>up.shape[0]-3:
            axes[row,col].set_ylim([0,15])
        row+=1
        if (row==2):
            row=0
            col+=1
    f.tight_layout(rect=[0, 0.08, 1, 1])
    # f.legend([leg1,patches[2],patches[1],patches[0]],["Prior distribution","SWMM posterior","Emulator (improved) posterior","Emulator (standard) posterior",], bbox_to_anchor=[0.5, 0.05],loc='center',ncol=2)
    f.legend([leg1,lines[2],lines[1],lines[0]],["Prior distribution","SWMM posterior","Emulator (improved) posterior","Emulator (standard) posterior",], bbox_to_anchor=[0.5, 0.05],loc='center',ncol=2)
    f.savefig("posteriors.pdf",dpi=500)
    plt.close()

def plot_posteriors_iter(eli,file_names):
    from scipy import stats
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams.update({'font.size': 20})
    import matplotlib.pyplot as plt
    no_files=np.array(file_names).shape[0]
    samples=list()
    for i in np.arange(no_files):
        samples.append(np.genfromtxt(file_names[i]))
    down=eli.lower_bounds
    up=eli.upper_bounds
    no_of_pars=up.shape[0]
    plt.clf()
    f,axes=plt.subplots(2,int(np.ceil(no_of_pars/2)),figsize=(24,12))
    row=0
    col=0
    colors='bgcymrkwbgrcmykw'
    linetypes=['--','-.',':','--','-','-','-']
    kolmog_stats=np.zeros((2*no_files-1,up.shape[0]))
    for i in np.arange(up.shape[0]):
        # patches=[]
        lines=[]
        x=np.arange(down[i],up[i],0.01)
        for j in np.arange(no_files):
            Ds,ps=stats.ks_2samp(samples[j][:,i],samples[-1][:,i])
            kolmog_stats[j*2,i]=Ds
            if j>0:
                D,p=stats.ks_2samp(samples[j][:,i],samples[j-1][:,i])
                kolmog_stats[j*2-1,i]=D
                # print(D)
            kde=stats.gaussian_kde(samples[j][:,i])
            kde.covariance_factor = lambda : .3
            kde._compute_covariance()
            line,=axes[row,col].plot(x,kde.evaluate(x),linestyle=linetypes[j],color=colors[j],lw=3)
            lines.append(line)
        axes[row,col].set_title(eli.names[i])
        print("_________")
        if i>up.shape[0]-3:
            axes[row,col].set_ylim([0,15])
        row+=1
        if (row==2):
            row=0
            col+=1
    f.tight_layout(rect=[0, 0.13, 1, 1])
    np.savetxt("kolmog.dat",kolmog_stats,delimiter=" & ", fmt="%.2f")
    # f.legend([leg1,patches[2],patches[1],patches[0]],["Prior distribution","SWMM posterior","Emulator (improved) posterior","Emulator (standard) posterior",], bbox_to_anchor=[0.5, 0.05],loc='center',ncol=2)
    f.legend([lines[0],lines[1],lines[2],lines[3],lines[4],lines[5]],["0$^{th}$ iteration","1$^{st}$ iteration","2$^{nd}$ iteration","3$^{rd}$ iteration","4$^{th}$ iteration","SWMM",], bbox_to_anchor=[0.5, 0.08],loc='center',ncol=2)
    f.savefig("posteriors.pdf",dpi=500)
    plt.close()


def plot_posteriors_cumulative(eli,file_name,amount_of_samples):
    from scipy import stats
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams.update({'font.size': 20})
    samples=np.genfromtxt(file_name)
    down=eli.lower_bounds
    up=eli.upper_bounds
    no_of_pars=up.shape[0]
    plt.clf()
    f,axes=plt.subplots(2,int(np.ceil(no_of_pars/2)),figsize=(24,12))
    row=0
    col=0
    colors='bgrcmykwbgrcmykw'
    linetypes=['--','-','-','-','-','-','-']
    for i in np.arange(up.shape[0]):
        # patches=[]
        lines=[]
        x=np.arange(down[i],up[i],0.01)
        pri=np.zeros(x.shape[0])
        for j in np.arange(x.shape[0]):
            pri[j]=eli.prior_dist(x[j],i)
        leg1,=axes[row,col].plot(x,pri,'k.-')
        D1,p1=stats.ks_2samp(samples[no_files-1][:,i],samples[no_files-3][:,i])
        D2,p2=stats.ks_2samp(samples[no_files-1][:,i],samples[no_files-2][:,i])
        print(D1)
        print(D2)
        print("_________")
        for j in np.arange(no_files):
            kde=stats.gaussian_kde(samples[j][:,i])
            kde.covariance_factor = lambda : .3
            kde._compute_covariance()
            line,=axes[row,col].plot(x,kde.evaluate(x),linestyle=linetypes[j],color=colors[j],lw=2)
            lines.append(line)
        axes[row,col].set_title(eli.names[i])
        if i>up.shape[0]-3:
            axes[row,col].set_ylim([0,15])
        row+=1
        if (row==2):
            row=0
            col+=1
    f.tight_layout(rect=[0, 0.08, 1, 1])
    # f.legend([leg1,patches[2],patches[1],patches[0]],["Prior distribution","SWMM posterior","Emulator (improved) posterior","Emulator (standard) posterior",], bbox_to_anchor=[0.5, 0.05],loc='center',ncol=2)
    f.legend([leg1,lines[2],lines[1],lines[0]],["Prior distribution","SWMM posterior","Emulator (improved) posterior","Emulator (standard) posterior",], bbox_to_anchor=[0.5, 0.05],loc='center',ncol=2)
    f.savefig("posteriors.pdf",dpi=500)
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

def plot_some_chains(input_file,total_chains,no_pars,chain_no):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    samples=np.genfromtxt(input_file)
    chain_length=samples.shape[0]/total_chains
    plt.clf()
    f,axes=plt.subplots(2,6,figsize=(24,8))
    i=chain_no
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
    filename=create_file_name(["chainz",input_file,chain_no,".pdf"])
    f.savefig(filename)

def separate_some_chains(input_file,total_chains,no_pars,chain_no):
    samples=np.genfromtxt(input_file)
    chain_length=samples.shape[0]/total_chains
    i=chain_no
    chain=samples[i*chain_length:(i+1)*chain_length,:]
    return chain



def plot_triangle(lower_par,upper_par,names,filename):
    import triangle
    inp_file=np.genfromtxt(filename)
    extent=[(1,1)]*lower_par.shape[0]
    for i in np.arange(lower_par.shape[0]):
        extent[i]=(lower_par[i],upper_par[i])
    fig = triangle.corner(inp_file,labels=names,exents=extent)
    filename=create_file_name(["triangle_",filename,".png"])
    fig.savefig(filename)

# calculates the distance of a vector to a closest one from
# a bunch of vectors
def closest_distance(vectors,vector):
    norm=np.inf
    for i in np.arange(vectors.shape[0]):
        if np.linalg.norm(vectors[i]-vector)<norm:
            norm=np.linalg.norm(vectors[i]-vector)
    return norm

def remove_duplicate_rows(matrix):
    a,b=matrix.T
    sorted=np.lexsort((a,b))
    matrix=matrix[sorted]
    matrix_d=np.diff(matrix,axis=0)
    sorted=np.any(matrix_d,axis=1)
    matrix=matrix[sorted]
    return(matrix)

def send_notification(value1):
    import requests
    requests.get('http://maker.ifttt.com/trigger/python_msg/with/key/WXvFafnfwS26tgInnIILn',params={"value1":value1})
