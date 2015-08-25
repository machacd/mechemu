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

def plot_posteriors(eli,file_names):
    from scipy import stats
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    no_files=np.array(file_names).shape[0]
    samples=list()
    for i in np.arange(no_files):
        samples.append(np.genfromtxt(file_names[i]))
    down=eli.lower_bounds
    up=eli.upper_bounds
    no_of_pars=up.shape[0]
    plt.clf()
    f,axes=plt.subplots(int(np.ceil(no_of_pars/3)),3,figsize=(16,16))
    row=0
    col=0
    colors='bgrcmykwbgrcmykw'
    for i in np.arange(up.shape[0]):
        for j in np.arange(no_files):
            axes[row,col].hist([samples[j][:,i]],
                                           bins=100,
                                           alpha=0.8,
                                           range=[down[i],up[i]],
                                           color=colors[j],
                                           histtype='stepfilled',
                                           normed=True)
        axes[row,col].set_title(eli.names[i])
        x=np.arange(down[i],up[i],0.01)
        pri=np.zeros(x.shape[0])
        for j in np.arange(x.shape[0]):
            pri[j]=eli.prior_dist(x[j],i)
        axes[row,col].plot(x,pri)
        row+=1
        if (row==np.ceil(no_of_pars/3)):
            row=0
            col+=1
    f.tight_layout()
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
    a,b,c,d,e,f,g,h=matrix.T
    sorted=np.lexsort((a,b,c,d,e,f,g,h))
    matrix=matrix[sorted]
    matrix_d=np.diff(matrix,axis=0)
    sorted=np.any(matrix_d,axis=1)
    matrix=matrix[sorted]
    return(matrix)
