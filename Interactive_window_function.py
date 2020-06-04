import numpy as np
from bqplot import pyplot as plt
from bqplot import marks as mk
from scipy import stats
from scipy import integrate
#load data to be used as global variables
pk=np.load('./data/pk.npy',allow_pickle = True)
M=np.load('./data/M.npy',allow_pickle = True)
c=np.load('./data/c.npy',allow_pickle = True)

#set the cosmology
h=0.7
omega_m=0.2793
omega_l=0.7207
N=1024

def log_norm(k,R,mu,sigma):
    y=np.exp(-0.5*(np.log10(k*R)-mu)**2/sigma**2)

    return(y)

def Critical_density(scale_factor,omega_matt,omega_lambda,hubble_const, little_h=True):
    #'hubble constant' is little h

    g=4.301*10**(-9) #big g in Mpc M_sun (km s^-1)^2
    h=hubble_const*np.sqrt(omega_matt*scale_factor**(-3)+omega_lambda)

    if little_h==True:
        return(3/(8*np.pi*g)*(100)**2)
    elif little_h==False:
        return(3/(8*np.pi*g)*(100*h)**2)

def density_rms_non_parametric(R,pk,window_function):

        rms=integrate.simps(pk[:,0]**2*pk[:,1]*window_function**2/(2*np.pi**2) , pk[:,0])

        return(rms)

def error_non_parametric(pk,M,c,window_function):
    #assume that *window_function is a 2d array with the the first dimension being k and the second being the window function
    #also assume the all of the ks from pk are the same
    #window function is assumed to be in units of kR not k


    peak_height=np.array([])
    conc=np.array([])
    for i in range(len(M)):
        #test if the length of the array is zero (this breaks the functions used)
        if len(M)==0:
            continue

        #calculate R200 from M200
        rho_crit=Critical_density(1.0,omega_m,omega_l,h)
        R=(M[i]/(4/3*rho_crit))**(1/3)


        
        for j in range(len(R)):
            window_sampled=np.interp(pk[0][:,0]*R[j],window_function[:,0],window_function[:,1])
            #edit so that anyting above the max k in window_function is set to zero
            #window_sampled[pk[0][:,0]*R[j]>np.max(window_function[:,0])]=0
            

            rms=density_rms_non_parametric(R,pk[i],window_sampled)

            peak_height=np.append(peak_height,1.68/rms**0.5)
            conc=np.append(conc,c[i][j])


    #fit c-peak height relation with a power law
    rho=stats.spearmanr(conc,peak_height)[0]

    return(rho)

def foo(change):
    #function to update when data is moved

    #grab x and y data
    x_dat=scatter.x
    y_dat=scatter.y

    #sort in terms of x
    y_dat=y_dat[np.argsort(x_dat)]
    x_dat=x_dat[np.argsort(x_dat)]

    #update line plot data
    line.x=x_dat
    line.y=y_dat   

    window_fun=np.empty((len(y_dat),2))
    window_fun[:,0]=10**x_dat
    window_fun[:,1]=y_dat

    rho=error_non_parametric(pk,M,c,window_fun)
    fig.title='log(1-rho)=%.2f'%(np.log10(1-np.abs(rho)))


def plotterer(k,W):

    global fig
    fig=plt.figure()
    
    #set scatter line and text to be global so can be updated by any function
    global scatter
    global line
    scatter=plt.scatter(np.log10(k),W,default_size=200)
    line=plt.plot(np.log10(k),W)

    window_fun=np.empty((len(k),2))
    window_fun[:,0]=k
    window_fun[:,1]=W
    rho=error_non_parametric(pk,M,c,window_fun)
    fig.title='log(1-rho)=%.2f'%(np.log10(1-np.abs(rho)))


    scatter.observe(foo,'y')
    scatter.enable_move=True
    plt.show()