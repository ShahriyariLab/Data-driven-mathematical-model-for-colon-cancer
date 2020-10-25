'''
QSPClusterAnalysis: script for analysis of  parameter sensitivity and dynamics of 
                    Mathematical model of immune response in colon cancer*

Author: Arkadz Kirshtein, https://sites.google.com/site/akirshtein/
(c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/

Conceptualization of sensitivity algorithm by Wenrui Hao, http://personal.psu.edu/wxh64/


*part of https://github.com/ShahriyariLab/Data-driven-mathematical-model-for-colon-cancer
 If using this or related code please cite 
 Kirshtein, A.; Akbarinejad, S.; Hao, W.; Le, T.; Aronow, R.A.; Shahriyari, L. 
  Data driven mathematical model of colon cancer progression. (Manuscript submitted for publication).
'''

from qspmodel import *
import pandas as pd
import csv
import os
import scipy as sp

# Checking or creating necessary output folders
if not os.path.exists('Data/'):
    os.makedirs('Data/Dynamic/')
    os.makedirs('Data/Sensitivity/')
    os.makedirs('Data/InitialDifference/') 
else:
    if not os.path.exists('Data/Dynamic/'):
        os.makedirs('Data/Dynamic/')
    if not os.path.exists('Data/Sensitivity/'):
        os.makedirs('Data/Sensitivity/')
    if not os.path.exists('Data/InitialDifference/'):
        os.makedirs('Data/InitialDifference/')    


# some global parameters
lmod=[0, 1, 2, 3, 4, 5, 6, 8, 9] #indices of variables in cell data (excluding M0)
clusters=5 #number of clusters

T=3500
t=np.linspace(0, T, 35001)

nvar=Colon_QSP_Functions().nvar # number of variables
nparam=Colon_QSP_Functions().nparam # number of parameters

# Sensitivity Analysis
print('Starting steady state sensitivity analysis')

# reading mean-based cell data for sensitivity analysis
clustercells=pd.read_csv('input/Large_Tumor_cell_data.csv').to_numpy()

# Read the parameter perturbation grid
gridlevel=3
gridname='grid63-level'+str(gridlevel)
data = pd.read_csv('input/'+gridname+'.csv').to_numpy()
w=data[:,0]
x=data[:,1:]
del data

# coefficients for variable sensitivity
lambda0=np.zeros((nvar,2))
lambda0[7,0]=1 # just cancer

# Sensitive parameters to modify for dynamics variability
parmodind=[[],[],[],[],[]]

for cluster in range(clusters):
    print('Starting computations for cluster '+str(cluster+1)+' of 5') 
    filename='V63-grid'+str(gridlevel)+'-cluster-'+str(cluster+1)+'-of-5-results-'

    lambda0[0:9,1]=clustercells[cluster,lmod]/np.sum(clustercells[cluster,lmod]) # all cells    

    QSP_=QSP.from_cell_data(clustercells[cluster])
    params=QSP_.par
    
    print(' Parameters set. Computing steady state sensitivity')

    dudp=np.zeros((nparam,4))
    dif=1e-4 #relative timestep for the finite difference
    sensitivity_radius=1 # perturbation percentage for sensitivity integration
    for k in range(w.size):
        # set parameter sample
        param_sample=params*(1+(sensitivity_radius*1e-2)*x[k,:])
        QSP_.set_parameters(param_sample);   
        #Calculate cancer and total cell sensitivity analytically
        dudp[:,:2]=dudp[:,:2]+w[k]*np.dot(QSP_.Sensitivity(),lambda0)
        #Compute eigenvalue sensitivity numerically
        for i in range(nparam):
            pert=np.zeros(nparam)
            pert[i]=1
            QSP_.set_parameters(param_sample*(1+dif*pert))
            u1=abs(np.real(np.linalg.eigvals(QSP_.Ju(0,QSP_.steady_state()))))
            l1=np.array([min(u1), max(u1)])
            QSP_.set_parameters(param_sample*(1-dif*pert))
            u2=abs(np.real(np.linalg.eigvals(QSP_.Ju(0,QSP_.steady_state()))))
            l2=np.array([min(u2), max(u2)])
            dudp[i,2:]=dudp[i,2:]+w[k]*(l1-l2)/(2*dif*param_sample[i])
        # print('finished node',k+1)
    
    # Updating list of most sensitive parameters
    parmodind[cluster]=np.unique(np.concatenate((np.argsort(np.abs(dudp[:,0]))[::-1][:4],np.argsort(np.abs(dudp[:,1]))[::-1][:4],np.argsort(dudp[:,2])[::-1][:4])))
        
    print(' Writing to file')

    c=csv.writer(open('Data/Sensitivity/'+filename+'sensitivity_steady.csv',"w"))
    c.writerows(dudp)
    del c
    
print('Sensitivity analysis complete')

print('-----------------------------')


#Computations of dynamics
print('Starting dynamics computations')
# reading initial conditions
IC=np.array(pd.read_csv('input/Initial_Conditions.csv'))

for cluster in range(clusters):
    print('Starting computations for cluster '+str(cluster+1)+' of 5') 
    filename='V63-cluster-'+str(cluster+1)+'-of-5-results-'

    lambda0[0:9,1]=clustercells[cluster,lmod]/np.sum(clustercells[cluster,lmod]) # all cells

    QSP_=QSP.from_cell_data(clustercells[cluster])
    params=QSP_.par
    
    print(' Parameters set. Computing the solution')

    u, _ = QSP_.solve_ode(t, IC[cluster], 'given')
    
    wr=np.empty((t.size, 15))
    wr[:,0]=t
    wr[:,1:]=u
    c=csv.writer(open('Data/Dynamic/'+filename+'dat.csv',"w"))  
    c.writerow(['time']+QSP_.variable_names())
    c.writerows(wr)
    del c
    
    print(' Computing dynamics with varying parameters')
    # here we record a minimal and maximal values of each variable across all perturbations of most sensitive parameters
    uvarmax=u
    uvarmin=u
    
    # using parameter indices create arrays for aprameter modification
    parmod=np.zeros((len(parmodind[cluster]),nparam))
    for k in range(len(parmodind[cluster])):
        parmod[k,parmodind[cluster][k]]=1
    # multiplier to perturb parameters by 10% in each direction 
    pert=0.9+0.2*np.array(range(2))
    
    #compute the solution for each perturbation
    for k in range(len(parmodind[cluster])):
        for j in range(pert.size):
            QSP_.set_parameters((1+(pert[j]-1)*parmod[k,:])*params)
            u, _ = QSP_.solve_ode(t, IC[cluster], 'given')
            uvarmax=np.maximum(u,uvarmax)
            uvarmin=np.minimum(u,uvarmin)
        
    print(' Writing to file')
    
    wr[:,1:]=uvarmax
    c=csv.writer(open('Data/Dynamic/'+filename+'varmaxdat.csv',"w"))  
    c.writerow(['time']+QSP_.variable_names()) 
    c.writerows(wr)
    wr[:,1:]=uvarmin
    c=csv.writer(open('Data/Dynamic/'+filename+'varmindat.csv',"w"))  
    c.writerow(['time']+QSP_.variable_names())
    c.writerows(wr)
    del c
    

print('-----------------------------')

print('Computing dynamics with varying initial conditions')  

# less time points to save to conserve disk space
t=np.linspace(0, T, 3501)

# reading initial conditions
IC=np.array(pd.read_csv('input/Small_Tumor_IC_bank.csv'))
patients=IC.shape[0] # number of patients

# creating a model object for each cluster
QSP_=[]
for cluster in range(clusters):
        QSP_.append(QSP.from_cell_data(clustercells[cluster]))

# Recording cluster indices for each patient
c=csv.writer(open('Data/InitialDifference/ClusterNumbers.csv',"w"))  
c.writerow(IC[:,nvar])
del c

for p in range(patients):
    print(' Starting computations for patient', (p+1), 'of', patients) 
    
    u, _ = QSP_[int(IC[p,nvar])-1].solve_ode(t, IC[p,:nvar], 'given')
    
    wr=np.empty((t.size, nvar+1))
    wr[:,0]=t
    wr[:,1:]=u
    c=csv.writer(open('Data/InitialDifference/V63-patient-'+str(p+1)+'-of-'+str(patients)+'-results.csv',"w")) 
    c.writerow(['time']+QSP_[0].variable_names()) 
    c.writerows(wr)
    
    del wr, u, c

print('done')