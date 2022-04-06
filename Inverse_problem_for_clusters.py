import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNetCV, enet_path
from scipy.io import loadmat






#for visualization, parameters from matlab
viz = True
dur_preceding = 20
win_length = 401
Wave_dur = 20


#parameters to change
ndir_minimization = 1
num_cl = 1
name = '0026'
dir_path = '/home/ksasha/projects/K_waves/results/'


data = np.loadtxt(dir_path+name+'/spike'+'.csv', delimiter=",")
waves = np.loadtxt(dir_path+name+'/wave'+'.csv', delimiter=",")
ndir = np.loadtxt(dir_path+name+'/ndir'+'.csv', delimiter=",").astype(int)
n_strts = np.loadtxt(dir_path+name+'/nstrts'+'.csv', delimiter=",").astype(int)
cluster_ind = loadmat(dir_path+name+'/time_and_sources'+'.mat')['cluster_ind'][0,:]

# speeds = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
nspikes = len(n_strts) # number of spikes
R = int(data.shape[1]/nspikes) # number of shifts in sliding window
S = int(waves.shape[1]/data.shape[0]) # number of propagation speeds
T = data.shape[0]

regression = ElasticNetCV(l1_ratio = 0.5, positive=True, normalize=True, cv=5, max_iter=10000) # elastic net regression
cumdir = 0

nspikes_in_cl = sum(cluster_ind == num_cl)


bestind = np.zeros([nspikes_in_cl,], dtype=int) # indices of best speeds
best_strt = np.zeros([nspikes_in_cl,],dtype = int)
numdir = np.zeros(nspikes_in_cl,) # number of nonzero directions in optimum
#bestcoef = np.zeros([nspikes, int(ndir.T.max())])
finalscore = np.zeros(nspikes_in_cl,)

bestshifts = np.zeros((nspikes_in_cl,S),dtype = int)

all_coefs = list()
cum_strt_vert =0

ind_sp_cl = 0
for ind_sp in range(0,nspikes):
    
    if cluster_ind[ind_sp] != num_cl:
       
        for var in range(round(n_strts[ind_sp])):
            cumdir += int(ndir[cum_strt_vert])
            cum_strt_vert += 1
    else:
        datasp = data[:, ind_sp*R:(ind_sp + 1)*R]  # wavedur* Nsensors x nshifts
        
        
        R2 = np.zeros(round(n_strts[ind_sp]))
        for var in range(round(n_strts[ind_sp])):
            Ndir = int(ndir[cum_strt_vert]) # number of propagation directions
            wavessp = waves[cumdir:cumdir + Ndir]  # waves for this spike Ndir x 10Speeds*102ch*20ms
    
            coefs = np.zeros([R,S,Ndir]) # regression coefficients
            intercept = np.zeros([R,S]) # regression intercept
            score = np.zeros([R,S]) # R-squared scores
            nzdir = np.zeros([R,S]) # number of nonzero directions
            y_pred = np.zeros([R,S,data.shape[0]]) # predicted spikes
            cumdir = cumdir+Ndir
            cum_strt_vert +=1
            for r in range(0,R):
                DataLin = datasp[:,r] # for each time slice
                for s in range(0,int(S)):
                    wavesspeed = wavessp[:,(s*T):((s+1)*T)]
                    regression.fit(wavesspeed.T, DataLin)
                    coefs[r,s,:] = regression.coef_
                    intercept[r,s] = regression.intercept_
                    score[r,s] = regression.score(wavesspeed.T, DataLin)
                    y_pred[r,s,:] = regression.predict(wavesspeed.T)
                    nzdir[r,s] = np.sum(coefs[r,s,:]!=0)
    
            bestshifts[ind_sp_cl] = score.argmax(axis = 0) # best shifts for each speed
            bestscore = score.max(axis = 0) # corresponding scores
            
            bestdir = np.zeros(S,) # corresponding number of nonzero directions (without intercept)
            for s in range(0,S):
                bestdir[s] = nzdir[bestshifts[ind_sp_cl,s],s]
            
            score_sort_ind = (-(bestscore-2*(bestscore==0))).argsort() # indices of sorted scores for all speeds
            dir_sort_ind = (bestdir[score_sort_ind[0:ndir_minimization]]).argsort() # indices of sorted number of nonzero directions for top-3 scores
            bestind[ind_sp_cl] = score_sort_ind[dir_sort_ind[0]]
            R2[var] = bestscore[bestind[ind_sp_cl]]
        n_highest_R2 = np.argmax(R2)
        best_strt[ind_sp_cl] = n_highest_R2
    
    
        if n_highest_R2 != len(R2)-1:
            Ndir = int(ndir[cum_strt_vert+n_highest_R2-round(n_strts[ind_sp])]) # number of propagation directions
            ndirs_back = round(sum(ndir[cum_strt_vert+n_highest_R2-round(n_strts[ind_sp]):cum_strt_vert]))
            wavessp = waves[cumdir-ndirs_back:cumdir-ndirs_back + Ndir]  # waves for this spike Ndir x 10Speeds*102ch*20ms
        
            coefs = np.zeros([R,S,Ndir]) # regression coefficients
            intercept = np.zeros([R,S]) # regression intercept
            score = np.zeros([R,S]) # R-squared scores
            nzdir = np.zeros([R,S]) # number of nonzero directions
            y_pred = np.zeros([R,S,data.shape[0]]) # predicted spikes
            for r in range(0,R):
                DataLin = datasp[:,r] # for each time slice
                for s in range(0,int(S)):
                    wavesspeed = wavessp[:,(s*T):((s+1)*T)]
                    regression.fit(wavesspeed.T, DataLin)
                    coefs[r,s,:] = regression.coef_
                    intercept[r,s] = regression.intercept_
                    score[r,s] = regression.score(wavesspeed.T, DataLin)
                    y_pred[r,s,:] = regression.predict(wavesspeed.T)
                    nzdir[r,s] = np.sum(coefs[r,s,:]!=0)
        
            bestshifts[ind_sp_cl] = score.argmax(axis = 0) # best shifts for each speed
            bestscore = score.max(axis = 0) # corresponding scores
            
            bestdir = np.zeros(S,) # corresponding number of nonzero directions (without intercept)
            for s in range(0,S):
                bestdir[s] = nzdir[bestshifts[ind_sp_cl,s],s]
        
            score_sort_ind = (-(bestscore-2*(bestscore==0))).argsort() # indices of sorted scores for all speeds
            dir_sort_ind = (bestdir[score_sort_ind[0:ndir_minimization]]).argsort() # indices of sorted number of nonzero directions for top-3 scores
                
            bestind[ind_sp_cl] = score_sort_ind[dir_sort_ind[0]] # index of best speed
        numdir[ind_sp_cl] = bestdir[bestind[ind_sp_cl]] # number of nonzero directions in optimum
       
        finalscore[ind_sp_cl] = bestscore[bestind[ind_sp_cl]] # R-squared value in optimum
        all_coefs.append(coefs[bestshifts[ind_sp_cl,bestind[ind_sp_cl]],bestind[ind_sp_cl],:])
        print(ind_sp_cl)
        ind_sp_cl += 1
    
        #plt.figure()
        #plt.plot(DataLin)
        #plt.plot(y_pred[bestshifts[bestind[ind_sp]], bestind[ind_sp], :])
        #plt.title(['R-squared = ', str(finalscore[ind_sp])])
        
    
    
    
    
    
    
    
    
    
    
    
    
    
########################################################
    
    
    
#####################################################
    
    
    
from matplotlib import pyplot as plt
import numpy as np


dataraw = np.loadtxt(dir_path+name+'/spike_raw'+'.csv', delimiter=",")
directions= np.loadtxt(dir_path+name+'/directions'+'.csv', delimiter=",")


step = 0.05
bins = np.arange(0,1.0,step)

R2_hist = np.histogram(finalscore,bins)


fig,ax= plt.subplots(figsize = [5.8,5],dpi = 600)


clrs = [[0.6350,    0.0780,    0.1840],
    [0.3010,    0.7450,    0.9330],
    [0.4660,    0.6740,    0.1880],
    [0.4940,    0.1840,    0.5560],
    [0.9290,    0.6940,    0.1250],
    [0.8500,    0.3250,    0.0980],
         [0,    0.4470,    0.7410],
         [0.1,0.1,0.1],
         [0, 0.5, 0.5], 
         [0.5,1.0,0.5],[0.6,0.1,0.5]]
         
clr =clrs[num_cl-1]




ax.bar(bins[:-1],R2_hist[0],width = step,align='edge',linewidth = 0.5,edgecolor = 'black',color=clr)
plt.xticks(np.arange(0.0,1.0,0.2),fontsize = 11)
plt.yticks(fontsize = 11)
plt.xlim(0.0,0.9)
plt.xlabel('$R^2$ values',fontsize = 14)
plt.ylabel('Count of spikes',fontsize = 14)

dir_bins = np.arange(1,7);
dir_hist = np.histogram(numdir,dir_bins)

fig,ax= plt.subplots(figsize = [5.8,5],dpi = 600)


ax.bar(dir_bins[:-1],dir_hist[0],width = 1,linewidth = 0.5,edgecolor = 'black',color=clr)

plt.xticks(np.arange(1,6,dtype = 'int'),fontsize = 11)
plt.yticks(fontsize = 11)
plt.xlabel('Number of used directions',fontsize = 14)
plt.ylabel('Count of spikes',fontsize = 14)
plt.text(5,0.8*ax.get_ylim()[1],str(num_cl),fontsize = 40)




  
if viz:
    Nch = data.shape[0]//Wave_dur
        
    slide = np.zeros((win_length,Nch))
    cumdir = 0
        
        
    ndir = ndir.astype(int)
 
    inds = np.where(cluster_ind == num_cl)[0]
    for ind_sp in range(0,len(finalscore)):
     
            
        fig,ax =plt.subplots()
           
        slide = dataraw[:,inds[ind_sp]].reshape((win_length,Nch))
                
            
        idx = np.argsort(-np.max(np.abs(slide[win_length//2-30:win_length//2+30,:]),axis = 0))
        ax.plot(np.arange(0,win_length),slide[:,idx[:30]],linewidth = 1)
            
            
        ax.fill_between((bestshifts[ind_sp,bestind[ind_sp]]+win_length//2-dur_preceding,
                             win_length//2-dur_preceding+Wave_dur+bestshifts[ind_sp,bestind[ind_sp]]),
                            ax.get_ylim()[0],ax.get_ylim()[1],color = 'green',alpha = 0.3)
        plt.title(str(ind_sp))
            
          

    
from scipy.io import savemat

#bestshifts[ind_sp,bestind[ind_sp]],bestind[ind_sp]


savemat(dir_path+name+'/wave_prop_cl'+str(num_cl)+'.mat', mdict ={'bestind':bestind,
                         'all_coefs':all_coefs,'bestshifts':bestshifts,'best_strt':best_strt})



savemat(dir_path+name+'/R2val_cl'+str(num_cl)+'.mat', mdict ={'finalscore':finalscore})


# plt.figure()
# ax = plt.axes(projection='3d')
# ax.quiver(0,0,0,nz_wave_dir[0], nz_wave_dir[1], nz_wave_dir[2])
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
