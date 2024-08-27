import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyemma
from pyemma.util.contexts import settings
import mdtraj as md
import glob
import pandas as pd
from itertools import product
import pickle
import copy
import seaborn as sns
from pandas import DataFrame


#=====================================Load Trajs===============================================

pdb = md.load("./TPS_22500.cp/gn_005_xpr2.noWAT.nc.start.pdb")

files_TPS_2500_cp = glob.glob("./TPS_2500.cp/*.noWAT.align.nc")
files_TPS_5000_cp = glob.glob("./TPS_5000.cp/*.noWAT.align.nc")
files_TPS_7500 = glob.glob("./TPS_7500/*.noWAT.align.nc")
files_TPS_10000 = glob.glob("./TPS_10000/*.noWAT.align.nc")
files_TPS_12500 = glob.glob("./TPS_12500/*.noWAT.align.nc")
files_TPS_15000_cp = glob.glob("./TPS_15000.cp/*.noWAT.align.nc")
files_TPS_17500_cp = glob.glob("./TPS_17500.cp/*.noWAT.align.nc")
files_TPS_20000_cp = glob.glob("./TPS_20000.cp/*.noWAT.align.nc")
files_TPS_22500_cp = glob.glob("./TPS_22500.cp/*.noWAT.align.nc")
files_TPS_25000 = glob.glob("./TPS_25000/*.noWAT.align.nc")
files_TPS_27500_cp = glob.glob("./TPS_27500.cp/*.noWAT.align.nc")
files_TPS_30000 = glob.glob("./TPS_30000/*.noWAT.align.nc")
files_TPS_35000_cp = glob.glob("./TPS_35000.cp/*.noWAT.align.nc")
files_TPS_37500_cp = glob.glob("./TPS_37500.cp/*.noWAT.align.nc")


files_TPS_2500_cp_2 = glob.glob("./TPS_2500.cp-2/*.noWAT.align.nc")
files_TPS_5000_cp_2 = glob.glob("./TPS_5000.cp-2/*.noWAT.align.nc")
files_TPS_7500_cp_2 = glob.glob("./TPS_7500.cp-2/*.noWAT.align.nc")
files_TPS_10000_cp_2 = glob.glob("./TPS_10000.cp-2/*.noWAT.align.nc")
files_TPS_12500_cp_2 = glob.glob("./TPS_12500.cp-2/*.noWAT.align.nc")
files_TPS_15000_cp_2 = glob.glob("./TPS_15000.cp-2/*.noWAT.align.nc")
files_TPS_17500_cp_2 = glob.glob("./TPS_17500.cp-2/*.noWAT.align.nc")
files_TPS_20000_cp_2 = glob.glob("./TPS_20000.cp-2/*.noWAT.align.nc")


files_TPS_2500_cp.sort() 
files_TPS_5000_cp.sort() 
files_TPS_7500.sort()
files_TPS_10000.sort()
files_TPS_12500.sort()
files_TPS_15000_cp.sort() 
files_TPS_17500_cp.sort()
files_TPS_20000_cp.sort()
files_TPS_22500_cp.sort() 
files_TPS_25000.sort()
files_TPS_27500_cp.sort()
files_TPS_30000.sort()
files_TPS_35000_cp.sort()
files_TPS_37500_cp.sort()


files_TPS_2500_cp_2.sort()
files_TPS_5000_cp_2.sort()
files_TPS_7500_cp_2.sort()
files_TPS_10000_cp_2.sort()
files_TPS_12500_cp_2.sort()
files_TPS_15000_cp_2.sort()
files_TPS_17500_cp_2.sort()
files_TPS_20000_cp_2.sort()



files = files_TPS_2500_cp + files_TPS_5000_cp + files_TPS_7500 + files_TPS_10000 + files_TPS_12500 + files_TPS_15000_cp + files_TPS_17500_cp + files_TPS_20000_cp + files_TPS_22500_cp + files_TPS_25000 + files_TPS_27500_cp + files_TPS_30000 + files_TPS_35000_cp + files_TPS_37500_cp + files_TPS_2500_cp_2 + files_TPS_5000_cp_2 + files_TPS_7500_cp_2 + files_TPS_10000_cp_2 + files_TPS_12500_cp_2 + files_TPS_15000_cp_2 + files_TPS_17500_cp_2 + files_TPS_20000_cp_2


coord_feat = pyemma.coordinates.featurizer(pdb) 
P_ind = [33,67,97,130,160,190,224,547,580,613,646,680,711,742,773,807,837,868,898,928,961,992,1025,1056,1087,1211,2052,2085,2115,2148,2181,2212,2246,2280,2314,2345,2375,2405,2435] 
P_index = list(map(lambda x:x-1,P_ind)) # P in resi 1-7+18-36+65-77
coord_feat.add_selection(P_index) #3*39
coord_feat.add_residue_COM(list(range(1,8))+list(range(18,37))+list(range(65,78))) #
coord_feat.add_distances(coord_feat.pairs(P_ind,excluded_neighbors=2),periodic=True) #975
coord_data = pyemma.coordinates.load(files, features=coord_feat)



race_files = ["./TPS_2500.cp/race_history.log","./TPS_5000.cp/race_history.log","./TPS_7500/race_history.log","./TPS_10000/race_history.log","./TPS_12500/race_history.log","./TPS_15000.cp/race_history.log","./TPS_17500.cp/race_history.log","./TPS_20000.cp/race_history.log","./TPS_22500.cp/race_history.log","./TPS_25000/race_history.log","./TPS_27500.cp/race_history.log","./TPS_30000/race_history.log","./TPS_35000.cp/race_history.log", "./TPS_37500.cp/race_history.log","./TPS_2500.cp-2/race_history.log", "./TPS_5000.cp-2/race_history.log", "./TPS_7500.cp-2/race_history.log", "./TPS_10000.cp-2/race_history.log", "./TPS_12500.cp-2/race_history.log", "./TPS_15000.cp-2/race_history.log","./TPS_17500.cp-2/race_history.log", "./TPS_20000.cp-2/race_history.log"]



scale = []
accept_traj = []


for race_log in race_files:
    f=open(race_log,"r")
    lines = f.readlines()
    #scale = []
    #accept_traj = []
    for line in lines:
        if "TRAJS ACCEPTED" in line:
            accept_traj.append(line.split()[5])
            traj_prev = lines[lines.index(line)-7].split()[1]
            if traj_prev.split("_")[2]=="xpr2":
                scale.append([-1,1])
            else:
                scale.append([1,-1])
    f.close()




traj_order = []
traj_order_byf = []
frame_order = []
coord_data_update=[]
for i in range(int(len(coord_data)/2)):
    if int(accept_traj[i]) < 10:
        traj_name ="gn_00"+accept_traj[i]+"_xpr"
    elif int(accept_traj[i]) >= 10 and int(accept_traj[i]) < 100:
        traj_name ="gn_0"+accept_traj[i]+"_xpr"
    else:
        traj_name ="gn_"+accept_traj[i]+"_xpr"
    if scale[i]==[1,-1]:
        coord_subset=np.append(coord_data[i*2+1][::-1],coord_data[i*2],axis=0)
        traj = [traj_name+"2",traj_name+"1"]
        traj_order = traj_order + traj
        frame_a = list(range(len(coord_data[i*2+1])))
        frame_a.reverse()
        frame_b = list(range(len(coord_data[i*2])))
        frame_order = frame_order + [fa+1 for fa in frame_a] + [fb+1 for fb in frame_b]
        traj_order_byf = traj_order_byf + [traj_name+"2"]*len(coord_data[i*2+1]) + [traj_name+"1"]*len(coord_data[i*2])
    else:
        coord_subset=np.append(coord_data[i*2][::-1],coord_data[i*2+1],axis=0)
        traj = [traj_name+"1",traj_name+"2"]
        traj_order = traj_order + traj
        frame_a = list(range(len(coord_data[i*2])))
        frame_a.reverse()
        frame_b = list(range(len(coord_data[i*2+1])))
        frame_order = frame_order + [fa+1 for fa in frame_a] + [fb+1 for fb in frame_b]
        traj_order_byf = traj_order_byf + [traj_name+"1"]*len(coord_data[i*2]) + [traj_name+"2"]*len(coord_data[i*2+1])
    coord_data_update.append(coord_subset)
    

coord_data_concatenated = np.concatenate(coord_data_update)

    
tica_coord = pyemma.coordinates.tica(coord_data_update, lag=100)
tica_coord_output = tica_coord.get_output()
tica_coord_concatenated = np.concatenate(tica_coord_output)


nclusters=75


cluster_tica_coord = pyemma.coordinates.cluster_kmeans(
    tica_coord_output, k=nclusters, max_iter=50, stride=1, fixed_seed=1)
dtrajs_tica_coord_concatenated = np.concatenate(cluster_tica_coord.dtrajs)
cluster_tica_coord_out = cluster_tica_coord.get_output()

#cluster_coord = pyemma.coordinates.cluster_kmeans(
#    coord_data_update, k=nclusters, max_iter=50, stride=1, fixed_seed=1)
#dtrajs_coord_concatenated = np.concatenate(cluster_coord.dtrajs)
#cluster_coord_out = cluster_coord.get_output()


cluster_tica_coord.save('coord.lag100.cluster{}.region3.tica.h5'.format(nclusters), save_streaming_chain=True)


coord_source_coord = pyemma.coordinates.source(files, features=coord_feat)

cluster_samples_coord = cluster_tica_coord.sample_indexes_by_cluster(list(range(0,nclusters)),100)

cluster_samples_update_coord=copy.deepcopy(cluster_samples_coord)
for i in range(len(cluster_samples_coord)):
    for j in range(len(cluster_samples_coord[i])):
        traj=cluster_samples_coord[i][j][0]
        struct=cluster_samples_coord[i][j][1]
        if scale[traj]==[1,-1]:
            trajname=traj_order[traj*2+1]
            if struct>(coord_data[traj*2+1].shape[0]-1):
                frame=struct-coord_data[traj*2+1].shape[0]
                cluster_samples_update_coord[i][j][0]=traj*2
                cluster_samples_update_coord[i][j][1]=frame
            else:
                frame=struct
                cluster_samples_update_coord[i][j][0]=traj*2+1
                cluster_samples_update_coord[i][j][1]=frame
        else:
            trajname=traj_order[traj*2]
            if struct>(coord_data[traj*2].shape[0]-1):
                frame=struct-coord_data[traj*2].shape[0]
                cluster_samples_update_coord[i][j][0]=traj*2+1
                cluster_samples_update_coord[i][j][1]=frame
            else:
                frame=struct
                cluster_samples_update_coord[i][j][0]=traj*2
                cluster_samples_update_coord[i][j][1]=frame        

pyemma.coordinates.save_trajs(
    coord_source_coord,
    cluster_samples_update_coord,
    outfiles=['Cluster{}_100samples'.format(n + 1)+'.coord.cluster{}.region3.tica.pdb'.format(nclusters)
              for n in range(0,nclusters)])
              


df = pd.DataFrame(list(zip(traj_order_byf, frame_order)), columns =['Traj', 'Original_frame']) 
df['Cluster'] = dtrajs_tica_coord_concatenated.tolist()   
df['IC1'] = tica_coord_concatenated[:,0]
df['IC2'] = tica_coord_concatenated[:,1]
df.to_csv("Cluster_Traj_ByFrame.lag=100.cluster={}.tica_coord.region3.tica.csv".format(nclusters))

c_center_coord = cluster_tica_coord.clustercenters
dc_center_coord = pd.DataFrame(c_center_coord)
dc_center_coord.to_csv("Cluster_Center.lag=100.cluster={}.coord.region3.tica.csv".format(nclusters))






############################## MSM ##########################################


msm_coord = pyemma.msm.bayesian_markov_model(cluster_tica_coord.dtrajs, lag=100, dt_traj='0.002 ns')
print('fraction of states used = {:.2f}'.format(msm_coord.active_state_fraction))
print('fraction of counts used = {:.2f}'.format(msm_coord.active_count_fraction))


#fraction of states used = 1.00
#fraction of counts used = 1.00


active_index=list(msm_coord.active_set)
print(active_index)

#nstates=len(active_index)
nstates=36
msm_coord.pcca(nstates)
pcca_samples_coord = msm_coord.sample_by_distributions(msm_coord.metastable_distributions, 100)


pcca_samples_update_coord=copy.deepcopy(pcca_samples_coord)
for i in range(len(pcca_samples_coord)):
    for j in range(len(pcca_samples_coord[i])):
        traj=pcca_samples_coord[i][j][0]
        struct=pcca_samples_coord[i][j][1]
        if scale[traj]==[1,-1]:
            trajname=traj_order[traj*2+1]
            if struct>(coord_data[traj*2+1].shape[0]-1):
                frame=struct-coord_data[traj*2+1].shape[0]
                pcca_samples_update_coord[i][j][0]=traj*2
                pcca_samples_update_coord[i][j][1]=frame
            else:
                frame=struct
                pcca_samples_update_coord[i][j][0]=traj*2+1
                pcca_samples_update_coord[i][j][1]=frame
        else:
            trajname=traj_order[traj*2]
            if struct>(coord_data[traj*2].shape[0]-1):
                frame=struct-coord_data[traj*2].shape[0]
                pcca_samples_update_coord[i][j][0]=traj*2+1
                pcca_samples_update_coord[i][j][1]=frame
            else:
                frame=struct
                pcca_samples_update_coord[i][j][0]=traj*2
                pcca_samples_update_coord[i][j][1]=frame        


pyemma.coordinates.save_trajs(
    coord_source_coord,
    pcca_samples_update_coord,
    outfiles=['PCCA_meta{}_100samples.nstate='.format(n + 1)+str(nstates)+'.lag=100.cluster={}.coord.region3.tica.pdb'.format(nclusters)
              for n in range(msm_coord.n_metastable)])




cktest_coord = msm_coord.cktest(nstates, mlags=10)
#cktest_coord = msm_coord.cktest(nstates)
fig, ax = pyemma.plots.plot_cktest(cktest_coord, figsize=(24,24), dt=0.002, units='ns')
fig.tight_layout()
fig.savefig("Transition_CK.large_fig.n={0}.lag=100.cluster={1}.coord.tica.region3.tica.png".format(nstates,nclusters))


for i in range(0,10):
    data = cktest_coord.estimates[i]
    # plotting the heatmap
    fig, axes = plt.subplots(figsize=(6, 5))
    hm = sns.heatmap(data=data,vmin=0, vmax=1)
    plt.savefig("Transition_estimate_CK.ni="+str(i+1)+".lag=100.cluster={}.coord.region3.tica.png".format(nclusters))
    plt.close()


msm_coord.save('msm_coord.lag100.cluster{}.n={}.coord.region3.tica.h5'.format(nclusters,nstates), model_name='lag100_cluster{}_n={}'.format(nclusters,nstates))


active_index=list(msm_coord.active_set)

list_ndx=[]
dtrajs_structure_index=np.array([]).astype(int)
for i in range(len(dtrajs_tica_coord_concatenated)):
    if dtrajs_tica_coord_concatenated[i] not in active_index:
        list_ndx.append(i)
    else:
        dtrajs_structure_index=np.append(dtrajs_structure_index,i)
        
dtrajs_coord_active_concatenated = np.delete(dtrajs_tica_coord_concatenated,list_ndx)

dtrajs_index=np.array([]).astype(int)

for i in dtrajs_coord_active_concatenated:
    dtrajs_index=np.append(dtrajs_index,active_index.index(i))




fig, axes = plt.subplots(1, nstates, figsize=(nstates*4, 3), sharex=True, sharey=True)
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_coord_concatenated[dtrajs_structure_index, :2].T,
        msm_coord.metastable_distributions[i][dtrajs_index],
        ax=ax,
        cmap='afmhot_r', 
        mask=True,
        cbar_label='metastable distribution {}'.format(i + 1))
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
fig.tight_layout()
fig.savefig("MSM_metastable.nstate="+str(nstates)+".lag=100.cluster={}.coord.region3.tica.png".format(nclusters))



metastable_traj_coord = msm_coord.metastable_assignments[dtrajs_index]

fig, ax = plt.subplots(figsize=(5, 4))
_, _, misc = pyemma.plots.plot_state_map(
    *tica_coord_concatenated[dtrajs_structure_index, :2].T, metastable_traj_coord, ax=ax)
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')
misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
                             for i in range(nstates)])
fig.tight_layout()
fig.savefig("MSM_meta_states.nstate="+str(nstates)+".lag=100.cluster={}.coord.region3.tica.png".format(nclusters))              





fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
pyemma.plots.plot_contour(
    *tica_coord_concatenated[dtrajs_structure_index, :2].T,
    msm_coord.pi[dtrajs_index],
    ax=axes[0],
    mask=True,
    cbar_label='stationary distribution')
pyemma.plots.plot_free_energy(
    *tica_coord_concatenated[:, :2].T,
    weights=np.concatenate(msm_coord.trajectory_weights()),
    ax=axes[1],
    legacy=False)
for ax in axes.flat:
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
axes[0].set_title('Stationary distribution', fontweight='bold')
axes[1].set_title('Reweighted free energy surface', fontweight='bold')
fig.tight_layout()
fig.savefig("FEL_stationary_reweighted.lag=100.cluster={}.n={}.coord.region3.tica.png".format(nclusters,nstates))


eigvec_coord = msm_coord.eigenvectors_right()
print('The first eigenvector is one: {} (min={}, max={})'.format(
    np.allclose(eigvec_coord[:, 0], 1, atol=1e-15), eigvec_coord[:, 0].min(), eigvec_coord[:, 0].max()))

fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_coord_concatenated[dtrajs_structure_index, :2].T,
        eigvec_coord[dtrajs_index, i + 1],
        ax=ax,
        cmap='PiYG',
        cbar_label='{}. right eigenvector'.format(i + 2),
        mask=True)
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
fig.tight_layout()
fig.savefig("MSM_eigenvec.lag=100.cluster={}.n={}.coord.region3.tica.png".format(nclusters,nstates))



print('state\tπ\t\tG/kT')
for i, s in enumerate(msm_coord.metastable_sets):
    p = msm_coord.pi[s].sum()
    print('{}\t{:f}\t{:f}'.format(i + 1, p, -np.log(p)))


mfpt_coord = np.zeros((nstates, nstates))
for i, j in product(range(nstates), repeat=2):
    if msm_coord.metastable_sets[i].shape[0]!=0 and msm_coord.metastable_sets[j].shape[0]!=0:
        mfpt_coord[i, j] = msm_coord.mfpt(msm_coord.metastable_sets[i],msm_coord.metastable_sets[j])


print('MFPT / ns:')
data_mfpt_coord=pd.DataFrame(np.round(mfpt_coord, decimals=2), index=range(1, nstates + 1), columns=range(1, nstates + 1))
data_mfpt_coord.to_csv("Data_MFPT.nstate="+str(nstates)+".lag=100.cluster={}.n={}.coord_only.region3.tica.csv".format(nclusters,nstates))


A = msm_coord.metastable_sets[0]
B = np.concatenate(msm_coord.metastable_sets[1:])
print('MFPT 1 -> other: ({:6.1f} ± {:5.1f}) ns'.format(
    msm_coord.sample_mean('mfpt', A, B), msm_coord.sample_std('mfpt', A, B)))
print('MFPT other -> 1: ({:.1f} ± {:5.1f}) ns'.format(
    msm_coord.sample_mean('mfpt', B, A), msm_coord.sample_std('mfpt', B, A)))


for start in list(range(nstates)):
    for final in list(range(start+1,nstates)):
        A = msm_coord.metastable_sets[start]
        B = msm_coord.metastable_sets[final]
        if len(A)==0 or len(B)==0:
            continue
        flux = pyemma.msm.tpt(msm_coord, A, B)
        cg, cgflux = flux.coarse_grain(msm_coord.metastable_sets)

        fig, ax = plt.subplots(figsize=(5, 4))
        pyemma.plots.plot_contour(*tica_coord_concatenated[dtrajs_structure_index, :2].T,flux.committor[dtrajs_index],cmap='brg',ax=ax,mask=True,cbar_label=r'committor $\mathcal{S}_%d \to \mathcal{S}_%d$' % (start + 1, final + 1))
        fig.tight_layout()
        fig.savefig("MSM_committor."+str(start)+"-"+str(final)+"-n"+str(nstates)+".lag=100.cluster={}.n={}.coord.region3.tica.png".format(nclusters,nstates))



def its_separation_err(ts, ts_err):
    """
    Error propagation from ITS standard deviation to timescale separation.
    """
    return ts[:-1] / ts[1:] * np.sqrt(
        (ts_err[:-1] / ts[:-1])**2 + (ts_err[1:] / ts[1:])**2)


nits = 20

timescales_mean_coord = msm_coord.sample_mean('timescales', k=nits)
timescales_std_coord = msm_coord.sample_std('timescales', k=nits)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].errorbar(
    range(1, nits + 1),
    timescales_mean_coord, 
    yerr=timescales_std_coord, 
    fmt='.', markersize=10)
axes[1].errorbar(
    range(1, nits),
    timescales_mean_coord[:-1] / timescales_mean_coord[1:], 
    yerr=its_separation_err(
        timescales_mean_coord, 
        timescales_std_coord), 
    fmt='.', 
    markersize=10,
    color='C0')

for i, ax in enumerate(axes):
    ax.set_xticks(range(1, nits + 1))
    ax.grid(True, axis='x', linestyle=':')
    
axes[0].axhline(msm_coord.lag * 0.002, lw=1.5, color='k')
axes[0].axhspan(0, msm_coord.lag * 0.002, alpha=0.3, color='k')
axes[0].set_xlabel('implied timescale index')
axes[0].set_ylabel('implied timescales / ns')
axes[1].set_xticks(range(1, nits))
axes[1].set_xticklabels(
    ["{:d}/{:d}".format(k, k + 1) for k in range(1, nits + 2)],
    rotation=45)
axes[1].set_xlabel('implied timescale indices')
axes[1].set_ylabel('timescale separation')
fig.tight_layout()
fig.savefig("MSM_ITS_separation.lag=100.cluster={}.coord.n={}.region3.tica.png".format(nclusters,nstates))



all_cluster=[]
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
for i, k in enumerate([25, 75, 200]):
    cluster = pyemma.coordinates.cluster_kmeans(tica_coord_output, k=k, max_iter=50, stride=1)
    all_cluster.append(cluster)
    pyemma.plots.plot_density(*tica_coord_concatenated[:,:2].T, ax=axes[0, i], cbar=False, alpha=0.1)
    axes[0, i].scatter(*cluster.clustercenters[:,:2].T, s=15, c='C1')
    axes[0, i].set_xlabel("IC1")
    axes[0, i].set_ylabel("IC2")
    axes[0, i].set_title('k = {} centers'.format(k))
    pyemma.plots.plot_implied_timescales(
        pyemma.msm.its(cluster.dtrajs,lags=[1, 5, 25, 100, 250], nits=10, errors='bayes'),
        ax=axes[1, i], units='ns', dt=0.002)
    #axes[1, i].set_ylim(1, 2000)
fig.tight_layout()
fig.savefig("Lag_from_MSMestimate.cluster_series.region3.tica.png")


fig, axes = plt.subplots(2, 3, figsize=(12, 6))
for i, k in enumerate([20, 50, 100]):
    cluster = pyemma.coordinates.cluster_kmeans(data, k=k, max_iter=50, stride=10)
    pyemma.plots.plot_density(*data_concatenated.T, ax=axes[0, i], cbar=False, alpha=0.1)
    axes[0, i].scatter(*cluster.clustercenters.T, s=15, c='C1')
    axes[0, i].set_xlabel('$\Phi$')
    axes[0, i].set_ylabel('$\Psi$')
    axes[0, i].set_title('k = {} centers'.format(k))
    pyemma.plots.plot_implied_timescales(
        pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 5, 10, 20, 50], nits=4, errors='bayes'),
        ax=axes[1, i], units='ps')
    axes[1, i].set_ylim(1, 2000)
fig.tight_layout()



              
