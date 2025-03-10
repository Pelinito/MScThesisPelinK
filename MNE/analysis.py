import mne 
import pandas as pd 
import numpy as np 
import scipy.io as sio
import h5py
import os
import matplotlib.pyplot as plt
plt.ion()
import matplotlib
#matplotlib.use('WXAgg') 
from mne.time_frequency import tfr_morlet
import pickle
import scipy
import random 
#import seaborn as sns

random.seed(19)


def weightedAverages(subjNo):

    head_epochs = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/subj{subjNo}_head-epo.fif')
    bgrd_epochs = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/subj{subjNo}_bgrd-epo.fif')

    info_head = head_epochs.info
    info_bgrd = bgrd_epochs.info

    data_head = head_epochs.get_data()
    data_bgrd = bgrd_epochs.get_data()

    weight = len(head_epochs) / len(bgrd_epochs) + len(head_epochs)

    weighted_head = weight * data_head
    weighted_bgrd = weight * data_bgrd

    weighted_epoch_head = mne.EpochsArray(weighted_head, info_head, tmin=-0.875)
    weighted_epoch_bgrd = mne.EpochsArray(weighted_bgrd, info_bgrd, tmin=-0.875)

    return weighted_epoch_head, weighted_epoch_bgrd


for subj in range(1,19):
    print('working on subject', subj)
    weighted_head, weighted_bgrd = weightedAverages(subj)
    weighted_head.save(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_head_{subj}-epo.fif', overwrite=True)
    weighted_bgrd.save(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_bgrd_{subj}-epo.fif', overwrite=True)

head_evo_list = []
bgrd_evo_list = []

for subj in range(1,19):
    head_epo = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_head_{subj}-epo.fif')
    bgrd_epo = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_bgrd_{subj}-epo.fif')

    head_evo = head_epo.average()
    bgrd_evo = bgrd_epo.average()   

    head_evo_list.append(head_evo)
    bgrd_evo_list.append(bgrd_evo)


grandAverage_head = mne.grand_average(head_evo_list)
grandAverage_bgrd = mne.grand_average(bgrd_evo_list)

fig1 = grandAverage_head.plot(spatial_colors=True)
grandAverage_head.plot_topomap(times=[0.1, 0.2, 0.3, 0.4, 0.5], ch_type='eeg', time_unit='s')
fig1.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_head_grandAverage.png')

fig2 = grandAverage_bgrd.plot(spatial_colors=True)
fig2.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_bgrd_grandAverage.png')

picks = ['Oz']#, 'O1', 'O2']
fig3 = grandAverage_head.plot(picks=picks, spatial_colors=True, ylim=dict(eeg=[-300, 300]), window_title='Grandaverage of weighted head conditions in occipital channels across subjects')
#fig3.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_head_grandAverage_occChan.png')
fig3.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_head_grandAverage_Oz.png')

fig4 = grandAverage_bgrd.plot(picks=picks, spatial_colors=True, ylim=dict(eeg=[-300, 300]), window_title='Grandaverage of weighted background in occipital channels across subjects')
#fig4.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_bgrd_grandAverage_occChan.png')
fig4.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_bgrd_grandAverage_Oz.png')

# Plot the difference between the grand averages
diff_evo_data = grandAverage_head.get_data() - grandAverage_bgrd.get_data()
diff_evo = mne.EvokedArray(diff_evo_data, grandAverage_head.info)

fig5 = diff_evo.plot(spatial_colors=True, picks=picks)
fig5.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_diff_grandAverage.png')

figx = grandAverage_head.plot_image(cmap='RdBu_r',show_names= True, show=True)
# y_tick_labels = head_evo.ch_names

# figx.yticks(y_tick_labels)
figx.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_head_grandAverage_Oz_image.png')

figx = grandAverage_bgrd.plot_image(cmap='RdBu_r', show_names=True, show=True)
# y_tick_labels = head_evo.ch_names
# figx.yticks(y_tick_labels)
figx.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_bgrd_grandAverage_Oz_image.png')


figx = grandAverage_head.plot_topomap(times=[0.0, 0.1, 0.2])
figx.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_head_grandAverage_topomap.png')

## Plot the Oz channel activity on one plot by extracting data from the evoked lists
import matplotlib.cm as cm
from matplotlib.colors import Normalize


head_Oz = []
bgrd_Oz = []

for head_evo, bgrd_evo in zip(head_evo_list, bgrd_evo_list):
    head_Oz.append(head_evo.data[head_evo.ch_names.index('Oz')])
    bgrd_Oz.append(bgrd_evo.data[bgrd_evo.ch_names.index('Oz')])

head_Oz = np.array(head_Oz)
bgrd_Oz = np.array(bgrd_Oz)

head_Oz_mean = np.mean(head_Oz, axis=0)
bgrd_Oz_mean = np.mean(bgrd_Oz, axis=0)

cmap = cm.get_cmap('RdBu_r')
norm = Normalize(vmin=0, vmax=1)
colors = cmap(norm([0.2, 0.8]))

plt.figure(figsize=(7, 4))

plt.plot(head_evo.times, head_Oz_mean, color=colors[0], label='Head condition', linewidth=1)
plt.plot(bgrd_evo.times, bgrd_Oz_mean, color=colors[1], label='Background condition', linewidth=1)

plt.xlabel('Time (ms)')
plt.ylabel('Amplitude (µV)')
#plt.xline(875, color='black', linestyle='--')


y_ticks = [-0.0002, -0.0001, 0, 0.0001, 0.0002]
y_tick_labels = ['-0.2', '-0.1', '0', '0.1', '0.2']
x_ticks = np.arange(-0.8, 1.2, 0.2)
x_tick_labels = ['-800', '-600', '-400', '-200', '0', '200', '400', '600', '800', '1000']
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)
plt.axvline(0, color='black', linestyle='--')    
plt.title('Oz channel activity in Head and Background conditions')
plt.legend()
plt.set_cmap("RdBu_r")
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/weighted_Oz_channel_activity.png')


# # TIME-FREQUENCY ANALYSIS # # 

## Concatenation

# concatenate the head epochs

# weighted_concat_head = []

# weighted_epoch_head_1 = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_head_1-epo.fif')
# weighted_concat_head.append(weighted_epoch_head_1)
# for subj in range(2,19):
#     weighted_epoch_head = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_head_{subj}-epo.fif')
#     weighted_concat_head.append(weighted_epoch_head)

# weighted_concat_head = mne.concatenate_epochs(weighted_concat_head)

# weighted_concat_head.save(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_concat_head-epo.fif', overwrite=True)

# concatenate the background epochs

# weighted_concat_bgrd = []
# weighted_epoch_bgrd_1 = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_bgrd_1-epo.fif')   
# weighted_concat_bgrd.append(weighted_epoch_bgrd_1)
# for subj in range(2,19):
#     weighted_epoch_bgrd = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_bgrd_{subj}-epo.fif')
#     weighted_concat_bgrd.append(weighted_epoch_bgrd)

# weighted_concat_bgrd = mne.concatenate_epochs(weighted_concat_bgrd)
# weighted_concat_bgrd.save(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_concat_bgrd-epo.fif', overwrite=True)

# **
# **
##### RUN FROM HERE ##### 
# **
# **


## LOAD CONCATENATED EPOCHS
# def loadConcatData():
#     weighted_concat_head = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_concat_head-epo.fif')
#     weighted_concat_bgrd = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_concat_bgrd-epo.fif')

#     return weighted_concat_head, weighted_concat_bgrd

# weighted_concat_head, weighted_concat_bgrds = loadConcatData()

## TFR FROM CONCATENATED EPOCHS
freqs = np.arange(4, 12.5, 0.5)
#freqs = np.arange(0.5, 45.5, 0.5)
n_cycles = 5 # or 3, 5

tfrs_head = []
tfrs_bgrd = []

for subject in range(1, 19):

    weighted_epoch_head = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_head_{subject}-epo.fif')
    weighted_epoch_bgrd = mne.read_epochs(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/weighted_bgrd_{subject}-epo.fif')

    tfr_head = tfr_morlet(weighted_epoch_head, freqs=freqs, n_cycles=n_cycles, return_itc=False, average=True, output='power')
    tfr_bgrd = tfr_morlet(weighted_epoch_bgrd, freqs=freqs, n_cycles=n_cycles, return_itc=False, average=True, output='power')

    tfrs_head.append(tfr_head)
    tfrs_bgrd.append(tfr_bgrd)

## Save the Arrays
import pickle

# Save the combined TFRs list to a single file using the standard pickle module
#combined_tfrs = {'head': tfrs_head, 'bgrd': tfrs_bgrd}
combined_tfrs_cyc1 = {'head': tfrs_head, 'bgrd': tfrs_bgrd}
combined_tfrs_cyc5 = {'head': tfrs_head, 'bgrd': tfrs_bgrd}

# with open('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/combined_tfrs.pkl', 'wb') as f:
#     pickle.dump(combined_tfrs, f)

with open('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/combined_tfrs_cyc1.pkl', 'wb') as f:
    pickle.dump(combined_tfrs_cyc1, f)

with open('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/combined_tfrs_cyc5.pkl', 'wb') as f:
    pickle.dump(combined_tfrs_cyc5, f)

    

## Take the Grand Average of the TFRs
    
grandAverage_bgrd = mne.grand_average(tfrs_bgrd)
grandAverage_head = mne.grand_average(tfrs_head)

## Plot the Grand Average TFRs
fig = grandAverage_bgrd.plot(picks=['Oz'], baseline=(-0.5, -0.2), mode='logratio', title='ERSP of Background Condition at Oz - Averaged across subjects')
#fig[0].savefig(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/grandAverage_bgrd_TFR_{n_cycles}cyc_{min(freqs)}-{max(freqs)}.png')
fig[0].savefig(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/grandAverage_bgrd_TFR_{n_cycles}cyc_{min(freqs)}-{max(freqs)}_Oz.png')

fig = grandAverage_head.plot(picks=['Oz'], baseline=(-0.5, -0.2), mode='logratio', title='ERSP of Head Condition at Oz - Averaged across subjects')
fig[0].savefig(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/grandAverage_head_TFR_{n_cycles}cyc_{min(freqs)}-{max(freqs)}_Oz.png')    

## LOAD the TFR data

with open('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/eeglabtomne/combined_tfrs.pkl', 'rb') as f:
    loaded_tfrs = pickle.load(f)

# Access the data
tfrs_head = loaded_tfrs['head']
tfrs_bgrd = loaded_tfrs['bgrd']
# Cut the time axis
for tfr in tfrs_head:
    tfr.crop(tmin=-0.8, tmax=0.975)
for tfr in tfrs_bgrd:
    tfr.crop(tmin=-0.8, tmax=0.975)

# Convert the data to numpy arrays
data_head = np.array([tfr.data for tfr in tfrs_head]) # Shape: (n_subjects, n_channels, n_freqs, n_times)
data_bgrd = np.array([tfr.data for tfr in tfrs_bgrd]) # Shape: (n_subjects, n_channels, n_freqs, n_times)

# Grand average
grandAverage_head = mne.grand_average(tfrs_head)
grandAverage_bgrd = mne.grand_average(tfrs_bgrd)

grAverage_head_data = grandAverage_head.data
grAverage_bgrd_data = grandAverage_bgrd.data

## ONE SUBJECT -- FIRST SUBJECT
## Cut the time axis
data_single_head = tfrs_head[0]
data_single_bgrd = tfrs_bgrd[0]

# data_single_head.crop(tmin=-0.8, tmax=0.975)
# data_single_bgrd.crop(tmin=-0.8, tmax=0.975)

# Pick channel
data_single_head = data_single_head.pick_channels(['Oz'])
data_single_bgrd = data_single_bgrd.pick_channels(['Oz'])
# Plot the data
fig = data_single_head.plot([0], baseline=(-0.5, -0.2), mode='logratio', title='ERSP of Head Condition - Single Subject')
#save the plot
fig[0].savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/single_subject_head_OZ.png')   
fig = data_single_bgrd.plot([0], baseline=(-0.5, -0.2), mode='logratio', title='ERSP of Background Condition - Single Subject')
fig[0].savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/single_subject_bgrd_OZ.png')

## Difference in one subject
difference_single = data_single_head
difference_single.data = data_single_head.data - data_single_bgrd.data
fig = difference_single.plot([0], baseline=(-0.5, -0.2), mode='logratio', title='Difference of activity between Head and Background conditions - Single Subject')
fig[0].savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/single_subject_difference_OZ.png')

# T-test
t_result_single = scipy.stats.ttest_ind(data_single_head.data, data_single_bgrd.data, axis=None)
t_result_single_channel = scipy.stats.ttest_ind(data_single_head.data, data_single_bgrd.data, axis=0)
t_result_single_time = scipy.stats.ttest_ind(data_single_head.data, data_single_bgrd.data, axis=2)
# Plot the T-test channel results
plt.clf()
plt.figure()
plt.imshow(t_result_single_channel.statistic, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('T-test between Head and Background conditions on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/T-test_single_channel.png')

# Variance
var_single = np.var(difference_single.data, axis=None)
var_single_time = np.var(difference_single.data, axis=2)
var_single_freq = np.var(difference_single.data, axis=1)
var_single_channels = np.var(difference_single.data, axis=0)

# Plot variance in single data time
plt.clf()
plt.figure()
plt.imshow(var_single_time, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the difference between Head and Background conditions on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.colorbar()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_single_time.png')

# plot variance in single data channel

plt.clf()
plt.figure()
plt.imshow(var_single_channels, extent=[np.min(data_single_head.times), np.max(data_single_head.times), np.min(data_single_head.freqs), np.max(data_single_head.freqs)], aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the difference between Head and Background conditions on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_single_channels.png')


# plot the T-test time results
plt.clf()
plt.figure()
plt.imshow(t_result_single_time.statistic, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('T-test between Head and Background conditions on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(data_single_head.ch_names), 1), data_single_head.ch_names)
plt.colorbar()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/T-test_single_time.png')


### T-test ###
from mne.stats import f_oneway
#from scipy.stats import f_oneway


# T-test BETWEEN conditions
#f_stats, p_value = scipy.stats.f_oneway(grAverage_head_data, grAverage_bgrd_data, axis=None)
t_result = scipy.stats.ttest_ind(grAverage_head_data, grAverage_bgrd_data, axis=None)
t_result_time = scipy.stats.ttest_ind(grAverage_head_data, grAverage_bgrd_data, axis=2)
t_result_freq = scipy.stats.ttest_ind(grAverage_head_data, grAverage_bgrd_data, axis=1)
t_result_channels = scipy.stats.ttest_ind(grAverage_head_data, grAverage_bgrd_data, axis=0)
# plot t_result_time

plt.figure(figsize=(8, 10))
plt.imshow(t_result_time.statistic, aspect='auto', origin='lower',cmap='RdBu_r',vmin=-80, vmax=80)
plt.title('T-test between Head and Background conditions on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/T-test_time.png')

# Plot p-values

plt.figure(figsize=(8, 10))
plt.imshow(t_result_time.pvalue, aspect='auto', origin='lower',cmap='RdBu_r')
plt.title('P-values of T-test between Head and Background conditions on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar(ticks=np.arange(0, 1, 0.05))
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/P-values_time.png')

# P-value zoomed
plt.figure(figsize=(8, 10))
plt.imshow(t_result_time.pvalue, aspect='auto', origin='lower',cmap='RdBu_r', vmin=0, vmax=0.1)
plt.title('P-values of T-test between Head and Background conditions on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/P-values_time-ZommedIN.png')


# plot t_result_channels
plt.figure()
plt.imshow(t_result_channels.statistic,extent=[np.min(tfrs_head[0].times), np.max(tfrs_head[0].times), np.min(tfrs_head[0].freqs), np.max(tfrs_head[0].freqs)], aspect='auto', origin='lower',cmap='RdBu_r')
plt.title('T-test between Head and Background conditions on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/T-test_channels.png')

# Plot p-values
plt.figure(figsize=(8, 6))
plt.imshow(t_result_channels.pvalue, extent=[np.min(tfrs_head[0].times), np.max(tfrs_head[0].times), np.min(tfrs_head[0].freqs), np.max(tfrs_head[0].freqs)], aspect='auto', origin='lower',cmap='RdBu_r')
plt.title('P-values of T-test between Head and Background conditions on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar(ticks=np.arange(0, 1, 0.05))
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/P-values_channels.png')

# P-value zoomed
plt.figure(figsize=(8, 6))
plt.imshow(t_result_channels.pvalue, extent=[np.min(tfrs_head[0].times), np.max(tfrs_head[0].times), np.min(tfrs_head[0].freqs), np.max(tfrs_head[0].freqs)], aspect='auto', origin='lower',cmap='RdBu_r', vmin=0, vmax=0.1)
plt.title('P-values of T-test between Head and Background conditions on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/P-values_channels-ZommedIN.png')


# plot t_result_freq
plt.figure(figsize=(8, 10))
plt.imshow(t_result_freq.statistic, aspect='auto', origin='lower',cmap='RdBu_r', vmin=-5, vmax=5)
plt.title('T-test between Head and Background conditions on frequency axis')
plt.xlabel('Time')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/T-test_freq.png')

# plot p-values
plt.figure(figsize=(8, 10))
plt.imshow(t_result_freq.pvalue, aspect='auto', origin='lower',cmap='RdBu_r')
plt.title('P-values of T-test between Head and Background conditions on frequency axis')
plt.xlabel('Time')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar(ticks=np.arange(0, 1, 0.05))
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/P-values_freq.png')


# T-test 

## Variances
## Raw difference
difference_raw = grAverage_head_data - grAverage_bgrd_data
overall_variance = np.var(difference_raw, axis=None)
overall_var_time = np.var(difference_raw, axis=2)
overall_var_freq = np.var(difference_raw, axis=1)
overall_var_channels = np.var(difference_raw, axis=0)

## Head condition
head_overall_variance = np.var(grAverage_head_data, axis=None)
head_overall_var_time = np.var(grAverage_head_data, axis=2)
head_overall_var_freq = np.var(grAverage_head_data, axis=1)
head_overall_var_channels = np.var(grAverage_head_data, axis=0)
## Background condition
bgrd_overall_variance = np.var(grAverage_bgrd_data, axis=None)
bgrd_overall_var_time = np.var(grAverage_bgrd_data, axis=2)
bgrd_overall_var_freq = np.var(grAverage_bgrd_data, axis=1)
bgrd_overall_var_channels = np.var(grAverage_bgrd_data, axis=0)

# Plot variance in overall data time
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(overall_var_time, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the difference between Head and Background conditions on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_time.png')

# plot variance in overall data frequency
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(overall_var_freq, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the difference between Head and Background conditions on frequency axis')
plt.xlabel('Time')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_freq.png')

# plot variance in overall data channels
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(overall_var_channels, extent=[np.min(tfrs_head[0].times), np.max(tfrs_head[0].times), np.min(tfrs_head[0].freqs), np.max(tfrs_head[0].freqs)],aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the difference between Head and Background conditions on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_channels.png')

# Plot variance in head data time
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(head_overall_var_time, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the Head condition on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_head_time.png')

# plot variance in head data frequency

plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(head_overall_var_freq, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the Head condition on frequency axis')
plt.xlabel('Time')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_head_freq.png')

# plot variance in head data channels
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(head_overall_var_channels,extent=[np.min(tfrs_head[0].times), np.max(tfrs_head[0].times), np.min(tfrs_head[0].freqs), np.max(tfrs_head[0].freqs)], aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the Head condition on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_head_channels.png')

# plot variance in background data time
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(bgrd_overall_var_time, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the Background condition on time axis')
plt.xlabel('Frequencies')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_bgrd_time.png')

# plot variance in background data frequency
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(bgrd_overall_var_freq, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the Background condition on frequency axis')
plt.xlabel('Time')
plt.ylabel('Channels')
plt.yticks(np.arange(0, len(tfrs_head[0].ch_names), 1), tfrs_head[0].ch_names)
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_bgrd_freq.png')

# plot variance in background data channels
plt.clf()
plt.figure(figsize=(8, 10))
plt.imshow(bgrd_overall_var_channels, aspect='auto', origin='lower', cmap='RdBu_r')
plt.title('Variance of the Background condition on channel axis')
plt.xlabel('Time')
plt.ylabel('Frequencies')
plt.colorbar()
plt.tight_layout()
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/Variance_bgrd_channels.png')



f_mne_score = f_oneway(data_head, data_bgrd)

## Variance within the group -- All participants

# All axes
t_result_head_all = scipy.stats.ttest_1samp(data_head, popmean=np.mean(data_head), axis=None)



#### TFCE ####
from mne.stats import combine_adjacency, permutation_cluster_test, permutation_cluster_1samp_test

# Set neighboring relations

spa_adjacency, ch_names = mne.channels.find_ch_adjacency(tfrs_head[0].info, ch_type='eeg')


# n_times = len(tfrs_head[0].times)
# temp_adjacency = mne.stats.time_adjacency(n_times)


# Reshape to (n_subjects, n_times, n_freqs, n_channels)
# data_head = np.transpose(data_head, (0, 3, 2, 1))
# data_bgrd = np.transpose(data_bgrd, (0, 3, 2, 1))

assert data_head.data.shape == (len(data_head), len(tfrs_head[0].ch_names), len(tfrs_head[0].freqs), len(tfrs_head[0].times))

adjacency = combine_adjacency(spa_adjacency, len(tfrs_head[0].freqs), len(tfrs_head[0].times))

assert (
    adjacency.shape[0]
    == adjacency.shape[1]
    == len(tfrs_head[0].ch_names) * len(tfrs_head[0].freqs) * len(tfrs_head[0].times)
)

# Combine the data for permutation testing
X = [data_head, data_bgrd]

# Set up parameters for TFCE
n_permutations = 1000
threshold = 0.2 #dict(start=0.4, step=0.2) #0.6 #0.05 #dict(start=0, step=0.2)

# pval = 0.05  # arbitrary
# dfn = 1 #n_conditions - 1  # degrees of freedom numerator
# dfd = 17 #n_observations - n_conditions  # degrees of freedom denominator
# threshold = scipy.stats.f.ppf(1 - pval, dfn=dfn, dfd=dfd)  # F distribution

# Run the permutation cluster test with TFCE
t_obs, clusters, cluster_p_values, H0 = permutation_cluster_test(
    X, n_permutations=n_permutations, threshold=threshold, tail=0, adjacency=adjacency, n_jobs=1, out_type='mask'
)

#cluster_p_values

def collect_cluster_pvalues(thresholds):

    acc_t_values = []
    acc_clusters = []
    acc_p_values = []
    H_zeros = []
    labels = []
    

    for threshold in thresholds:

        labels.append(str(threshold))

        t_obs, clusters, cluster_p_values, H0 = permutation_cluster_test(
            X, n_permutations=n_permutations, threshold=threshold, tail=0, adjacency=adjacency, n_jobs=1, out_type='mask'
        )

        print(f'Threshold: {threshold} \n Number of clusters: {len(clusters)} \n Cluster p-values: {cluster_p_values}')
        acc_p_values.append(cluster_p_values)
        acc_clusters.append(clusters)
        acc_t_values.append(t_obs)
        H_zeros.append(H0)

    x_positions = range(len(acc_p_values))
    
    # Plot individual points
    for i, p_values in enumerate(acc_p_values):
        plt.scatter([i] * len(p_values), p_values, label=labels[i], color='b')

    # Customize x-axis
    plt.xticks(x_positions, labels)
    plt.ylabel("Cluster p-values")     
    plt.title("Threshold values")

    plt.savefig(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFCE_cluster_pvalues_{labels}.png')

    return acc_t_values, acc_clusters, acc_p_values, H_zeros

thresholds = [0.0]#[0.2, 0.4, 0.6]
random.seed(19)
acc_t_values, acc_clusters, acc_p_values, H_zeros = collect_cluster_pvalues(thresholds)

x_positions = range(len(acc_p_values))
labels =  [0.0]#['0.2', '0.4', '0.6']

# Print the cluster p-values for each threshold value

# Print report
for threshold, p_value in zip(labels, acc_p_values):
    
    print(f"Cluster p-values for the threshold value {threshold}: {p_value}")


plt.figure(figsize=(8, 6))
# Plot individual points
for i, p_values in enumerate(acc_p_values):
    plt.scatter([i] * len(p_values), p_values, label=labels[i], color='#0343DF')

    # Customize x-axis
plt.xticks(x_positions, labels)

plt.ylabel("Cluster p-values") 
plt.xlabel("Threshold values")    
plt.title("Cluster p-values for different threshold values")
plt.ylim([0.9, 1])

plt.savefig(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFCE_cluster_pvalues_{str(thresholds)}.png')


# Plot the null distribution from H_zeros in subplots

fig, ax = plt.subplots(1, len(thresholds), figsize=(5*(len(thresholds)), 5))
for i, H0 in enumerate(H_zeros):

    if len(thresholds) > 1:
        ax[i].hist(H0, bins=50, color='gray', alpha=0.75, label='Null distribution')
        ax[i].set_title(f'Threshold: {labels[i]}')
        ax[i].set_xlabel('T-scores')
        ax[i].set_ylabel('Count of the clusters')
        ax[i].legend()
        ax[i].axvline(np.max(acc_t_values[i]), color='red', linestyle='--', label='Observed max T-value')
    else:
        ax.hist(H0, bins=50, color='gray', alpha=0.75, label='Null distribution')
        ax.set_title(f'Threshold: {labels[i]}')
        ax.set_xlabel('T-scores')
        ax.set_ylabel('Count of the clusters')
        ax.legend()
        ax.axvline(np.max(acc_t_values[i]), color='red', linestyle='--', label='Observed max T-value')


fig.savefig(f'/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/CBP_null_distribution_{str(thresholds)}.png')



plt.hist(H0, bins=50, color='gray', alpha=0.75, label='Null distribution')

# Add a vertical line for the maximum observed T-value
plt.axvline(np.max(t_obs), color='red', linestyle='--', label='Observed max T-value')

plt.xlabel('TFCE statistic (T-value)')
plt.ylabel('Count of the clusters')
plt.title('Null Distribution vs Observed Maximum T-value')
plt.legend()
plt.show()
#
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/CPT_statistic_0.05.png')


## Map the clusters on the spectrogram

# Subtract Background condition from Head condition



difference_tfr = grandAverage_bgrd
difference_tfr.data = grandAverage_head.data - grandAverage_bgrd.data

fig = difference_tfr.plot([0], baseline=(-0.5, -0.2), mode='logratio', title='Difference of activity between Head and Background conditions') 
fig[0].savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFR_difference.png')


picks = ['Oz']
difference_tfr_Oz = difference_tfr.copy().pick_channels(picks)


fig = difference_tfr_Oz.plot([0], baseline=(-0.5, -0.2), mode='logratio', title='Difference of activity at Oz')
fig[0].savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFR_difference_Oz.png')

plt.figure()
plt.imshow(difference_tfr, cmap='RdBu_r')
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFCE_difference.png')

cmap = matplotlib.colors.ListedColormap(['none', 'green'])

plt.imshow(clusters[0], cmap=cmap, alpha=1)
plt.title('Clusters mapped on the spectogram of difference between Head and Background conditions')
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFCE_clusters_mapped.png')



plt.figure()
plt.imshow(data_head.data[0], aspect='auto', origin='lower',
           extent=[tfrs_head[0].times[0], tfrs_head[0].times[-1], freqs[0], freqs[-1]],
           cmap='RdBu_r', interpolation='bilinear')
plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFCE_clusters_.png')


## TFCE at Oz
tfrs_head_Oz = [tfr.copy().pick_channels(['Oz']) for tfr in tfrs_head]
tfrs_bgrd_Oz = [tfr.copy().pick_channels(['Oz']) for tfr in tfrs_bgrd]

data_head_Oz = np.array([tfr.data for tfr in tfrs_head_Oz])
data_bgrd_Oz = np.array([tfr.data for tfr in tfrs_bgrd_Oz])

grandAverage_head_Oz = grandAverage_head.copy().pick_channels(['Oz'])
grandAverage_bgrd_Oz = grandAverage_bgrd.copy().pick_channels(['Oz'])

X_Oz = [data_head_Oz, data_bgrd_Oz]

random.seed(19)
n_permutations = 1000
threshold = 5 
t_obs_Oz, clusters_Oz, cluster_p_values_Oz, H0_Oz = permutation_cluster_test(
    X_Oz, n_permutations=n_permutations, threshold=None, tail=0, adjacency=None, n_jobs=1, out_type='mask'
)
plt.figure()
plt.hist(H0_Oz, bins=50, color='gray', alpha=0.75, label='Null distribution', range=(0, 10*np.max(t_obs_Oz)))

# Add a vertical line for the maximum observed T-value
plt.axvline(np.max(t_obs_Oz), color='red', linestyle='--', label='Observed max T-value')

plt.xlabel('TFCE statistic (T-value)')
plt.ylabel('Count of the clusters')
plt.title('Null Distribution vs Observed Maximum T-value')
plt.legend()
plt.show()
#

plt.savefig('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots/MNEPlots/TFCE_statistic_0.05_Oz.png')


### ------------------------- ### ________________ ### ------------------------- ###
# baselined_weighted_head_tfr = weighted_head_tfr.apply_baseline((-0.5, -0.2), mode='logratio')
# baselined_weighted_bgrd_tfr = weighted_bgrd_tfr.apply_baseline((-0.5, -0.2), mode='logratio')

# head_tfr = baselined_weighted_head_tfr.data[:, :, :]
# bgrd_tfr = baselined_weighted_bgrd_tfr.data[:, :, :]

# threshold_tfce = dict(start=0, step=0.2)

# f_scores, clusters, p_values, H0 = permutation_cluster_1samp_test(
#     [head_tfr, bgrd_tfr],
#     n_jobs=None,
#     threshold=threshold_tfce,
#     adjacency=None,
#     n_permutations=1000,
#     out_type="mask",
# )

