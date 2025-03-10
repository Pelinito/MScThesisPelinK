% The subjects are brought together in corresponding conditions under one big subject to keep all trials for that condition

clear all;
close all;
clc;

%% Define the paths

tfce_path = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE';
plots_tfce = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots';
savesets_here = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets';
analysis_results = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/Analysis';

addpath(savesets_here);
addpath(analysis_results);

addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');
cd(tfce_path)

%% Load the .mat files
cd(savesets_here);
save('allSubjallCond.mat', 'allSubjallCond', '-v7.3');
load("allSubjallCond.mat");

%% Prepare data for the "surrogate" subject 
% surrogateAll: Stores the data from all subjects in both conditions
% surrogateHead: Stores the data from all subjects in the Head condition
% surrogateBground: Stores the data from all subjects in the Bground
% condition

% SurrogateAll
surrogateAll = allSubjallCond{1,1};
NSubj = size(allSubjallCond,2);

for subj = 1:(NSubj-1)
    size1 = size(surrogateAll.trial,2);
    size2 = size(allSubjallCond{1,subj+1}.trial,2);
    surrogateAll.trial(size1+1:size1+size2) = allSubjallCond{1,subj+1}.trial;
end

surrogateAll.hdr.nTrials = size(surrogateAll.trial,2);
surrogateAll.time(1:size1+size2) = surrogateAll.time(1,1);


NTrials = size(surrogateAll.trial,2);

for trial = 1:NTrials-1
    surrogateAll.sampleinfo(trial+1, 1) = surrogateAll.sampleinfo(trial,1) + 1050;
    surrogateAll.sampleinfo(trial+1, 2) = surrogateAll.sampleinfo(trial,2) + 1050;
end

% SurrogateHead

surrogateHead = allSubjHead{1,1};
NSubj = size(allSubjHead,2);

for subj = 1:(NSubj-1)
    size1 = size(surrogateHead.trial,2);
    size2 = size(allSubjHead{1,subj+1}.trial,2);
    surrogateHead.trial(size1+1:size1+size2) = allSubjHead{1,subj+1}.trial;
end

surrogateHead.hdr.nTrials = size(surrogateHead.trial,2);
surrogateHead.time(1:size1+size2) = surrogateHead.time(1,1);


NTrials = size(surrogateHead.trial,2);

for trial = 1:NTrials-1
    surrogateHead.sampleinfo(trial+1, 1) = surrogateHead.sampleinfo(trial,1) + 1050;
    surrogateHead.sampleinfo(trial+1, 2) = surrogateHead.sampleinfo(trial,2) + 1050;
end

% Surrogate Bground

surrogateBground = allSubjBgrd{1,1};
NSubj = size(allSubjBgrd,2);

for subj = 1:(NSubj-1)
    size1 = size(surrogateBground.trial,2);
    size2 = size(allSubjBgrd{1,subj+1}.trial,2);
    surrogateBground.trial(size1+1:size1+size2) = allSubjBgrd{1,subj+1}.trial;
end

surrogateBground.hdr.nTrials = size(surrogateBground.trial,2);
surrogateBground.time(1:size1+size2) = surrogateBground.time(1,1);


NTrials = size(surrogateBground.trial,2);

for trial = 1:NTrials-1
    surrogateBground.sampleinfo(trial+1, 1) = surrogateBground.sampleinfo(trial,1) + 1050;
    surrogateBground.sampleinfo(trial+1, 2) = surrogateBground.sampleinfo(trial,2) + 1050;
end
%% Channel(s)
% If multiple channels are desired, input as a cell array, e.g {'Oz','POz'}
channel_selection = 'Oz';
%% Change directory to save the results of the analyses
cd(analysis_results);
%% Time Analysis Preparation


cfg_time = [];
cfg_time.channel = channel_selection;
cfg_time.trials = 'all';

time_surrogateAll = ft_timelockanalysis(cfg_time, surrogateAll);
time_surrogateHead = ft_timelockanalysis(cfg_time, surrogateHead);
time_surrogateBground = ft_timelockanalysis(cfg_time, surrogateBground);

%% Plot

figure;
hold on
title(sprintf('ERP of all trials of all subjects at channel(s): %s', channel_selection));
xlabel('Time[s]');
ylabel(sprintf('Potential[%s]', time_surrogateAll.elec.chanunit{1}));
xline(0)

plot(time_surrogateAll.time, time_surrogateAll.avg)
plot(time_surrogateHead.time, time_surrogateHead.avg)
plot(time_surrogateBground.time, time_surrogateBground.avg)
legend('Both Conditions', 'Head', 'Background');
%% Time-Frequency Analysis

cfg_tf = [];
cfg_tf.method = 'wavelet';
cfg_tf.output = 'pow';

cfg_tf.toi =  linspace(min(allSubjallCond{1,1}.time{1,1}),max(allSubjHead{1,1}.time{1,1}), length(allSubjHead{1,1}.time{1,1}));
cfg_tf.foi = 1:0.5:45;
cfg_tf.channel = channel_selection;
cfg_tf.pad = 'nextpow2';


cfg_tf.width = 1;
tf_surrogateAll_cyc1 = ft_freqanalysis(cfg_tf, surrogateAll);
save('tf_surrogateAll_cyc1.mat', 'tf_surrogateAll_cyc1');
tf_surrogateHead_cyc1= ft_freqanalysis(cfg_tf, surrogateHead);
save('tf_surrogateHead_cyc1.mat', 'tf_surrogateHead_cyc1');
tf_surrogateBground_cyc1 = ft_freqanalysis(cfg_tf, surrogateBground);
save('tf_surrogateBground_cyc1.mat', 'tf_surrogateBground_cyc1');

cfg_tf.width = 3;
tf_surrogateAll_cyc3 = ft_freqanalysis(cfg_tf, surrogateAll);
save('tf_surrogateAll_cyc3.mat', 'tf_surrogateAll_cyc3');
tf_surrogateHead_cyc3 = ft_freqanalysis(cfg_tf, surrogateHead);
save('tf_surrogateHead_cyc3.mat', 'tf_surrogateHead_cyc3');
tf_surrogateBground_cyc3 = ft_freqanalysis(cfg_tf, surrogateBground);
save('tf_surrogateBground_cyc3.mat', 'tf_surrogateBground_cyc3');

cfg_tf.width = 5;
tf_surrogateAll_cyc5 = ft_freqanalysis(cfg_tf, surrogateAll);
save('tf_surrogateAll_cyc5.mat', 'tf_surrogateAll_cyc5');
tf_surrogateHead_cyc5 = ft_freqanalysis(cfg_tf, surrogateHead);
save('tf_surrogateHead_cyc5.mat', 'tf_surrogateHead_cyc5');
tf_surrogateBground_cyc5 = ft_freqanalysis(cfg_tf, surrogateBground);
save('tf_surrogateBground_cyc5.mat', 'tf_surrogateBground_cyc5');
%% Load TF
cd(analysis_results);

load('tf_surrogateAll_cyc1.mat');
load('tf_surrogateHead_cyc1.mat');

%% Plot TF

cfg_plotTF = [];
cfg_plotTF.parameter = 'powspctrm';
cfg_plotTF.baselinetype = 'db';
cfg_plotTF.baseline = [-0.5 -0.2];
cfg_plotTF.colormap = 'RdBu';
% cfg_plotTF.title = sprintf('ERSPs for the Head condition at the channel %s using %d cycles',cfg_tf.channel, cfg_tf.width);
% cfg_plotTF.zlim = [0 9000];

cfg_plotTF.title = sprintf('ERSPs for the Background condition at the channel %s using %d cycles',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, tf_surrogateBground_cyc1);

cfg_plotTF.title = sprintf('ERSPs for the Head condition at the channel %s using %d cycles',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, tf_surrogateHead_cyc1);

cfg_plotTF.title = sprintf('ERSPs for all conditions at the channel %s using %d cycles',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, tf_surrogateAll_cyc1);
