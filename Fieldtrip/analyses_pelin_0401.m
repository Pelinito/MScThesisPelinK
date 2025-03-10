% Attempting TFCE in Fieltrip

clear all;
close all;
clc;

%% Define the paths

tfce_path = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE';
plots_tfce = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots';
savesets_here = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets';
analysis_results = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/Analysis';
plots_each_subject = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/PlotsEachSubject';

addpath(savesets_here);
addpath(analysis_results);

addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');
cd(tfce_path)

%% Load the .mat files
cd(savesets_here);

load("allSubjHead.mat");
load("allSubjBgrd.mat");

%% Prep for Analyses

NSubj = 18;

% Empty cell arrays for time and time-freq analyses
% All Subjects
time_head = cell(1, NSubj);
time_bground = cell(1, NSubj);
tf_head = cell(1, NSubj);
tf_bground = cell(1, NSubj);

% One Subject

% oneSubj_time_head = cell(1,1);
% oneSubj_time_bground = cell(1,1);
% oneSubj_tf_head = cell(1,1);
% oneSubj_tf_bground = cell(1,1);

%% Time Analysis

% Define the channel(s) here
channel_selection = 'all';

% PREP

cfg_time = [];
cfg_time.channel = channel_selection;
cfg_time.trials = 'all';
%cfg_time.keeptrials = 'yes';

% 2 conditions
for subj = 1:NSubj
    
    time_head{1,subj} = ft_timelockanalysis(cfg_time, allSubjHead{1,subj});
    time_bground{1,subj} = ft_timelockanalysis(cfg_time, allSubjBgrd{1,subj});
end

% % All in one
% 
% for subj =1:NSubj*2
%     
%     time_all{1,subj} = ft_timelockanalysis(cfg_time, allSubjallCond{1,subj});
%     
% end
cd(analysis_results)
save('time_head_allchan.mat', 'time_head');
save('time_bground_allchan.mat', 'time_bground');
%% Plot the potential for each subject in each condition

cd(plots_each_subject);

for subj = 1:NSubj
   
    plotname = sprintf('ERP_head_Subj_%d.png',subj);
    
    ax = figure;
    hold on
    title(sprintf('ERP of subject %d for head condition averaged across all channels', subj));
    xlabel('Time[s]');
    ylabel(sprintf('Potential[%s]', time_head{1,subj}.elec.chanunit{1}));
    xline(0)
    
    plot(time_head{1,subj}.time, mean(time_head{1,subj}.avg))
    
    %ax = gca;
    saveas(ax, plotname);
end

for subj = 1:NSubj
   
    plotname = sprintf('ERP_bground_Subj_%d.png',subj);
    
    ax = figure;
    hold on
    title(sprintf('ERP of subject %d for background condition averaged across all channels', subj));
    xlabel('Time[s]');
    ylabel(sprintf('Potential[%s]', time_bground{1,subj}.elec.chanunit{1}));
    xline(0)
    
    plot(time_bground{1,subj}.time, mean(time_bground{1,subj}.avg))
    
    %ax = gca;
    saveas(ax, plotname);
end


% Subplot for each

for subj = 1:NSubj
   
    plotname = sprintf('ERP_together_Subj_%d.png',subj);
    
    ax = figure;
    subplot(1,2,1)
    plot(time_head{1,subj}.time, mean(time_head{1,subj}.avg))
    title(sprintf('ERP of subject %d for head condition averaged across all channels', subj));
    xlabel('Time[s]');
    ylabel(sprintf('Potential[%s]', time_head{1,subj}.elec.chanunit{1}));
    xline(0)
    ylim([-0.2 0.2])
    
    subplot(1,2,2)
    plot(time_bground{1,subj}.time, mean(time_bground{1,subj}.avg))
    title(sprintf('ERP of subject %d for background condition averaged across all channels', subj));
    xlabel('Time[s]');
    ylabel(sprintf('Potential[%s]', time_bground{1,subj}.elec.chanunit{1}));
    xline(0)
    ylim([-0.2 0.2])

    set(ax, 'Position', [100 100 1500 600]');
    saveas(ax, plotname);
end


%% Plot the potenatial across the trials for each subject
% HEAD
figure;
hold on
title(sprintf('Potential of each subject for head condition at channel(s): %s', channel_selection));
xlabel('Time[s]');
ylabel(sprintf('Potential[%s]', time_head{1,1}.elec.chanunit{1}));
xline(0)

for subj = 1:NSubj
    plot(time_head{1,1}.time, time_head{1,subj}.avg)
end

% HEAD - with and w/o baseline

figure;
%hold on
title(sprintf('Potential of each subject for head condition at channel(s): %s', channel_selection));
xlabel('Time[s]');
ylabel(sprintf('Potential[%s]', time_head{1,1}.elec.chanunit{1}));
xline(0)

plot(time_head{1,1}.time, time_head{1,subj}.avg)
hold on
plot(bl_time_head{1,1}.time, bl_time_head{1,subj}.avg)
legend('without baseline', 'with baseline')

% BACKGROUND
figure;
hold on
title(sprintf('Potential of each subject for background condition at channel(s): %s', channel_selection));
xlabel('Time[s]');
ylabel(sprintf('Potential[%s]', time_bground{1,1}.elec.chanunit{1}));
xline(0)

for subj = 1:NSubj
    plot(time_bground{1,1}.time, time_bground{1,subj}.avg)
end

% The head condition is noisy. Take the mean across the subjects

avg_time_head = mean(cell2mat((cellfun(@(s) s.avg, time_head, 'UniformOutput', false)')),1);
avg_bl_time_head = mean(cell2mat((cellfun(@(s) s.avg, bl_time_head, 'UniformOutput', false)')),1);
% Explaining the obnoxious line above:
% (1) cellfun takes the .avg field of from the structs that keeps the
% activation of each subject
% (2) Cellfun results in a cell array of 1x18, transpose it to 18x1
% (3) Convert the cell array to a matrix with cell2mat
% (4) Calculate the mean across the rows

% Plot the average of all subjects

figure;
hold on
title(sprintf('Average activation across all subjects for the head condition at the channel(s): %s', channel_selection));
ylabel(sprintf('Potential[%s]', time_head{1,1}.elec.chanunit{1}));
xlabel('Time[s]');
xline(0);
plot(time_head{1,1}.time, avg_time_head);
plot(time_head{1,1}.time, avg_bl_time_head);
legend('without baseline', 'with baseline')

% Might do the same for the background, doesn't hurt anyone

avg_time_bground = mean(cell2mat((cellfun(@(s) s.avg, time_bground, 'UniformOutput', false)')),1);
avg_bl_time_bground = mean(cell2mat((cellfun(@(s) s.avg, bl_time_bground, 'UniformOutput', false)')),1);

figure;
hold on
title(sprintf('Average activation across all subjects for the background condition at the channel(s): %s', channel_selection));
ylabel(sprintf('Potential[%s]', time_bground{1,1}.elec.chanunit{1}));
xlabel('Time[s]');
xline(0);
plot(time_bground{1,1}.time, avg_time_bground);
plot(time_bground{1,1}.time, avg_bl_time_bground);
legend('without baseline', 'with baseline')

% All in one

% figure;
% hold on
% title(sprintf('Potential of each subject for background condition at channel(s): %s', channel_selection));
% xlabel('Time[s]');
% ylabel(sprintf('Potential[%s]', time_all{1,1}.elec.chanunit{1}));
% xline(0)

avg_time_all = mean(cell2mat((cellfun(@(s) s.avg, time_all, 'UniformOutput', false)')),1);
avg_bl_time_all = mean(cell2mat((cellfun(@(s) s.avg, bl_time_all, 'UniformOutput', false)')),1);


figure;
hold on
title(sprintf('Average activation across all subjects for the background condition at the channel(s): %s', channel_selection));
ylabel(sprintf('Potential[%s]', time_all{1,1}.elec.chanunit{1}));
xlabel('Time[s]');
xline(0);
plot(time_bground{1,1}.time, avg_time_all);
%plot(time_bground{1,1}.time, avg_bl_time_bground);
legend('without baseline', 'with baseline')

%% TF Analysis

%channel_selection = 'Oz';
% widths = [3 5 7]; % Maybe use this list to loop in it for subplotting.
% Idea is not finalized!

cfg_tf = [];
cfg_tf.method = 'wavelet';
cfg_tf.output = 'pow';
cfg_tf.width = 3;
cfg_tf.toi =  linspace(min(allSubjHead{1,1}.time{1,1}),max(allSubjHead{1,1}.time{1,1}), length(allSubjHead{1,1}.time{1,1}));
cfg_tf.foi = 1:0.5:45;
cfg_tf.channel = channel_selection;
cfg_tf.pad = 'nextpow2';
%cfg_tf.keeptrials = 'yes';       %---------DO I NEED THIS?


for subj = 1:NSubj
    
    tf_head{1,subj} = ft_freqanalysis(cfg_tf, allSubjHead{1,subj});
    tf_bground{1,subj} = ft_freqanalysis(cfg_tf, allSubjBgrd{1,subj});
end

for subj = 1:NSubj*2
    
    sprintf('Subject %d/%d', subj,NSubj*2);
    tf_all{1,subj} = ft_freqanalysis(cfg_tf, allSubjallCond{1,subj});
    
end

%% Apply baseline
cfg_baseline.baselinetype = 'db';
for subj = 1:NSubj
    
    bl_tf_head{1,subj} = ft_freqbaseline(cfg_baseline, tf_head{1,subj});
    bl_tf_bground{1,subj} = ft_freqbaseline(cfg_baseline, tf_bground{1,subj});
end

%% Grandaverage -- BASELINED

concatenated_powspec_head = bl_tf_head{1,subj}.powspctrm;
concatenated_powspec_bground = bl_tf_bground{1,subj}.powspctrm;
for subj = 1:NSubj-1
    
    concatenated_powspec_head = cat(1, concatenated_powspec_head, tf_head{1,subj+1}.powspctrm);
    concatenated_powspec_bground = cat(1, concatenated_powspec_bground, tf_bground{1,subj+1}.powspctrm);
end

avg_tf_head = mean(concatenated_powspec_head, 1);
avg_tf_bground = mean(concatenated_powspec_bground, 1);

% avg_time_head = mean(cell2mat((cellfun(@(s) s.avg, time_head, 'UniformOutput', false)')),1);
% avg_bl_time_head = mean(cell2mat((cellfun(@(s) s.avg, bl_time_head, 'UniformOutput', false)')),1);

%% Surrogate struct for the averaged data

avg_tf_head_st = tf_head{1,1};
avg_tf_head_st.powspctrm = avg_tf_head;

avg_tf_bground_st = tf_bground{1,1};
avg_tf_bground_st.powspctrm = avg_tf_bground;
%% Plot the spectrogram

cfg_plotTF = [];
cfg_plotTF.parameter = 'powspctrm';
% cfg_plotTF.baselinetype = 'db';
% cfg_plotTF.baseline = [-500 -200];
cfg_plotTF.colormap = 'RdBu';
cfg_plotTF.title = sprintf('ERSPs for the Head condition at the channel %s using %d cycles',cfg_tf.channel, cfg_tf.width);
%cfg_plotTF.zlim = [0 9000];
ft_singleplotTFR(cfg_plotTF, tf_head{1,1});

cfg_plotTF.title = sprintf('ERSPs for the Head condition at the channel %s using %d cycles with baseline',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, bl_tf_head{1,1});

cfg_plotTF.title = sprintf('ERSPs for the Head condition at the channel %s using %d cycles averaged across subjects',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, avg_tf_head_st);

cfg_plotTF.title = sprintf('ERSPs for the Background condition at the channel %s using %d cycles',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, tf_bground{1,1});

cfg_plotTF.title = sprintf('ERSPs for the Background condition at the channel %s using %d cycles with baseline',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, bl_tf_bground{1,1});

cfg_plotTF.title = sprintf('ERSPs for the Background condition at the channel %s using %d cycles averaged across subjects',cfg_tf.channel, cfg_tf.width);
ft_singleplotTFR(cfg_plotTF, avg_tf_bground_st);