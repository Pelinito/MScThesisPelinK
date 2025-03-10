%% Goal: Apply TFCE to ONE subject for HEAD vs. NON-HEAD conditions on TIME DOMAIN - 25.09.2023
% Steps: 
% - Create the pathways
% - Load the subject data for one subject - DONE AND NOW DOING MULTIPLE
% - Seperate the events
% - Apply time analysis with FieldTrip
% - Apply TFCE

%% Progress: 04.10.2023
% Created the design
% Prepared the neighbours
% Extended to 3 subjects
% Plotted the TFCE for 3 subjects

%% Goal: 05.10.2023
%Fix the code so it will not give errors about the cell output after
%looping one subject
% Save the transformed data for all conditions - This takes the longest,
% get it done for ease in the future
%
% Fix the error: "Expected one output from a curly brace or dot indexing
% expression, but there were 68 results."

%% Progress: 06.10.2023
% Fixed the above-mentioned error
% The data is saved under respective file names
% Now exists: allSubjHead.mat, allSubjBody.mat, allSubjBgrd.mat,
% allSubjnonHead.mat.
% The timelock analysis and timelock statistics are conducted.

%% Further goals: From 11.10.2023
% This script is going to be sectioned into shorter scripts to only conduct
% what is necessary to avoid the clutter in the code.
% Refer further to the script, tfce_v1.

%% PROBLEM - realized on 11.10.2023

% The .mat files for allSubjBgrd.mat, allSubjHead.mat, allSubjnonHead.mat
% are not saved. They are empty. 

% The script should be investigated. 
% Write an additional script to preprocess all those data. 
% Now, you do not need any file outside of the .../Pelin/TFCE folder.
% The sets are epoched and seperated are saved under
% .../Pelin/TFCE/ConditionSets
% Load 'em, loop 'em, save 'em.
% Then, continue with analysis and statistics.


%% ------------------------------------------------------------------------

%% 
% clear variables
% close all;
% clc;

%% 
% 
% savepath = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE';
% save_plots = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots';
% 
% cd(savepath);
% 

% addpath('/net/store/nbp/projects/wd_ride_village/processedData/village');
% addpath('/net/store/nbp/projects/wd_ride_village/processedData/village/preprocessed');
% addpath('/net/store/nbp/projects/wd_ride_village/Analysis');
% addpath(save_plots);
% 
% % load EEGlab
% addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
% eeglab;
% % load Fieldtrip
% addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');

%% 
%% to restart everything
% % do not run when running it automatically
clear variables
close all;
clc;
%% Load necessary code
basepath='/net/store/nbp/projects/wd_ride_village/recordedData/wd_village/';
cd(basepath);

savepath = '/net/store/nbp/projects/wd_ride_village/processedData/village';
addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE');
addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets')
addpath('/net/store/nbp/projects/wd_ride_village/processedData/village');
addpath('/net/store/nbp/projects/wd_ride_village/Analysis');

tfce_path = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE';
%transformedtoFT = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets/TransformedtoFT';

% load EEGlab
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
eeglab;
% load Fieldtrip
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');
ft_defaults

savedata_all = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE';
savesets_here = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets';

plotting = false;

% Directory to save the plots
plots_tfce = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots';

if exist(plots_tfce, 'dir')
    addpath(genpath(plots_tfce));
else
    mkdir /net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots;
    addpath(genpath(plots_tfce));
end

% %% Test the recursive dir
% 
% rootdir = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK/';
% 
% %filelist = dir(fullfile(rootdir, '**/*.*'));  %get list of files and folders in any subfolder
% len_channels_test = length(filelist);
% %filelist = filelist(~[filelist.isdir]);  %remove folders from list


%% get subjects to include
%subjects = [1,2,5];%12,15,16,17,19,20,21,22,27,29,32,33,34,37,38];
subjects = [1,2,5,12,15,16,17,19,20,21,22,27,29,32,33,34,37,38];
%% variables to save the data for all subjects
% variables to store the data for all subjects

%data_all = []; ---- NOT NEEDED for ONE (1) SUBJECT

% Create the condition tables
%Nsubj = length(subjects);

% allSubjx cell structures to save the data transformed into FT, appends in
% the loop

allSubjHead = {};
allSubjBody = {};
allSubjBgrd = {};
allSubjnonHead = {};


%% for loop: load individual subjects
for s = 1:length(subjects)

cd(basepath);
  
%sub = subjects(s);
if subjects(s)<10
    subStr = '0'+string(subjects(s));
else
    subStr = string(subjects(s));
end

% get the data path
if subjects(s) < 40
    cntfile = dir(subStr+'_v_'+'*'+'.xdf');
else
    cntfile = dir(subStr+'_v1_'+'*'+'.xdf');
end
cntpath = fullfile(basepath,cntfile.name);


% get the proper file name
if exist(cntpath,'file')
    [filepath,filename,ext] = fileparts(cntpath); % define variables
    % to keep the same as the manual selection procedure
    filepath = [filepath filesep];
    filename = [filename ext];
else
    cntpath = fullfile(basepath,sprintf('VP%u',subjects(s)),'*.*');
    [filename, filepath]=uigetfile(cntpath);
end

% load the csv file with the participant overview
rec_vill = readtable('recordings_village_corr.csv');
idx = -1;
% find the corresponding uid: firt get index 
for i = 1:height(rec_vill(:,2))
    if strcmp(rec_vill{i,2},filename)==1 
        idx = i;
    end
end

% now get the corresponding name
if  idx >= 0
    uidname = rec_vill{idx,1};
else
    %throw error if the uid could not be found
    error('No uid found')
end

%convert uidname and cut the ending 
uidname = uidname{1,1};

% create a target folder & change permissions so everyone
% in nbp can create and delete files
str = [savepath, '/preprocessed/', uidname];
savedata = join(str,''); % create new path for each individual subject  

% clear all variables that are no longer needed
%clear savepath cntfile ext filename i idx rec_vill str sub subStr
clc;

%% set cfg
cfg         = [];
cfg.avgref  = 1; % set it to 1 to rereference the data to avg ref
cfg.subj    = char(uidname); % get the name of the current subject 

cfg.mainfolder  = sprintf(savedata); %subject folder
cfg.ReData    = sprintf('new_full_data_%s.set',cfg.subj); %get the correct file of the subject


savesets_here = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets';
%%
EEG     = pop_loadset(cfg.ReData,cfg.mainfolder);
EEG     = eeg_checkset(EEG);


%% Change the type in the events ------------------------------- last part of Debbies code
%code taken from ERSP_analysis.m
%2 = 'bgrd'; 1 = 'body'; 0 = 'face'
head_count = 0;
body_count = 0;
bgrd_count = 0;

for evt = 1:length(EEG.event)
    if EEG.event(evt).type     == '0'
        EEG.event(evt).type     = 'head';
        head_count = head_count + 1;
    elseif EEG.event(evt).type == '1'
        EEG.event(evt).type     ='body';
        body_count = body_count + 1;
    elseif EEG.event(evt).type == '2'
        EEG.event(evt).type     ='bgrd';
        bgrd_count = bgrd_count + 1;
    end
end
disp('changing type done')



%% Seperate the events here

%events = {'head', 'body', 'bgrd'};
%EEG_epoched = pop_epoch(EEG, events, [-0.875 1.175], 'epochinfo', 'yes');

%% code from ERSP_analysis.m

event_codes = {'head', 'body', 'bgrd'};
%epochs all data of all conditions
EEG_epoched = pop_epoch(EEG, event_codes, [-0.875 1.175], 'epochinfo', 'yes');

%creating data file for face condition with epoched data
condition_event = find(strcmp({EEG.event.type}, 'head'));
EEG_condition_head = EEG;
EEG_condition_head.event = EEG_condition_head.event(condition_event);

EEG_epoched_head = pop_epoch(EEG_condition_head, event_codes, [-0.875 1.175], 'epochinfo', 'yes');
EEG_epoched_head = eeg_checkset(EEG_epoched_head);


%creating data file for body condition with epoched data
condition_event = find(strcmp({EEG.event.type}, 'body'));
EEG_condition_body = EEG;
EEG_condition_body.event = EEG_condition_body.event(condition_event);

EEG_epoched_body = pop_epoch(EEG_condition_body, event_codes, [-0.875 1.175], 'epochinfo', 'yes');
EEG_epoched_body = eeg_checkset(EEG_epoched_body);

%creating data file for background condition with epoched data
condition_event = find(strcmp({EEG.event.type}, 'bgrd'));
EEG_condition_bgrd = EEG;
EEG_condition_bgrd.event = EEG_condition_bgrd.event(condition_event);

EEG_epoched_bgrd = pop_epoch(EEG_condition_bgrd, event_codes, [-0.875 1.175], 'epochinfo', 'yes');
EEG_epoched_bgrd = eeg_checkset(EEG_epoched_bgrd);


%% create one file with all subject data while looping
% in the first loop the data_all is still empty 
% we go into the first statement and also save the data for subj 1
% if isempty(data_all)
%     data_all = EEG_epoched;    
% else
%     data_all = pop_mergeset(data_all, EEG_epoched);
% end


%% Save the datasets for each epoched condition - this saves .fdt files too, 
%that's why save fnc from matlab itself does not work
cd(savesets_here);

head_var = sprintf('EEG_epoched_head_%s.set',cfg.subj);
body_var = sprintf('EEG_epoched_body_%s.set',cfg.subj);
bgrd_var = sprintf('EEG_epoched_bgrd_%s.set',cfg.subj);

pop_saveset(EEG_epoched_head, head_var, savesets_here);
pop_saveset(EEG_epoched_body, body_var, savesets_here);
pop_saveset(EEG_epoched_bgrd, bgrd_var, savesets_here);

%% Load aka Preprocess the data for FieldTrip

cfg = [];

% load head condition
filename = head_var;
cfg.dataset = [savesets_here,'/',filename];
EEG_FT_head = ft_preprocessing(cfg);

% load body condition
filename = body_var;
cfg.dataset = [savesets_here,'/',filename];
EEG_FT_body = ft_preprocessing(cfg);

% load bgrd condition
filename = bgrd_var;
cfg.dataset = [savesets_here,'/',filename];
EEG_FT_bgrd = ft_preprocessing(cfg);

%% Merge the body and bgrd conditions to have the NON-HEAD condition

cfg = [];
cfg.keepsampleinfo = 'yes';

EEG_FT_nonHead = ft_appenddata(cfg, EEG_FT_body, EEG_FT_bgrd);

%% Save each data to the respective cell structure

allSubjHead{1,s} = EEG_FT_head;
allSubjBody{1,s} = EEG_FT_body;
allSubjBgrd{1,s} = EEG_FT_bgrd  ;
allSubjnonHead{1,s} = EEG_FT_nonHead ;



end

%% Save the cell structures as .mat to the folder
save('allSubjHead.mat', 'allSubjHead');
save('allSubjBody.mat', 'allSubjBody');
save('allSubjBgrd.mat', 'allSubjBgrd');
save('allSubjnonHead.mat', 'allSubjnonHead');

%% Time-lock analysis - Make a loop (05.10)

head_timelock = {};
nonHead_timelock = {};
body_timelock = {};
bgrd_timelock = {};



cfg = [];
cfg.channel = 'all';
%cfg.keeptrials 
cfg.trials = 'all';

for timelock=1:length(allSubjHead)
    
    head_timelock{1,timelock} = ft_timelockanalysis(cfg, allSubjHead{1,timelock});
end

for timelock=1:length(allSubjnonHead)
    
    nonHead_timelock{1,timelock} = ft_timelockanalysis(cfg, allSubjnonHead{1,timelock});
end

for timelock=1:length(allSubjBody)
    
    body_timelock{1,timelock} = ft_timelockanalysis(cfg, allSubjBody{1,timelock});
end

for timelock=1:length(allSubjBgrd)
    
    bgrd_timelock{1,timelock} = ft_timelockanalysis(cfg, allSubjBgrd{1,timelock});
end

% head_timelock = ft_timelockanalysis(cfg, EEG_FT_head);
% nonHead_timelock = ft_timelockanalysis(cfg, EEG_FT_nonHead);

%% Wavelets - for whatever reason - actually, to understand the example on FieldTrip website - Make a loop (05.10)

% cfg = [];
% cfg.channel    = 'all';
% cfg.method     = 'wavelet';
% cfg.output     = 'pow';
% 
% %cfg.gwidth     = 3;
% cfg.toi        = linspace(min(EEG_FT_head.time{1,1}),max(EEG_FT_head.time{1,1}), length(EEG_FT_head.time{1,1})); %matching data points in time depending on min and max in time in equal length
% %cfg.foilim     = [1 45];
% cfg.foi        = 1:0.5:45;
% cfg.pad = 'nextpow2';
% 
% % 3 cycles
% 
% cfg.width = 3;
% 
% head_tfr = ft_freqanalysis(cfg, EEG_FT_head);
% nonHead_tfr = ft_freqanalysis(cfg, EEG_FT_nonHead);


%% ----- TFCE -----
%% Design and Neighbours

% This is taken from the tfce tutorial on FieldTrip

Nsubj = 18;

design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];    % this is the uvar
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];  % this is the ivar

cfg = [];
cfg.method = 'distance';
cfg.neighbourdist = 40;

nb = ft_prepare_neighbours(cfg, head_timelock{1});


%% TFCE

cfg = [];
cfg.channel = 'all';
cfg.method = 'montecarlo';
cfg.correctm = 'tfce';
%cfg.frequency = [1 45];
cfg.avgoverchannel = 'no';
cfg.avgovertime = 'no';
%cfg.avgoverfreq = 'no';
cfg.clusteralpha = 0.05;
cfg.tfce_H = 2;
cfg.tfce_E = 0.5;
cfg.statistic = 'indepsamplesT';
%cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;   % 2^Nsubj possible randomizations

cfg.design = design;
%cfg.uvar = 1;
cfg.ivar = 2;

%stat_tfce = ft_freqstatistics(cfg, head_tfr, nonHead_tfr);
stat_tfce_timelock = ft_timelockstatistics(cfg, head_timelock{:}, nonHead_timelock{:});

stat_tfce_timelock2 = ft_timelockstatistics(cfg, nonHead_timelock{:}, head_timelock{:});

stat_tfce_headBgrd = ft_timelockstatistics(cfg, head_timelock{:}, bgrd_timelock{:});

%% Visualize?

figure(1), clf, hold on,
set(gcf, 'units','centimeters','position',[0 0 18 10] );
subplot(2,2,1), hold on, grid on,

fig = figure;
plot( stat_tfce_timelock.time, stat_tfce_timelock.stat_tfce, 'r');
%plot( stat_tfce_timelock.time(stat_tfce_timelock.mask), 400*ones(sum(stat_tfce_timelock.mask),1), 'r', 'linewidth',2 );
ylabel( 'TFCE' );
ylim( [-20, 420] );
xticks( 0:0.2:1 );

title( 'TFCE: H=2, E=0.5' );
set(gca, 'TickDir','out' );


fig = figure;
plot( stat_tfce_timelock2.time, stat_tfce_timelock2.stat_tfce, 'r');
%plot( stat_tfce_timelock.time(stat_tfce_timelock.mask), 400*ones(sum(stat_tfce_timelock.mask),1), 'r', 'linewidth',2 );
ylabel( 'TFCE' );
ylim( [-20, 420] );
xticks( 0:0.2:1 );

title( 'TFCE: H=2, E=0.5' );
set(gca, 'TickDir','out' );

fig = figure;
plot( stat_tfce_headBgrd.time, stat_tfce_headBgrd.stat_tfce, 'r');
%plot( stat_tfce_headBgrd.time(stat_tfce_headBgrd.mask), 400*ones(sum(stat_tfce_headBgrd.mask),1), 'r', 'linewidth',2 );
ylabel( 'TFCE' );
ylim( [-20, 420] );
xticks( 0:0.2:1 );

title( 'TFCE: H=2, E=0.5' );
set(gca, 'TickDir','out' );