%% Goal: Create ERSPs with different number of cycles - 23.08.2023
% Steps: 
% - Create the pathways
% - Load the subject data
% - TF Analysis

%% The code is written for the poster presentation in NGK conference, Ankara, Turkey.


%% ------------------------------------------------------------------------

%% 
% clear variables
% close all;
% clc;

%% 
% 
% savepath = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK';
% save_plots = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK/plots';
% 
% cd(savepath);
% 
% addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK/decoyPlots');
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
%close all;
clc;
%% Load necessary code
basepath='/net/store/nbp/projects/wd_ride_village/recordedData/wd_village/';
cd(basepath);

savepath = '/net/store/nbp/projects/wd_ride_village/processedData/village';
addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK');
addpath('/net/store/nbp/projects/wd_ride_village/processedData/village');
addpath('/net/store/nbp/projects/wd_ride_village/Analysis');

% load EEGlab
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
eeglab;
% load Fieldtrip
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');
ft_defaults

savedata_all = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK';

plotting = false;

% Directory to save the plots
plots_ersp = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK/plots';

if exist(plots_ersp, 'dir')
    addpath(genpath(plots_ersp));
else
    mkdir /net/store/nbp/projects/wd_ride_village/ERSP/paper/plots_ersp/;
    addpath(genpath(plots_ersp));
end

% %% Test the recursive dir
% 
% rootdir = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK/';
% 
% %filelist = dir(fullfile(rootdir, '**/*.*'));  %get list of files and folders in any subfolder
% len_channels_test = length(filelist);
% %filelist = filelist(~[filelist.isdir]);  %remove folders from list


%% get subjects to include
subjects = [1,2,5,12,15,16,17,19,20,21,22,27,29,32,33,34,37,38];
%% variables to save the data for all subjects
% variables to store the data for all subjects
data_all = [];

% each row will be filled with one subject data
%% for loop: load individual subjects
for s = 1:length(subjects)
    
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

% to get the name without the extension
setname=dir(fullfile(filepath,filename));
setname=setname.name(1:end-4); 

% load the csv file with the participant overview
rec_vill = readtable('recordings_village_corr.csv');
idx = -1;
% find the corresponding uid: firt get index 
for i = 1:height(rec_vill(:,2))
    if strcmp(rec_vill{i,2},[setname, '.xdf'])==1 
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



%%
EEG     = pop_loadset(cfg.ReData,cfg.mainfolder);
EEG     = eeg_checkset(EEG);

EEG_epoched = pop_epoch(EEG, {'0', '1', '2'}, [-0.875 1.175], 'epochinfo', 'yes');

%% create one file with all subject data while looping
% in the first loop the data_all is still empty 
% we go into the first statement and also save the data for subj 1
if isempty(data_all)
    data_all = EEG_epoched;    
else
    data_all = pop_mergeset(data_all, EEG_epoched);
end


end

%% Save the all subj data

pop_saveset(data_all, 'all_subjects', savedata_all); 
%% Transform to FT data

% filepath_alldata = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/NGK';
% filename_alldata = 'all_subjects.set';
% 
% cfg = [];
% cfg.dataset = [filepath_alldata,'/',filename_alldata];
% EEG_FT = ft_preprocessing(cfg);

EEG_FT = eeglab2fieldtrip(data_all, 'raw', 'coord_transform');

%% Wavelets!

cfg = [];
cfg.channel    = 'Oz';
cfg.method     = 'wavelet';
cfg.output     = 'pow';

%cfg.gwidth     = 3;
cfg.toi        = linspace(min(EEG_FT.time{1,1}),max(EEG_FT.time{1,1}), length(EEG_FT.time{1,1})); %matching data points in time depending on min and max in time in equal length
%cfg.foilim     = [1 45];
cfg.foi        = 1:0.5:45;
cfg.pad = 'nextpow2';

% 3 cycles

cfg.width = 3;
TFR_cyc3 = ft_freqanalysis(cfg, EEG_FT);

% 5 cycles
cfg.width = 5;
TFR_cyc5 = ft_freqanalysis(cfg, EEG_FT);

% 1 cycle
cfg.width = 1;
TFR_cyc1 = ft_freqanalysis(cfg, EEG_FT);

% % all channels ----
% cfg.channel    = 'all';
% cfg.width = 3;
% TFR_cyc3_all = ft_freqanalysis(cfg, EEG_FT);
% 
% % 5 cycles
% cfg.width = 5;
% TFR_cyc5_all = ft_freqanalysis(cfg, EEG_FT);
% 
% % 1 cycle
% cfg.width = 1;
% TFR_cyc1_all = ft_freqanalysis(cfg, EEG_FT);
% % 
%% Single plot

cfg          = [];
cfg.baseline     = [-0.5 -0.2]; %(keeping the average saccade duration out of the baseline)
cfg.baselinetype = 'db';
cfg.marker = 'on';
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.xlim     = [-0.5 0.8];
cfg.ylim     = [2 45];
cfg.zlim     =  [-0.8 0.8];
cfg.layout   = 'ordered' 
%cfg.channel  = 'Oz'; 
cfg.colormap = '*RdBu';

% Save plot - 1 - Oz
cfg.title = 'ERSPs at channel Oz - 1 cycle wavelet';
figure;
ft_singleplotTFR(cfg, TFR_cyc1);

xlabel('time (s)');
ylabel('frequency (Hz)');


file_name = 'TFR_cyc1_Oz.png';
file_path = fullfile(plots_ersp, file_name);
saveas(gcf, file_path);

% Save plot - 3 - Oz
cfg.title = 'ERSPs at channel Oz - 3 cycle wavelet';
figure;
ft_singleplotTFR(cfg, TFR_cyc3);

xlabel('time (s)');
ylabel('frequency (Hz)');

file_name = 'TFR_cyc3_Oz.png';
file_path = fullfile(plots_ersp, file_name);
saveas(gcf, file_path);

% Save plot - 5 - Oz
cfg.title = 'ERSPs at channel Oz - 5 cycle wavelet';
figure;
ft_singleplotTFR(cfg, TFR_cyc5);

xlabel('time (s)');
ylabel('frequency (Hz)');

file_name = 'TFR_cyc5_Oz.png';
file_path = fullfile(plots_ersp, file_name);
saveas(gcf, file_path);
% 
% 
% % Save plot - 1 - all
% cfg.channels =  'all';
% cfg.title = 'ERSP at all channels - 1 cycle wavelet';
% figure;
% ft_singleplotTFR(cfg, TFR_cyc1_all);
% 
% xlabel('time (s)');
% ylabel('frequency (Hz)');
% 
% 
% file_name = 'TFR_cyc1_allchannels.png';
% file_path = fullfile(save_plots, file_name);
% saveas(gcf, file_path);
% 
% % Save plot - 3 - all
% cfg.channels =  'all';
% cfg.title = 'ERSP at all channels - 3 cycle wavelet';
% figure;
% ft_singleplotTFR(cfg, TFR_cyc3_all);
% 
% xlabel('time (s)');
% ylabel('frequency (Hz)');
% 
% file_name = 'TFR_cyc3_allchannels.png';
% file_path = fullfile(save_plots, file_name);
% saveas(gcf, file_path);
% 
% % Save plot - 5 - all
% cfg.channels =  'all';
% cfg.title = 'ERSP at all channels - 5 cycle wavelet';
% figure;
% ft_singleplotTFR(cfg, TFR_cyc5_all);
% 
% xlabel('time (s)');
% ylabel('frequency (Hz)');
% 
% file_name = 'TFR_cyc5_allchannels.png';
% file_path = fullfile(save_plots, file_name);
% saveas(gcf, file_path);

%%
