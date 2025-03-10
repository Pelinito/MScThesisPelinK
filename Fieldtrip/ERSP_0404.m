% This code is written with Eva Marie von Butler

% fb = cwtfilterbank --> wavelet toolbox: costs something ....
% freqz(fb)
% --> trying it with fieldtrip toolbox 


% including fieltrip toolbox
restoredefaultpath
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');
ft_defaults


%% Load necessary code
cd /net/store/nbp/projects/wd_ride_village/ERSP/Pelin %folder to save matlab steps
savepath = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin';
savepath_old = '/net/store/nbp/projects/wd_ride_village/processedData/village';

addpath('/net/store/nbp/projects/wd_ride_village/processedData/village');
addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin');
addpath('/net/store/nbp/projects/wd_ride_village/Analysis');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
% load EEGlab (opens the blue eeglab window)
eeglab;

%% Choose which data to work on/ get participants number
% choose which steps you want to do
sub = input('Subject: ');
if sub<10
    subStr = '0'+string(sub);
else
    subStr = string(sub);
end


% add all necesarry tools, may not use all of them
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/preprocessing_helpers');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/NoiseTools');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0/plugins/xdfimport1.18/xdf');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources');


filtType   ='acausal';	% how do you want to filter? Options: acausal or causal
filtFreq   = 0.5;%1.0; %0.1;       % which frequency do you want to use for filtering? Recomended: 0.1Hz
eyetracking= 0;         % did you do eyetracking and want to integrate that data? This is not implemented in the EEG_Training but IntoTheWild

% Where is your data stored?
basepath='/net/store/nbp/projects/wd_ride_village/recordedData/wd_village/';

% Where is your trigger file?
%trgpath='/net/store/nbp/projects/wd_ride_village/processedData/village/TriggerInfo/';
trgpath='/net/store/nbp/projects/wd_ride_village/repos/wd-pilot-pipeline/data/village/processed/TriggerInfo_MAD';


%%% LOAD RAW DATA %%%
%select data file to work on

cd(basepath)
cntfile = dir(subStr+'_v_'+'*'+'.xdf');
cntpath = fullfile(basepath,cntfile.name);

if exist(cntpath,'file')
    [filepath,filename,ext] = fileparts(cntpath); % define variables
    % to keep the same as the manual selection procedure
    filepath = [filepath filesep];
    filename = [filename ext];
else
    cntpath = fullfile(basepath,sprintf('VP%u',sub),'*.*'); %MB probably relevant only when the EEG_training code chunk is there? Possibly can be deleted
    [filename, filepath]=uigetfile(cntpath);
end

%% delete? ------------------------
% % Construct file for the 2nd part if existing
% cntfile2 = dir(subStr+'_v2_'+'*'+'.xdf');
% cntpath2 = fullfile(basepath,cntfile2.name);
% 
% if exist(cntpath2,'file')
%     [filepath2,filename2,ext2] = fileparts(cntpath2); % define variables
%     % to keep the same as the manual selection procedure
%     filepath2 = [filepath2 filesep];
%     filename2 = [filename2 ext2]; 
% end
% -------------------------

% to get the name without the extension, used later for saving checkpoints
setname=dir(fullfile(filepath,filename));
setname=setname.name(1:end-4); 

% load the csv file with the participant overview
% load the actual name (the decrypted name, so no subject numbers are
% visible in the saved data?)
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
%mkdir(fullfile(savepath,'preprocessed/'))
str = [savepath, '/preprocessed/', uidname];
%mkdir(fullfile(str))
savedata = join(str,''); % create new path for each individual subject  (!!! warning !!! because this directorz alreadz exists from the other file as it creates the exact same file)

%% loading all preprocessing steps (include step 7) --> code from Debbies script: ---------------------------- I still need to change the filenames!
%% set cfg
cfg         = [];
cfg.avgref  = 1; % set it to 1 to rereference the data to avg ref
cfg.subj    = char(uidname); % get the name of the current subject 

cfg.mainfolder  = sprintf(savedata); % subject folder
% data = ['Raw,Rereference,Highpass0.1,Resampled'];
cfg.ChanAll    = sprintf('2_%s_bandpass_resample_deblank_ERSP.set',cfg.subj);
% data = ['Raw,Rereference,Highpass0.5,Resampled,ChanPos,ChannelReject'];
cfg.ChanFull    = sprintf('3_%s_channelrejTriggersXensor_ERSP.set',cfg.subj);
% data = ['Raw,Rereference,Highpass0.1,Resampled,ChanPos,ChannelReject,Cleaning,];
cfg.EEGfile     = sprintf('4_%s_Clean_ERSP.set', cfg.subj); % adjust for interchangable 
% mat file with components to reject
cfg.ICAfile     = sprintf('6_%s_ICAcleancont_ERSP.mat',cfg.subj);

%% Load data
if exist('EEG') == 1
    clear EEG
end
    
EEG     = pop_loadset(cfg.ChanFull,cfg.mainfolder);
EEG=eeg_checkset(EEG);
%% Reload trigger file
%trgpath='/net/store/nbp/projects/wd_ride_village/repos/wd-pilot-pipeline/data/village/processed/TriggerInfo_MAD';
EEG = pop_importevent(EEG,'event',fullfile(trgpath,strcat('TriggerFile_newTSdd_',uidname,'.csv')),'fields', {'latency', 'type'}, 'skipline', 1,'append','no'); % due to the 'append', 'no' the new trials overwrite the old trials and not just appended
EEG=eeg_checkset(EEG);

%% Exclude data segements
load(fullfile(savedata,sprintf('4_%s_cleaningTimes_ERSP.mat',uidname)),'tmprej','rej');
% reject the data
EEG = eeg_eegrej(EEG,tmprej(:,[3 4]));
EEG=eeg_checkset(EEG);

clear tmprej rej;

%% Exclude ICs
load(fullfile(cfg.mainfolder,cfg.ICAfile));
cfg.badcomponents = comps_to_rej;
mod = loadmodout12(fullfile(cfg.mainfolder,'amica')); % Load the output of AMICA from outdir

EEG.icaweights  = mod.W;
EEG.icasphere   = mod.S;
EEG.icawinv     = [];
EEG.icaact      = [];
EEG.icachansind = [];
EEG             = eeg_checkset(EEG);
% get the bad components out of the saved ICA file and re-reject them 
EEG             = pop_subcomp(EEG,comps_to_rej);

clear ic mod;

%% Adding Cz so it can be interpolated
% remove VEOG as it can not be interpolated
EEG = pop_select(EEG, 'nochannel', {'VEOG'}); 
tmp = pop_loadset('filename',cfg.ChanAll,'filepath',cfg.mainfolder,'loadmode','info');

% add Cz to the channel list 
if sum(strcmp({EEG.chanlocs.labels},'Cz')) == 0
    tmp.chanlocs(end+1)         = tmp.chanlocs(end);
    tmp.chanlocs(end).labels    = 'Cz';
    tmp.chanlocs(end).urchan    = 75;
    tmp.nbchan                  = tmp.nbchan + 1;
elseif sum(strcmp({tmp.chanlocs.labels},'Cz')) == 1
    error('Cz exists!')
end

tmp.chaninfo = [];
% apply the channelposition: BESA standard mapping     
tmp = pop_chanedit(tmp, 'lookup','/net/store/nbp/projects/EEG_Training/Analysis/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
% add this info now to cfg
cfg.urchanlocs= tmp.chanlocs;
%% Avg ref
if cfg.avgref
    EEG = pop_reref( EEG, []);
end
%% Interpolating missing channels
% all channels that are not needed
alldel = {'BIP1' 'BIP2' 'BIP3' 'BIP4' 'BIP5' 'BIP6' 'BIP7' 'BIP8' 'AUX1' 'AUX2' 'AUX3' 'AUX4' 'AUX5' 'AUX6' 'AUX7' 'AUX69' 'AUX70' 'AUX71' 'AUX72' 'IZ' 'TIME' 'L-GAZE-X' 'L-GAZE-Y' 'L-AREA' 'R-GAZE-X' 'R-GAZE-Y' 'R-AREA' 'L_GAZE_X' 'L_GAZE_Y' 'L_AREA', 'R_GAZE_X' 'R_GAZE_Y' 'R_AREA' 'INPUT' 'BIP65' 'BIP66' 'BIP67' 'BIP68' 'Reference'};

% get rid of all the channels that are not needed 
idxs = ~ismember({cfg.urchanlocs.labels},alldel);
%EEG= pop_interp(EEG,cfg.urchanlocs(idxs),'spherical');
EEG = pop_interp(EEG, eeg_mergelocs(cfg.urchanlocs(idxs)), 'spherical');

EEG = eeg_checkset(EEG);
 % check if duplicate channel label
if isfield(EEG.chanlocs, 'labels')
    if length( { EEG.chanlocs.labels } ) > length( unique({ EEG.chanlocs.labels } ) )
        disp('Warning: some channels have the same label');
    end
end

%% Filter EEG Data
EEG = pop_eegfiltnew(EEG, [], 40);
%save(fullfile(savedata,sprintf('EEG_forPython.mat',uidname)),'EEG');

%% Change the type in the events ------------------------------- last part of Debbies code
% 2 = 'bgrd'; 1 = 'body'; 0 = 'face'
for evt = 1:length(EEG.event)
    if EEG.event(evt).type     == '0'
        EEG.event(evt).type     = 'head';
    elseif EEG.event(evt).type == '1'
        EEG.event(evt).type     ='body';
    elseif EEG.event(evt).type == '2'
        EEG.event(evt).type     ='bgrd';
    end
end
disp('changing type done')

% ---------------------------

%% loading EEG data set
%EEG = pop_loadset(sprintf('4_%s_Clean_ERSP.set',uidname),fullfile(savedata));

% load ICA results ---------- ERRORS ---------------------------
%icapath = fullfile(savedata,'amica');
%mod = loadmodout12(icapath); 

%%% APPLY ICA WEIGHTS TO DATA %%%
%EEG.icasphere = mod.S;
%EEG.icaweights = mod.W;


%% EPOCH DATA %%% epoch only head trials

event_codes = {'head', 'body', 'bgrd'};
%events = EEG.event(strcmp({EEG.event.type}, 'head'));
EEG_epoched = pop_epoch(EEG, event_codes, [-0.875 1.175], 'epochinfo', 'yes');

[EEG_head, EEG_body, EEG_bground] = separateevents(EEG_epoched);

% Select the channels of interest for each condition -- in this case PO7
% and PO8
EEG_head_PO7 = pop_select(EEG_head, 'channel', {'PO7'});
EEG_head_PO8 = pop_select(EEG_head, 'channel', {'PO8'});

EEG_body_PO7 = pop_select(EEG_body, 'channel', {'PO7'});
EEG_body_PO8 = pop_select(EEG_body, 'channel', {'PO8'});

EEG_bground_PO7 = pop_select(EEG_bground, 'channel', {'PO7'});
EEG_bground_PO8 = pop_select(EEG_bground, 'channel', {'PO8'});


% eeglab;
% EEG_epoch_try = pop_selectevent(EEG, 'event', 1);
% EEG_epoched2 = pop_epoch(EEG, {'body', 'boundary'}, [-0.875 1.175], 'epochinfo', 'yes');

%save files
%EEG = pop_editset(EEG, 'setname', sprintf('5_%s_ICAEpoched_ERSP',uidname));
%EEG = pop_saveset(EEG, 'filename',sprintf('5_%s_ICAEpoched_ERSP',uidname),'filepath',fullfile(savedata));

%% load EEGLAB data into fieldtrip
% converting EEGLAB data into fieldtrip
%data_fieldtrip = eeglab2fieldtrip(EEG, 'raw', 'coord_transform');

% data used in eeglab2fieldtrip is the merged sets from
% all_subj_computation

data_ft_head_PO7 = eeglab2fieldtrip(EEG_merged_head_PO7, 'raw', 'coord_transform');
data_ft_head_PO8 = eeglab2fieldtrip(EEG_merged_head_PO8, 'raw', 'coord_transform');

data_ft_body_PO7 = eeglab2fieldtrip(EEG_merged_body_PO7, 'raw', 'coord_transform');
data_ft_body_PO8 = eeglab2fieldtrip(EEG_merged_body_PO8, 'raw', 'coord_transform');

data_ft_bground_PO7 = eeglab2fieldtrip(EEG_merged_bground_PO7, 'raw', 'coord_transform');
data_ft_bground_PO8 = eeglab2fieldtrip(EEG_merged_bground_PO8, 'raw', 'coord_transform');

%saving the fieldtrip data so we do not have to load it everytime
%save('/net/store/nbp/projects/wd_ride_village/ERSP/Eva/fieldtrip_data', 'data_fieldtrip')
%load('/net/store/nbp/projects/wd_ride_village/ERSP/Eva/fieldtrip_data.mat')

%cfg = [];
%cfg.dataset = data_fieldtrip  --> not necessary
%shows us the dataset in fieldtrip (data scroll)
%data = ft_databrowser(cfg,data_fieldtrip)

%% ------------------------------------
% % Hanning tapering - just to visualise the power spectrum
% cfg = [];
% cfg.output = 'pow';
% cfg.channel = 'FPz';
% cfg.method = 'mtmconvol';
% cfg.foi   = 1:30;
% %cfg.t_ftmwin = ones(length(cfg.foi),1).*0.5;
% cfg.toi = -1:0.05:1;
% cfg.taper = 'hanning';
% 
% hanning_ = ft_freqanalysis(cfg,data_fieldtrip);
% 
% cfg. = [];
% cfg.baselinetzpe = 'absolute';
% cfg.maskstyle = 'saturation';
% %figure;
% %hold on;
% %plot(hanning_.freq(2:19),hanning_.powspctrm(2:19));
% %ft_singleplotER(cfg, hanning_);
% ft_singleplotTFR(cfg, hanning_);
% ylim([0 20])
% xlim([0 20])
% xlabel('Frequency (Hz)')
% ylabel('Power (µV)')

%% ------------------------

%% Using morlet wavelets for time-frequency analyis
cfg = [];
cfg.channel    = 'all';
cfg.method     = 'wavelet';
cfg.output     = 'pow';
cfg.width      = 3;
cfg.gwidth     = 3;
cfg.foilim     = [1 45]; %instead of cfg.foi
%cfg.foi        = 1:2:45;
cfg.pad = 'nextpow2';
%cfg.keeptrials = 'no';
%cfg.polyremoval = 1;
%cfg.keeptapers = 'no';

cfg.toi = linspace(min(data_fieldtrip_head.time{1,1}),max(data_fieldtrip_head.time{1,1}), length(data_fieldtrip_head.time{1,1})); %matching data points in time depending on min and max in time in equal length
wavelet_head = ft_freqanalysis(cfg, data_fieldtrip_head);

cfg.toi = linspace(min(data_fieldtrip_body.time{1,1}),max(data_fieldtrip_body.time{1,1}), length(data_fieldtrip_body.time{1,1})); %matching data points in time depending on min and max in time in equal length
wavelet_body = ft_freqanalysis(cfg, data_fieldtrip_body);

cfg.toi = linspace(min(data_fieldtrip_bground.time{1,1}),max(data_fieldtrip_bground.time{1,1}), length(data_fieldtrip_bground.time{1,1})); %matching data points in time depending on min and max in time in equal length
wavelet_bground = ft_freqanalysis(cfg, data_fieldtrip_bground);

%----------------------------------------not necessary (for-loop try for different widths)
%freq = [];
%cfg.parameter = 'powspctrm';
% 
% for f = cfg.foi
%     cfg.width = max(3, round(3*f/10));
%     tmpfreq = ft_freqanalysis(cfg, data_fieldtrip);
%     freq = ft_appenddata(freq, tmpfreq);
% end
%---------------------------------------

%% ERP

cfg = [];
cfg.latency = [min(data_fieldtrip_head.time{1,1}) max(data_fieldtrip_head.time{1,1})];
cfg.channel = 'all';

% Latency is the same for all - no need to change again for each one

ERP_head_PO7 = ft_timelockanalysis(cfg, data_ft_head_PO7);
ERP_head_PO8 = ft_timelockanalysis(cfg, data_ft_head_PO8);


ERP_body_PO7 = ft_timelockanalysis(cfg, data_ft_body_PO7);
ERP_body_PO8 = ft_timelockanalysis(cfg, data_ft_body_PO8);

ERP_bground_PO7 = ft_timelockanalysis(cfg, data_ft_bground_PO7);
ERP_bground_PO8 = ft_timelockanalysis(cfg, data_ft_bground_PO8);

cfg =[];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'absolute';
cfg.title = 'ERP at all channels'; % Change if the channel selection is 
                                   % changed, might want to mention
                                   % the reference channel - No idea if
                                   % needed

figure;
hold on;
ft_singleplotER(cfg, ERP_head_PO7);
ft_singleplotER(cfg, ERP_body_PO7);
ft_singleplotER(cfg, ERP_bground_PO7);
xlabel('Time (ms)');
ylabel('Potential (μV)');


% Plotting PO7 activity for all conditions
x = ERP_head_PO7.time;
y1 = mean(ERP_head_PO7.var, 1);
y2 = mean(ERP_body_PO7.var, 1);
y3 = mean(ERP_bground_PO7.var, 1);
figure;
plot(x,y1,'g',x,y2,'r',x,y3,'b');
xlim([-0.5 0.8])
xline(0)
legend('Head', 'Body', 'Background', 'Gaze onset'); % face or head?
xlabel('Time (ms)');
ylabel('Potential (μV)');
title('ERPs at PO7 from 2 subjects for all conditions');


% Plotting PO8 activity for all conditions
x = ERP_head_PO8.time;
y1 = mean(ERP_head_PO8.var, 1);
y2 = mean(ERP_body_PO8.var, 1);
y3 = mean(ERP_bground_PO8.var, 1);
figure;
plot(x,y1,'g',x,y2,'r',x,y3,'b');
xlim([-0.5 0.8])
xline(0)
legend('Head', 'Body', 'Background', 'Gaze onset'); % face or head?
xlabel('Time (ms)');
ylabel('Potential (μV)');
title('ERPs at PO8 from 2 subjects for all conditions');

% Plotting HEAD at PO7 and PO8 to compare
x = ERP_head_PO8.time;
y1 = ERP_head_PO7.var;
y2 = ERP_head_PO8.var;
figure;
plot(x,y1,'r',x,y2,'b');
xlim([-0.5 0.8])
xline(0)
legend('PO7', 'PO8', 'Gaze onset'); % face or head?
xlabel('Time (ms)');
ylabel('Potential (μV)');
title('ERPs at PO7 and PO8 from 2 subjects for the HEAD condition');

% Plotting BODY at PO7 and PO8 to compare
x = ERP_body_PO8.time;
y1 = ERP_body_PO7.var;
y2 = ERP_body_PO8.var;
figure;
plot(x,y1,'r',x,y2,'b');
xlim([-0.5 0.8])
xline(0)
legend('PO7', 'PO8', 'Gaze onset'); % face or head?
xlabel('Time (ms)');
ylabel('Potential (μV)');
title('ERPs at PO7 and PO8 from 2 subjects for the BODY condition');

x = ERP_bground_PO8.time;
y1 = ERP_bground_PO7.var;
y2 = ERP_bground_PO8.var;
figure;
plot(x,y1,'r',x,y2,'b');
xlim([-0.5 0.8])
xline(0)
legend('PO7', 'PO8', 'Gaze onset'); % face or head?
xlabel('Time (ms)');
ylabel('Potential (μV)');
title('ERPs at PO7 and PO8 from 2 subjects for the BACKGROUND condition');


%xline(0);
pop_timtopo(EEG_merged_head_PO7, [-500 800]); % save as 'timtopo'





 %% Using Superlets (morlet wavelets) for time-frequency analyis
% cfg = [];
% cfg.channel    = 'all';
% cfg.method     = 'superlet';
% cfg.output     = 'pow';
% %cfg.foilim     = [3 40]; %instead of cfg.foi
% cfg.foi        = 1:30;
% cfg.width      = 7;
% cfg.gwidth     = 3;
% cfg.superlet.smooth_factor = 10;
% cfg.toi        = -0.5:0.01:1.5 %linspace(min(data_fieldtrip.time{1,1}),max(data_fieldtrip.time{1,1}), length(data_fieldtrip.time{1,1})); %matching data points in time depending on min and max in time in equal length
% cfg.combine    = 'multiplicative';
% %cfg.order      = ft_getopt(cfg, 'order', ones(1,numel(cfg.foi))); %linspace(1,5,numel(cfg.foi)); --> it runs but the plots are not working anymore
% TFRwave_sup = ft_freqanalysis(cfg, data_fieldtrip);
% 
% %plotting the superlet analysis:
% cfg = [];
% cfg.baseline = [-0.5 -0.2];
% cfg.baslinetype = 'absolute';
% cfg.channel  = 'Oz'; 
% ft_singleplotTFR(cfg, TFRwave_sup);


% %% Plotting topo plot left 
% cfg          = [];
% %cfg.baseline     = [-0.5 -0.2]; --> baseline gives error
% %cfg.baselinetype = 'absolute';
% cfg.colorbar = 'yes';
% cfg.parameter = 'powspctrm';
% cfg.xlim     = [-0.5 : 0.5];
% cfg.ylim     = [8 12];
% cfg.channel  = {'all','-FP2','-AF4','-AF8','-F2','-F4','-F6','-F8','-FC2','-FC4','-FC6','-FT8','-C2','-C4','-C6','-T8','-CP2','-CP4','-CP6','-TP8','-M2','-P2','-P4','-P6','-P8','-PO4','-PO6','-PO8','O2','FPz','Fz','FCz','Cz','CPz','Pz','POz','Oz'}; %selecting all apart from the electrodes on the right
% cfg.title    = 'Test topo left'
% 
% ft_topoplotTFR(cfg, TFRwave);
% 
% 
% 
% %% Plotting topo plot right 
% cfg          = [];
% %cfg.baseline     = [-0.5 -0.2];
% %cfg.baselinetype = 'absolute';
% cfg.colorbar = 'yes';
% cfg.parameter = 'powspctrm';
% cfg.xlim     = [-0.5 : 0.5];
% cfg.ylim     = [8 12];
% cfg.channel  = {'all','-FP1','-AF3','-AF7','-F1','-F3','-F5','-F7','-FC1','-FC3','-FC5','-FT7','-C1','-C3','-C5','-T7','-CP1','-CP3','-CP5','-TP7','-M1','-P1','-P3','-P5','-P7','-PO3','-PO5','-PO7','O1','FPz','Fz','FCz','Cz','CPz','Pz','POz','Oz'};
% xl
% figure;
% ft_topoplotTFR(cfg, TFRwave);untitled
% title('Test topo right');
% 
% 

%% Plotting a singleplot 
cfg          = [];
cfg.baseline     = [-0.5 -0.2]; %[-0.4 -0.1]
cfg.baselinetype = 'db';
cfg.marker = 'on';
%cfg.showlabels   = 'yes';
%cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.colormap = 'jet';
cfg.zlim     = [-1.5 1.5];
cfg.layout   = 'ordered' %'natmeg_customized_eeg1005.lay';
cfg.channel  = 'all'; %what channel should we use
%cfg.interactive = 'no';

figure;
ft_singleplotTFR(cfg, wavelet_head);
title('ERSP power spectrum of HEAD CONDITION at all channels');
xlim([-0.5 0.8])
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

ft_singleplotTFR(cfg, wavelet_body);
title('ERSP power spectrum of BODY CONDITION at all channels');
xlim([-0.5 0.8])
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

ft_singleplotTFR(cfg, wavelet_bground);
title('ERSP power spectrum of BACKGROUND CONDITION at all channels');
xlim([-0.5 0.8])
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

% %% Plotting a multiplot for all channels for frequency 3-40
% cfg          = [];
% cfg.baseline     = [-0.5 -0.2];
% cfg.baselinetype = 'absolute';
% cfg.showlabels   = 'yes';
% cfg.parameter = 'powspctrm';
% cfg.colorbar = 'yes';
% cfg.zlim     = 'maxabs';
% 
% figure;
% ft_multiplotTFR(cfg, TFRwave);
% title('Test ERSP');
% xlabel('Time (ms)');
% ylabel('Frequency (Hz)');
% 
