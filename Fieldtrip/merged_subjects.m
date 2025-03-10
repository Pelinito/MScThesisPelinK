%% SCRIPT:
% Loop all subjects to get their latest preprocessed file - final one
% and separate the conditions to append into a merged data of all subjects
% for that condition.

% Epoching and channel selection will follow the merged data.

% Script mostly written by me to be used with Eva Marie von Butler

%% Beginning

clear variables
close all;
clc;

%% Add data paths

basepath='/net/store/nbp/projects/wd_ride_village/ERSP/Pelin';
cd(basepath);

savepath = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin';
addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin');
addpath('/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/preprocessedfromDebbie');
%addpath('/net/store/nbp/projects/wd_ride_village/processedData/village');
%addpath('/net/store/nbp/projects/wd_ride_village/Analysis');
% load EEGlab
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
eeglab;
% load Fieldtrip
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');
ft_defaults

savedata_all = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/preprocessed/';


plotting = false;

%% Initiate the empty arrays to merge at the end of the loop

% Change the number of channels to 70, the merging sets need to have the
% same number of channels and sampling rates

% % EEG_merged_head = eeg_emptyset();
% % EEG_merged_head.nbchan = 70;
% % EEG_merged_head.srate = 512;
% % 
% % EEG_merged_body = eeg_emptyset();
% % EEG_merged_body.nbchan = 70;
% % EEG_merged_body.srate = 512;
% % 
% % EEG_merged_bground = eeg_emptyset();
% % EEG_merged_bground.nbchan = 70;
% % EEG_merged_bground.srate = 512;


%% get subjects to include --> change!!!
subjects = [1,2,5,12,15,16,19,20,21,22,27,29,32,33,34,37,38];%,42,43,44,46];

% Initiate empty tables for the Summary_EEG and fill them

% All channel info
summary_eeg_head_all = zeros(length(subjects), 1050); % 1050 bcs of the time points we have
summary_eeg_body_all = zeros(length(subjects), 1050);
summary_eeg_bground_all = zeros(length(subjects), 1050);


    
% Only for POz
summary_eeg_head_poz = zeros(length(subjects), 1050);
summary_eeg_body_poz = zeros(length(subjects), 1050);
summary_eeg_bground_poz = zeros(length(subjects), 1050);



%[1,2,5,12,15,16,17,19,20,21,22,27,29,32,33,34,37,38]
for s = 1:length(subjects)  % Subj name turns into string
    
    %sub = subjects(s);
    if subjects(s)<10
        subStr = '0'+string(subjects(s));
    else
        subStr = string(subjects(s)); 
    end
    
    % Read the participant data
    
    rec_vill = readtable('recordings_village_corr.csv');
    
    % Find the corresponding uidname from the subStr
%     setname = dir(subStr+"*"+".xdf");
%     setname = setname.name;
    %setname = subStr+'*'+'.xdf';
    
    
    idx = -1;
    % find the corresponding uid: firt get index 
    
    for i = 1:height(rec_vill(:,2))
        
        cell = rec_vill.file(i);
        cell = char(cell);
        
        
        if strcmp(cell(1:2), subStr)==1 
            idx = i;
        end
    end
    
    
    if  idx >= 0
        uidname = rec_vill{idx,1};
    else
    %throw error if the uid could not be found
        error('No uid found')
    end

    %convert uidname and cut the ending 
    uidname = uidname{1,1};
    
    % Reach the folder with the uidname and get the 7th step set
    
       
    str = [savepath, '/preprocessedfromDebbie/', uidname];
    if exist(str,'dir') == 0
        continue
    end
    cd(str); % We are in the folder for this uidname
    
    selected_file = dir('7_'+"*"+".set"); % Get the 7th step set
    % file_path = fullfile(selected_file.name);

    % Load the data 
    EEG_new = pop_loadset(selected_file.name, selected_file.folder); 


    
      % Change the type in the events
    % 2 = 'bgrd'; 1 = 'body'; 0 = 'face'
    for evt = 1:length(EEG_new.event)
        if EEG_new.event(evt).type     == '0'
            EEG_new.event(evt).type     = 'head';
        elseif EEG_new.event(evt).type == '1'
            EEG_new.event(evt).type     ='body';
        elseif EEG_new.event(evt).type == '2'
            EEG_new.event(evt).type     ='bgrd';
        end
    end
    
    disp('Changing type done')
    
    
    % Epoch the data
    
    event_codes = {'head', 'body', 'bgrd'};
    
    EEG_epoched = pop_epoch(EEG_new, event_codes, [-0.875 1.175], 'epochinfo', 'yes');
        
    % Separate the conditions to different datasets

    [EEG_head, EEG_body, EEG_bground] = separateevents(EEG_epoched);
    
    % Take ALL CHANNELS mean for the particular subject
    
    % Take the mean across trials 
    mean_sum_head_all = mean(EEG_head.data, 3); % results in 70x1050 table
    
    % Take the mean across channels
    mean_sum_head_all = mean(mean_sum_head_all, 1); % results in 1x1050 vector
    
    
    % Repeat it for the other conditions
    mean_sum_body_all = mean(EEG_body.data, 3);
    mean_sum_body_all = mean(mean_sum_body_all, 1);
    
    mean_sum_bground_all = mean(EEG_bground.data, 3);
    mean_sum_bground_all = mean(mean_sum_bground_all, 1);    
    
    % Take POz CHANNEL from the particular subject - Channel location is 28
    % For the specific channels, filter the channel you want then take
    % the mean across trials
    
    mean_sum_head_poz = EEG_head.data(28, :, :); % Extract the channel
    mean_sum_head_poz = mean(mean_sum_head_poz, 3); % Mean across trials
    
    
    mean_sum_body_poz = EEG_body.data(28, :, :);
    mean_sum_body_poz = mean(mean_sum_body_poz, 3);
    
    
    mean_sum_bground_poz = EEG_bground.data(28, :, :);
    mean_sum_bground_poz = mean(mean_sum_bground_poz, 3);    
            
    % Put the summarized values to the tables
    summary_eeg_head_all(s,:) = mean_sum_head_all;
    summary_eeg_body_all(s,:) = mean_sum_body_all;
    summary_eeg_bground_all(s,:) = mean_sum_bground_all;
    
    summary_eeg_head_poz(s,:) = mean_sum_head_poz;
    summary_eeg_body_poz(s,:) = mean_sum_body_poz;
    summary_eeg_bground_poz(s,:) = mean_sum_bground_poz;
        
        
        
        
    
    
    
    % ---------------------------------------------------------------------
    % SKIP MERGING
   % %Merge here
%     
%     if exist('EEG_merged','var') == 0
%         EEG_merged = EEG_new;
%         
%     else
%         EEG_merged = pop_mergeset(EEG_merged, EEG_new, 1);
%     end
    
end                                % End of the subject loop
   %% ---------------------------------------------------------------------
    
%     % Change the type in the events
%     % 2 = 'bgrd'; 1 = 'body'; 0 = 'face'
%     for evt = 1:length(EEG_merged.event)
%         if EEG_merged.event(evt).type     == '0'
%             EEG_merged.event(evt).type     = 'head';
%         elseif EEG_merged.event(evt).type == '1'
%             EEG_merged.event(evt).type     ='body';
%         elseif EEG_merged.event(evt).type == '2'
%             EEG_merged.event(evt).type     ='bgrd';
%         end
%     end
%     
%     disp('Changing type done')
%     
%     
%     % Epoch the data
%     
%     event_codes = {'head', 'body', 'bgrd'};
%     
%     EEG_epoched = pop_epoch(EEG_merged, event_codes, [-0.875 1.175], 'epochinfo', 'yes');
%         
%     % Separate the conditions to different datasets
% 
%     [EEG_head, EEG_body, EEG_bground] = separateevents(EEG_epoched);
    
    %% --------------------------------------------------------------------
    
    %pop_erpimage(EEG_new_head,1);
    
%     data_FT_head = eeglab2fieldtrip(EEG_new_head, "raw", "coordtransform");
%     data_FT_body = eeglab2fieldtrip(EEG_new_body, "raw", "coordtransform");
%     data_FT_bground = eeglab2fieldtrip(EEG_new_bground, "raw", "coordtransform");
    
    % ------------- REST IS FIELDTRIP, DON'T USE IT NOW (12.04)------------
%     cfg = [];
%     cfg.latency = [min(data_FT_head.time{1,1}) max(data_FT_head.time{1,1})];
%     cfg.channel = 'all';
% 
%     ERP_head = ft_timelockanalysis(cfg, data_FT_head);
%     ERP_body = ft_timelockanalysis(cfg, data_FT_body);
%     ERP_bground= ft_timelockanalysis(cfg, data_FT_bground);

% Statistical analysis ??


    % Plot
    
%      ymin = -2 * 10e-8;
%      ymax = 2 * 10e-8;
%     cfg = [];
%     %cfg.parameter = 'powspctrm';
%     cfg.channel = 'PO8';
%     cfg.baseline = 'yes';
%     cfg.baselinetype = 'absolute';
%     cfg.xlim = [-0.5 0.8];
%     %cfg.ylim = [ymin ymax];

% 
%     cfg.title = 'ERPs - all channels - HEAD STIMULI - ' + subStr;
%     ft_singleplotER(cfg, ERP_head);
%     xlabel('Time (ms)');
%     ylabel('Potential (μV)');

%     cfg.title = 'ERPs - all channels - BODY STIMULI - '+ subStr;
%     ft_singleplotER(cfg, ERP_body);
%     xlabel('Time (ms)');
%     ylabel('Potential (μV)');
% 
%     cfg.title = 'ERPs - all channels - BACKGROUND STIMULI - '+ subStr;
%     ft_singleplotER(cfg, ERP_bground);
%     xlabel('Time (ms)');
%     ylabel('Potential (μV)');
    
    
    
%     % Merge the merged set with the new set - 
%     
%     if exist('EEG_merged_head','var') == 1
%         
%         EEG_merged_head = pop_mergeset(EEG_merged_head, EEG_new_head, 1);
%     else
%         EEG_merged_head = EEG_new_head;
%     end
% % %     
% % %     if exist ('EEG_merged_body', 'var') == 1
% % %         
% % %         EEG_merged_body = pop_mergeset(EEG_merged_body, EEG_new_body, 1);
% % %     else
% % %         EEG_merged_body = EEG_new_body;
% % %     end
% % %     
% % %     if exist('EEG_merged_bground', 'var') == 1
% % %         
% % %         EEG_merged_bground = pop_mergeset(EEG_merged_bground, EEG_new_bground, 1);
% % %     else
% % %         EEG_merged_bground = EEG_new_bground;
% % %     end



%disp('MERGING COMPLETED - PROCEED WITH FIELDTRIP');

%% Data transformation to Fieldtrip

data_FT_head = eeglab2fieldtrip(EEG_head, "raw", "coordtransform");
data_FT_body = eeglab2fieldtrip(EEG_body, "raw", "coordtransform");
data_FT_bground = eeglab2fieldtrip(EEG_bground, "raw", "coordtransform");

%% ERP in Fieldtrip

% Time analysis
cfg = [];
cfg.latency = [min(data_FT_head.time{1,1}) max(data_FT_head.time{1,1})];
cfg.channel = 'all';

ERP_head = ft_timelockanalysis(cfg, data_FT_head);
ERP_body = ft_timelockanalysis(cfg, data_FT_body);
ERP_bground= ft_timelockanalysis(cfg, data_FT_bground);

% Statistical analysis ??


% Plot
figure;
hold on

cfg = [];
%cfg.parameter = 'average';
cfg.channel = 'PO8';
cfg.baseline = 'yes';
cfg.xlim = [-0.2 0.5];
cfg.baselinetype = 'absolute';


cfg.title = 'ERPs - PO8 - HEAD STIMULI' ;
ft_singleplotER(cfg, ERP_head);
xlabel('Time (ms)');
ylabel('Potential (μV)');

cfg.title = 'ERPs - PO8 - BODY STIMULI' ;
ft_singleplotER(cfg, ERP_body);
xlabel('Time (ms)');
ylabel('Potential (μV)');

cfg.title = 'ERPs - PO8 - BACKGROUND STIMULI';
ft_singleplotER(cfg, ERP_bground);
xlabel('Time (ms)');
ylabel('Potential (μV)');


%% ERSP in Fieldtrip

cfg = [];
cfg.channel    = 'all';
cfg.method     = 'wavelet';
cfg.output     = 'pow';
cfg.width      = 3;
cfg.gwidth     = 3;
cfg.foilim     = [1 45]; %instead of cfg.foi
%cfg.foi        = 1:2:45;
cfg.pad = 'nextpow2';

cfg.toi = linspace(min(data_FT_head.time{1,1}),max(data_FT_head.time{1,1}), length(data_FT_head.time{1,1})); %matching data points in time depending on min and max in time in equal length
wavelet_head = ft_freqanalysis(cfg, data_FT_head);
wavelet_body = ft_freqanalysis(cfg, data_FT_body);
wavelet_bground = ft_freqanalysis(cfg, data_FT_bground);

cfg          = [];
cfg.baseline     = [-0.5 -0.2]; %[-0.4 -0.1]
cfg.baselinetype = 'db';
cfg.marker = 'on';
%cfg.showlabels   = 'yes';
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.colormap = 'jet';
%cfg.zlim     = [-1.5 1.5];
%cfg.layout   = 'ordered' %'natmeg_customized_eeg1005.lay';
cfg.channel  = 'PO8'; %what channel should we use
%cfg.interactive = 'no';

figure;
ft_singleplotTFR(cfg, wavelet_head);
title('ERSP power spectrum of HEAD CONDITION at PO8');
xlim([-0.5 0.8])
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

figure;
ft_singleplotTFR(cfg, wavelet_body);
title('ERSP power spectrum of BODY CONDITION at PO8');
xlim([-0.5 0.8])
xlabel('Time (ms)');
ylabel('Frequency (Hz)');


%% EEGLAB Plottings 
pop_erpimage(EEG_head, 1);



pop_timtopo(EEG_merged_head, [-200 500]); % save as 'timtopo'

pop_topoplot(EEG_merged_head, 1, [0, 76, 100, 200]); % save as

plot(EEG_merged_head.times,ort_all)

%%

%% SUMMARY EEG INFO



