% path itself
%cd C:\Users\pelin\Desktop\Tez;
%datapath = "C:\Users\pelin\Desktop\Tez";

cd /net/store/nbp/projects/wd_ride_village/ERSP/Pelin
datapath = "C:\Users\pelin\Desktop\Tez";



addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/eeglab2020_0/plugins/xdfimport1.18/xdf');
addpath('/net/store/nbp/projects/wd_ride_village/Matlab-resources/fieldtrip-20230215');


% path for the images - tidy up
plots_ft = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/plots_ft';

% Note: save_plots_here is unnecessary but different namings help me

if exist(plots_ft, 'dir')
    save_plots_here = plots_ft;
    addpath(genpath(save_plots_here));
else
    mkdir /net/store/nbp/projects/wd_ride_village/ERSP/Pelin/plots_ft;
    save_plots_here = plots_ft;
    addpath(genpath(save_plots_here));

end

%% EEGLAB - Load Data

% subjID = input("Enter Subject ID: "); % Only option for now is 2
% 
% step = input("Enter which data you want to see: ");
% % Steps are: 
% % 0 - Raw datadatapath = 
% % 2 - Bandpass, resample
% % 3 - Channel reject
% % 4 - Cleaned
% 
% subjID = string(subjID);
% step = string(step);
% %preprocessed_folder = "C:\Users\pelin\Desktop\Tez\preprocessed";
% 
% uidname = "391543f6-7ff0-4ed3-aefd-c0ef68ebbbf3";

% TO DO: Fix the backslashes to the strings above
% selected_file = dir(preprocessed_folder+"\"+subjID+"*"+"\"+step+"*"+".set");
% file_path = fullfile(selected_file.name);

% Load the data - For the only clean data: Step selection is 4
% EEG = pop_loadset(selected_file.name, selected_file.folder); 
EEG = pop_loadset(391543f6-7ff0-4ed3-aefd-c0ef68ebbbf3_Clean_ERSP.set);

%eegplot(EEG.data);

%% Epoch the data
EEG_epoched = pop_epoch(EEG, {0,1,2}, [-0.2 0.5], 'epochinfo', 'yes');

%% EEGLAB to Fieldtrip

data_FT = eeglab2fieldtrip(EEG_epoched, "raw", "coordtransform");

%% Wavelets - Alpha band for now - changed to all freq for testing (03.03)

cfg = [];
cfg.channel = 'all';
cfg.method = 'superlet';
cfg.output = 'pow';
%cfg.foilim = [1 30]; % Frequency limits for the alpha band
cfg.foi = 1:0.1:40;
cfg.toi = linspace(min(data_FT.time{1,1}), max(data_FT.time{1,1}), length(data_FT.time{1,1}));
cfg.width = 3;
cfg.gwidth = 3;
cfg.combine = 'multiplicative';
cfg.order = linspace(1,5,numel(cfg.foi)); %default is ones(1, numel(cfg.foi))
cfg.pad = 'nextpow2';

%wavelet_ = ft_freqanalysis(cfg, data_FT);


% variables for specifications - for plot names

% if isfield(cfg, "foilim")
%     frq_str = sprintf("_%0.f_%0.f", cfg.foilim(1,1), cfg.foilim(1,end));
% elseif isfield(cfg, "foi")
%     frq_str = sprintf("_%0.f_%0.f", cfg.foi(1,1), cfg.foi(1,end));
% end
% 
% frq_str = frq_str + "_";


wavelet_ = ft_freqanalysis(cfg, data_FT);

%% Plot - Multiplot

cfg = [];
cfg.baseline     = [-0.1 0];
cfg.baselinetype = 'absolute';
%cfg.marker       = 'on';
%cfg.layout       = 'easycapM10.mat'; - Maybe we have a match to one of the options?
cfg.parameter = 'powspctrm';
cfg.interactive  = 'no';
cfg.colorbar = 'yes';
%cfg.zlim = [0 20];
cfg.trials = "all";
cfg.masknans = 'yes';
cfg.showlabels = 'yes';
cfg.title = "Frequency: 1:0.1:40. Order: 5";
%cfg.title = "Subj_" + subjID + "_" + "multiplot" + frq_str + baseline_str + trials_str + ".png";

ft_multiplotTFR(cfg, wavelet_);

% Save the plot
% Needs: Plot name - spritntf with the specifications, file directory
baseline_str = sprintf("_%0.1f_%0.1f_", cfg.baseline(1,1), cfg.baseline(1,2));
trials_str = "trials_" + cfg.trials;

% plot name

plot_name = "multiplicative_09.03_superlet_Subj_" + subjID + "_" + "multiplot" + frq_str + baseline_str + trials_str + ".png";

plot_fullname = fullfile(save_plots_here, plot_name);

plot = gcf;
exportgraphics(plot, plot_fullname);

%% Plot - Singleplot

cfg = [];
cfg.paramete = 'powspctrm';
cfg.baseline     = [-0.1 0];
cfg.baselinetype = 'absolute';
cfg.trials = 'all';
cfg.masknans = 'yes';
cfg.colorbar = 'yes';
cfg.channel ='all'; 
cfg.xlim = 'maxmin';
cfg.ylim = 'maxmin';
cfg.zlim = 'maxmin';
cfg.marker = 'on';
cfg.title = "Frequency: 1:0.1:40. Order: 5";
%cfg.title = "Subj_" + subjID + "_" + "all_channels" ...
    %+ "_" + "singleplot" + frq_str + baseline_str + trials_str;


ft_singleplotTFR(cfg, wavelet_);

% Save the plot
% Needs: Plot name - spritntf with the specifications, file directory
baseline_str = sprintf("_%0.1f_%0.1f_", cfg.baseline(1,1), cfg.baseline(1,2));
trials_str = "trials_" + cfg.trials;


% plot name

if cfg.channel == "all"
    
    plot_name = "multiplicative_09.03_superlet_Subj_" + subjID + "_" + "all_channels" ...
    + "_" + "singleplot" + frq_str + baseline_str + trials_str + ".png";
else

    plot_name = "multiplicative_09.03_superlet_Subj_" + subjID + "_" + string(wavelet_.label(1, cfg.channel)) ...
    + "_" + "singleplot" + frq_str + baseline_str + trials_str + ".png";
end

plot_fullname = fullfile(save_plots_here, plot_name);

plot = gcf;
exportgraphics(plot, plot_fullname);

% cfg = [];
% cfg.baseline     = [-0.5 -0.3];
% cfg.baselinetype = 'absolute';
% cfg.marker       = 'on';
% cfg.layout       = 'easycapM10.mat';
% cfg.channel      = '1';
% cfg.interactive  = 'no';

%% Plot topography

cfg = [];
cfg.parameter = 'powspctrm';
cfg.baseline     = [-0.1 0];
cfg.baselinetype = 'absolute';
cfg.xlim = [min(data_FT.time{1,1}), max(data_FT.time{1,1})]; % Time
cfg.ylim = [1 29]; % Frequency
cfg.channel = 'all';
%cfg.channel = {'all', '-T8', '-C6', '-C4', '-C2', '-FC6', '-FC4', '-FC2', '-FT8'};
cfg.colorbar = 'yes';
cfg.title = "Frequency: 1:0.1:40. Order: 5";

ft_topoplotTFR(cfg, wavelet_);

% Save the plot
% Needs: Plot name - spritntf with the specifications, file directory
baseline_str = sprintf("_%0.1f_%0.1f_", cfg.baseline(1,1), cfg.baseline(1,2));
%trials_str = "trials_" + cfg.trials;


% plot name

if cfg.channel == "all"
    
    plot_name = "multiplicative_09.03_superlet_Subj_" + subjID + "_" + "all_channels" ...
    + "_" + "topography" + frq_str + baseline_str + trials_str + ".png";
else

    plot_name = "multiplicative_09.03_superlet_Subj_" + subjID + "_" + string(wavelet_.label(1, cfg.channel)) ...
    + "_" + "topography" + frq_str + baseline_str + trials_str + ".png";
end

plot_fullname = fullfile(save_plots_here, plot_name);

plot = gcf;
exportgraphics(plot, plot_fullname);
