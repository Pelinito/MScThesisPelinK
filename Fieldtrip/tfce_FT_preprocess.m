%% Preprocess for the FieldTrip

% Start with the Holy Trinity

clear all;
close all;
clc;

%% Define the paths

tfce_path = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE';
plots_tfce = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/Plots';
savesets_here = '/net/store/nbp/projects/wd_ride_village/ERSP/Pelin/TFCE/ConditionSets';

%cd(savesets_here);  % files are loaded from here
%% The MAIN GOAL: Fill these cells

% Give the number of subjects - In our specific case, it is 18. 
% Probably, will not need to change

NSubj = 18;

% allSubjx cell structures to save the data transformed into FT

allSubjHead = cell(1,Nsubj);
allSubjBody = cell(1,Nsubj);
allSubjBgrd = cell(1,Nsubj);
allSubjnonHead = cell(1,Nsubj);



% Struct to read data
%cfg = [];
%cfg.folderlocation = sprintf(savesets_here);
%cfg.filename_head = sprintf()

%head_str = 

%% Read the matching files. dir makes them to be stored in a struct

bgrd_file_info = dir('EEG_epcohed_bgrd*.set');
body_file_info = dir('EEG_epcohed_body*.set');
head_file_info = dir('EEG_epcohed_head*.set');

%% Loop the file names in the struct

cfg = [];


for i=1:length(head_file_info)

    cfg.dataset = [head_file_info.name '/' head_file_info.folder];
    EEG_head = ft_preprocessing(cfg);
    allSubjHead{1,i} = EEG_head;

end


for i=1:length(body_file_info)

    cfg.dataset = [body_file_info.name '/' body_file_info.folder];
    EEG_body = ft_preprocessing(cfg);
    allSubjBody{1,i} = EEG_body;

end

for i=1:length(bgrd_file_info)

    cfg.dataset = [bgrd_file_info.name '/' bgrd_file_info.folder];
    EEG_bgrd = ft_preprocessing(cfg);
    allSubjBgrd{1,i} = EEG_bgrd;

end

%% Append the data for the nonHead condition



