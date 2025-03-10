% Not a function.

% Purpose of the script:
% 
% This script is run in between the analysis steps of the main analysis
% script.
% Without looping, we add the extracted subject data to one single array
% The variables created here will be used in the main script.



% Concatenate the subjects all together in one array


% Create empty arrays. We overwrite the data to the zeros. The number of
% time points are always the same due to epoching, that's why the
% EEG_head.times are used as the number of columns throughout the arrays.


% all_subj_head_po7 = zeros(3, length(EEG_head.times));
% all_subj_head_po8 = zeros(3, length(EEG_head.times));
% 
% all_subj_body_po7 = zeros(3, length(EEG_head.times));
% all_subj_body_po8 = zeros(3, length(EEG_head.times));
% 
% all_subj_bground_po7 = zeros(3, length(EEG_head.times));
% all_subj_bground_po8 = zeros(3, length(EEG_head.times));



%%  sdfdffds
% Run each section after their preprocessing completed

%% RUN FIRST
% Subj1 PO7
EEG_subj1_head_PO7 = EEG_head_PO7;
EEG_subj1_body_PO7 = EEG_body_PO7;
EEG_subj1_bground_PO7 = EEG_bground_PO7;

% Subj1 PO8
EEG_subj1_head_PO8 = EEG_head_PO8;
EEG_subj1_body_PO8 = EEG_body_PO8;
EEG_subj1_bground_PO8 = EEG_bground_PO8;

%% RUN SECOND
% Subj2 PO7
EEG_subj2_head_PO7 = EEG_head_PO7;
EEG_subj2_body_PO7 = EEG_body_PO7;
EEG_subj2_bground_PO7 = EEG_bground_PO7;

% Subj2 PO8
EEG_subj2_head_PO8 = EEG_head_PO8;
EEG_subj2_body_PO8 = EEG_body_PO8;
EEG_subj2_bground_PO8 = EEG_bground_PO8;

%% RUN THIRD

% Subj03 PO7
EEG_subj3_head_PO7 = EEG_head_PO7;
EEG_subj3_body_PO7 = EEG_body_PO7;
EEG_subj3_bground_PO7 = EEG_bground_PO7;

% Subj03 PO8

EEG_subj3_head_PO8 = EEG_head_PO8;
EEG_subj3_body_PO8 = EEG_body_PO8;
EEG_subj3_bground_PO8 = EEG_bground_PO8;


 


%% MERGE

EEG_merged_head_PO7 = pop_mergeset(EEG_subj1_head_PO7,EEG_subj2_head_PO7, 1);
EEG_merged_head_PO8 = pop_mergeset(EEG_subj1_head_PO8,EEG_subj2_head_PO8, 1);

EEG_merged_body_PO7 = pop_mergeset(EEG_subj1_body_PO7,EEG_subj2_body_PO7, 1);
EEG_merged_body_PO8 = pop_mergeset(EEG_subj1_body_PO8,EEG_subj2_body_PO8, 1);


EEG_merged_bground_PO7 = pop_mergeset(EEG_subj1_bground_PO7,EEG_subj2_bground_PO7, 1);
EEG_merged_bground_PO8 = pop_mergeset(EEG_subj1_bground_PO8,EEG_subj2_bground_PO8, 1);
