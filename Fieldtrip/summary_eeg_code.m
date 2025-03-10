%% SUMMARY EEG TABLES

% This script is written with Eva Marie von Butler

% Follow the instructions. They need to be placed to the different
% sections of the script.

% POz shows the best activation. But if the mean of all channels is needed, 
% it is included in the code. It can be taken out as well.
% Or, once a channel is located, then the code can be customized for that 
% channel too.
% Then the variables should be created.
% 
% This code uses the chanlocs from the EEG data, a very manual approach.
% 23 is P7
% 27 is P8
% 28 is POz --- currently used
% 62 is PO7
% 63 is PO8
% 64 is PO7

%--------------------------------------------------------------------------

% Place this section under the subjects array
% 
% HERE: subjects = [1,2,5,12,15,16,19,20,21,22,27,29,32,33,34,37,38];%,42,43,44,46];

% Initiate empty tables for the Summary_EEG and fill them

% All channel info
summary_eeg_head_all = zeros(length(subjects), 1050); % 1050 bcs of the time points we have
summary_eeg_body_all = zeros(length(subjects), 1050);
summary_eeg_bground_all = zeros(length(subjects), 1050);


    
% Only for POz
summary_eeg_head_poz = zeros(length(subjects), 1050);
summary_eeg_body_poz = zeros(length(subjects), 1050);
summary_eeg_bground_poz = zeros(length(subjects), 1050);

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% 
% This part is placed right after the separation of the events. 
% Epoch, separate to EEG_head etc. then next step is the following code

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
        
%--------------------------------------------------------------------------


%% PLOT THE ERPs
% 
% 
% We plot the first one to set the x andy labels as well as the title,
% then while the figure is hold on, we plot the rest.

% HEAD
mean_head_poz_allsubj = mean(summary_eeg_head_poz, 1);
figure;
plot(EEG_head.times, summary_eeg_head_poz(1,:), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.9)
title('Activation at POz in HEAD condition');
%xlim([-200 500]);
xlabel('Time (ms)');
ylabel('Potential(μV)');
hold on


for sub = 2:size(summary_eeg_head_poz,1)
    
    plot(EEG_head.times, summary_eeg_head_poz(sub,:), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.9)
    plot(EEG_head.times, mean_head_poz_allsubj, 'Color', [0 0.4470 0.7410 1], 'LineWidth', 1.2)
    
    %xlim([-200 500]);
    
end

% BODY

mean_body_poz_allsubj = mean(summary_eeg_body_poz, 1);
figure;
plot(EEG_head.times, summary_eeg_body_poz(1,:), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.9)
title('Activation at POz in BODY condition');
%xlim([-200 500]);
xlabel('Time (ms)');
ylabel('Potential(μV)');
hold on


for sub = 2:size(summary_eeg_body_poz,1)
    
    plot(EEG_head.times, summary_eeg_body_poz(sub,:), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.9)
    plot(EEG_head.times, mean_body_poz_allsubj, 'Color', [0 0.4470 0.7410 1], 'LineWidth', 1.2)
    
    %xlim([-200 500]);
    
end

% BACKGROUND

mean_bground_poz_allsubj = mean(summary_eeg_bground_poz, 1);
figure;
plot(EEG_head.times, summary_eeg_head_poz(1,:), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.9)
title('Activation at POz in BACKGROUND condition');
%xlim([-200 500]);
xlabel('Time (ms)');
ylabel('Potential(μV)');
hold on


for sub = 2:size(summary_eeg_bground_poz,1)
    
    plot(EEG_head.times, summary_eeg_bground_poz(sub,:), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.9)
    plot(EEG_head.times, mean_bground_poz_allsubj, 'Color', [0 0.4470 0.7410 1], 'LineWidth', 1.2)
    
    %xlim([-200 500]);
    
end

