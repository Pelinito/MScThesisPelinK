%% Separate the events

% SEPARATEEVENTS - Extracts the given events from the EEG data
%
%
% Usage: [EEG_head, EEG_body, EEG_bground] = separateevents(EEG, min, max);
%
% Input: - EEG: EEG data to separate the events from 
%        - min: the beginning of the interval. 
%        - max: the end of the interval.
     


function [EEG_head, EEG_body, EEG_bground] = separateevents(EEG)

   EEG = eeg_checkset(EEG);

%     if exist(EEG.event) == 0
%     
%         error('EEG events do not exist!');
%      
%     else
%         
%         disp('Events are there :)');
%         
%     end
    
%     min = string(min);
%     max = string(max);
%     
%     interval = string(min)+'<='+string(max);
%     
%     
    % separate the head events
     
    
     EEG_head = pop_selectevent(EEG, 'latency','-0.875<=1.175','type','head','renametype','head','deleteevents','on');
     EEG_head = eeg_checkset(EEG_head);
    
     EEG_body = pop_selectevent(EEG, 'latency','-0.875<=1.175','type','body','renametype','body','deleteevents','on');
     EEG_body = eeg_checkset(EEG_body);
     
     EEG_bground = pop_selectevent(EEG, 'latency','-0.875<=1.175','type','bgrd','renametype','bground','deleteevents','on');
     EEG_bground = eeg_checkset(EEG_bground);

     disp('Separation completed :)');
    
    

end