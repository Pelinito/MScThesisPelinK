%% Toy Struct to trick FieldTrip - Timelock Analysis

toy = [];

toy.time = timelockSubj1.time;
toy.label = timelockSubj1.label;
toy.elec = timelockSubj1.elec;
toy.sampleinfo = timelockSubj1.sampleinfo(1:18, :);
toy.trial = zeros(18,64,1050); % should be general: NSubj x elec x time
toy.trialinfo = ones(18,1);
toy.dimord = 'rpt_chan_time';   % Changed in case the functions check for it

for avg=1:NSubj

    toy.trial(avg,:,:) = head_timelock_avg{1,avg}.avg;

end

