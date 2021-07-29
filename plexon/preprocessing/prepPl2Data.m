function plexonStructure = prepPl2Data(dataPath, savePath)

plexonStructure = [];
prestimTime = 1;
% dataPath = 'G:\OneDrive - UAB - The University of Alabama at Birmingham\Victor\sorted\Victor_12_22_2020 1_4_MS.pl2';
%% Get metadata
plexonStructure = PL2GetFileIndex(dataPath);


% convert file time into useful date
plexonStructure = convertPlexonDate(plexonStructure); 
%% Get Spike Data

% for all spike channels
sortedCellCount = zeros(plexonStructure.TotalNumberOfSpikeChannels,1);
for i = 1:plexonStructure.TotalNumberOfSpikeChannels
    % for all units in each channel
    for x = 1:plexonStructure.SpikeChannels{i}.NumberOfUnits
    [plexonStructure.SpikeChannels{i}.spike{x}] = PL2Waves(dataPath, i, x); % get spike data for each sorted channel/unit
    
    sortedCellCount(i)= plexonStructure.SpikeChannels{i}.NumberOfUnits; % keep count of which channel have what no. of units
    end
end

plexonStructure.sortedCellCount = sortedCellCount;

% %% Get events
% stimON_Events = PL2EventTs(dataPath, 'EVT05');
% stimOFF_Events = PL2EventTs(dataPath, 'EVT09');


% %% Split up spikes into trials
% count =0;
% for z = 1:length(plexonStructure.sortedCellCount)
%    % if channel has sorted cells
%     if plexonStructure.sortedCellCount(z) > 0
%         count = count+1;
%        % run through sorted cells
%        for d = 1:plexonStructure.sortedCellCount(z)
%            % for each trial
%            for trialNo = 2:length(stimON_Events.Ts)
%                spikesTrial{count,d}{trialNo-1,1} =  plexonStructure.SpikeChannels{z}.spike{d}.Ts(plexonStructure.SpikeChannels{z}.spike{d}.Ts > stimON_Events.Ts(trialNo) - prestimTime  &  plexonStructure.SpikeChannels{z}.spike{d}.Ts < stimOFF_Events.Ts(trialNo));
%            end
%        end
%     end
% end
% 
% plexonStructure.spikesTrial = spikesTrial;

%% Get LFP Data

% for all LFP channels
plexonStructure.TotalNumberOfLFPChannels = length(sortedCellCount)/2;
for q = 1:plexonStructure.TotalNumberOfLFPChannels*2
    [plexonStructure.AnalogChannels{q}.trace] = PL2Ad(dataPath, q); % get CSC data for each  channel
end

end