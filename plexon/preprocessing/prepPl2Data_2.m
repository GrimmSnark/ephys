function plexonStructure = prepPl2Data_2(dataPath)
% This function does the basic read in and preprocessing of plexon pl2 data
% including getting spike data and calculating LFP data based on the
% butterworth settings below
%
% Input- dataPath: fullfile path to the plexon pl2 file you want to process
%
% Output- plexonStructure- structure containing all the information in the
%                          pl2 file in a matlab format. In short the fields
%                          of particular interest are:
%
%        spikes are in:
%           plexonStructure.SpikeChannels{channelNo}.spike{1,cellNo}.Ts  
%
%        LFP are in:
%           plexonStructure.LFP.tracesLowPass{channelNo}  
%
% Example Usage:
%  plexStrt= prepPl2Data_2(G:\Victor\sorted\Victor_12_22_2020 1_4_MS.pl2');
% 
% Written by MA Savage (8/17/2021)
%% defaults

% Butterworth LFP settings
Fs=1000; %Sampling frequency
n=3; %Controls the order of the filter
cuttOffFreq = 200; % Cuttoff freq in Hz

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

%% Get LFP Data

% for all LFP channels
plexonStructure.TotalNumberOfLFPChannels = length(sortedCellCount)/2;
for q = 1:plexonStructure.TotalNumberOfLFPChannels*2
    [plexonStructure.AnalogChannels{q}.trace] = PL2Ad(dataPath, q); % get CSC data for each  channel
end

% reconstruct time vector for CSC
timeVect = (1: plexonStructure.AnalogChannels{1, 1}.trace.FragCounts) / plexonStructure.TimestampFrequency;

% downsample timeVector to 1000Hz
downsampleInterval = plexonStructure.TimestampFrequency/1000;
timeVectDownSam = downsample(timeVect, downsampleInterval);

% downsample CSCs to 1000Hz
for w = 1:plexonStructure.TotalNumberOfLFPChannels
    plexonStructure.LFP.traces{w} = downsample(plexonStructure.AnalogChannels{w}.trace.Values,downsampleInterval);
end

% butterworth cut off 
[b,a]=butter(n,cuttOffFreq/Fs,'low'); %Low-pass filter

% apply filter to LFP channels
for w = 1:plexonStructure.TotalNumberOfLFPChannels
    plexonStructure.LFP.tracesLowPass{w} = filtfilt(b,a, plexonStructure.LFP.traces{w});
end

% transfer params into structure
plexonStructure.LFP.timeVect = timeVect;
plexonStructure.LFP.timeVectDownSam = timeVectDownSam;
plexonStructure.LFP.downsampleInterval = downsampleInterval;

plexonStructure.LFP.butterParams.Fs = Fs;
plexonStructure.LFP.butterParams.order = n;
plexonStructure.LFP.butterParams.cuttOffFreq = cuttOffFreq;


end