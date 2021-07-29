function processVictorMelanpsin(dataPath)

%% extract the data
plexonStructure = prepPl2Data(dataPath);

%% get the matchnig PTB files
[parentDirectory, filename]= returnParentFolder(dataPath);
PTB_Folder = dir([parentDirectory '\PTB*']);

% get PTB event folders
eventFolders = dir([PTB_Folder.folder '\' PTB_Folder.name '\*.mat']);

% get time date strings from all files
for i =1:length(eventFolders)
    tempString = split(eventFolders(i).name, '_');
    timeDates{i} = tempString{end}(1:end-4);
end


uniqueTimeString = unique(timeDates);

% split filepaths into unique pairs
for x = 1:length(uniqueTimeString)
    indexOfCells = strfind({eventFolders.name}, uniqueTimeString(x));
    index = find(~cellfun(@isempty, indexOfCells));
    
    for q = 1:length(index)
        pairedFilepaths{x}{q} = fullfile(eventFolders(index(q)).folder, eventFolders(index(q)).name);
    end
end

% trim to only include entries that have pairs (ie both
% stimParams and events
count = 1;
for c = 1:length(pairedFilepaths)
    if length(pairedFilepaths{c})>1
        pairedFilepathsTrimmed{count} = pairedFilepaths{c};
        uniqueTimeStringTrimmed(count) = str2double(uniqueTimeString{c});
        count = count +1;
    end
end

% find closest time to event files
timeDiff = uniqueTimeStringTrimmed - plexonStructure.DateTime  ;

timeDiff(timeDiff<0) = NaN;
[~, matchIndx]= min(timeDiff);

PTBFilePaths = pairedFilepaths{matchIndx};

% load in stim parameters
load(PTBFilePaths{2});
plexonStructure.stimParams = stimParams;

plexonStructure =  loadPTPEventData(PTBFilePaths{1}, plexonStructure);

%% Get events
stimON_Events = PL2EventTs(dataPath, 'EVT05');
stimOFF_Events = PL2EventTs(dataPath, 'EVT09');

% remove identical event (can be the first one due to serial port flush)
[~, match ] = intersect(stimON_Events.Ts, stimOFF_Events.Ts);

stimON_Events.Ts(match) = [];
stimOFF_Events.Ts(match) = [];

%% Filter LFPs

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

Fs=1000; %Sampling frequency
n=4; %Controls the order of the filter
cuttOffFreq = 200; % Cuttoff freq in Hz
[b,a]=butter(n,cuttOffFreq/(Fs/2)); %Low-pass filter

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

%% Split up spikes into trials
count =0;
prestimTime = plexonStructure.stimParams.preStimTime;
poststimTime = plexonStructure.stimParams.postStimTime;

for z = 33:length(plexonStructure.sortedCellCount)
    % if channel has sorted cells
    if plexonStructure.sortedCellCount(z) > 0
        count = count+1;
        % run through sorted cells
        for d = 1:plexonStructure.sortedCellCount(z)
            
            % for each trial
            for trialNo = 1:length(stimON_Events.Ts)
                spikesTrial{count,d}{trialNo,1} =  plexonStructure.SpikeChannels{z}.spike{d}.Ts(plexonStructure.SpikeChannels{z}.spike{d}.Ts > stimON_Events.Ts(trialNo) - prestimTime  &  plexonStructure.SpikeChannels{z}.spike{d}.Ts < stimOFF_Events.Ts(trialNo) + poststimTime)';
                
                if ~isempty(spikesTrial{count,d}{trialNo,1})
                    spikesTrialZeroed{count,d}{trialNo,1} = spikesTrial{count,d}{trialNo,1} -  spikesTrial{count,d}{trialNo,1}(1); % zero the spike trial
                else
                    spikesTrialZeroed{count,d}{trialNo,1} = 0; % zero the spike trial
                end
            end
        end
    end
end

plexonStructure.spikesTrial = spikesTrial;
plexonStructure.spikesTrialZeroed = spikesTrialZeroed;

%% Split spike trials into conditions

% for each spike channel
for z = 1:length(spikesTrial)
    % if channel has sorted cells
    if ~isempty(plexonStructure.spikesTrial(z))
        % run through sorted cells
        for d = 1:plexonStructure.sortedCellCount(z+32)
            if ~isempty(spikesTrialZeroed{z,d})
                % for each condition
                for cc = 1:max(plexonStructure.cnd)
                    % for each trial of that condition
                    cndTrials = find(plexonStructure.cnd == cc);
                    
                    blueTrials(cc,:) = cndTrials(rem(cndTrials,2)~=0);
                    redTrials(cc,:) = cndTrials(rem(cndTrials,2)==0);
                    
                    spikesTrialZeroedPerCnd{z,d}{cc,1} = {spikesTrialZeroed{z,d}{blueTrials(cc,:),1}};
                    spikesTrialZeroedPerCnd{z,d}{cc,2} = {spikesTrialZeroed{z,d}{redTrials(cc,:),1}};
                end
            end
        end
    end
end

plexonStructure.spikesTrialZeroedPerCnd = spikesTrialZeroedPerCnd;

%% Split LFP into trials

prestimTime = plexonStructure.stimParams.preStimTime;
poststimTime = plexonStructure.stimParams.postStimTime;

% get the trial length in samples
for trialNo = 1:length(stimON_Events.Ts)
    timeStart = stimON_Events.Ts(trialNo) - prestimTime;
    timeStop = stimOFF_Events.Ts(trialNo) + poststimTime;
    
    startIn = find(plexonStructure.LFP.timeVectDownSam > timeStart, 1, 'first');
    stopIn = find(plexonStructure.LFP.timeVectDownSam < timeStop, 1, 'last');
    
    trialLengthSamples(trialNo) = stopIn - startIn  +1;
end

% take mode value for trail extraction 
trialLength = mode(trialLengthSamples);

for z = 1:plexonStructure.TotalNumberOfLFPChannels
     
    % for each trial
    for trialNo = 1:length(stimON_Events.Ts)
        timeStart = stimON_Events.Ts(trialNo) - prestimTime;
        timeStop = stimOFF_Events.Ts(trialNo) + poststimTime;
        
        startIn = find(plexonStructure.LFP.timeVectDownSam > timeStart, 1, 'first');
        stopIn = startIn + trialLength-1;
        
        LFP_PerTrial{z}(:,trialNo) = plexonStructure.LFP.tracesLowPass{z}(startIn:stopIn); 
     end
end


% for z = 33:plexonStructure.TotalNumberOfLFPChannels*2
%      
%     % for each trial
%     for trialNo = 1:length(stimON_Events.Ts)
%         timeStart = stimON_Events.Ts(trialNo) - prestimTime;
%         timeStop = stimOFF_Events.Ts(trialNo) + poststimTime;
%         
%         startIn = find(plexonStructure.LFP.timeVectDownSam > timeStart, 1, 'first');
%         stopIn = startIn + trialLength-1;
%         
%         LFP_PerTrial{z}(:,trialNo) = plexonStructure.AnalogChannels{z}.trace.Values(startIn:stopIn); 
%      end
% end


plexonStructure.LFP.LFP_PerTrial = LFP_PerTrial;

%% Split LFP into conditions

for z = 1:plexonStructure.TotalNumberOfLFPChannels
%  for z = 1:plexonStructure.TotalNumberOfLFPChannels*2
    % for each condition
    for cc = 1:max(plexonStructure.cnd)
        % for each trial of that condition
        LFP_PerCnd{z}(:,:,cc,1) = LFP_PerTrial{z}(:,blueTrials(cc,:));
        LFP_PerCnd{z}(:,:,cc,2) = LFP_PerTrial{z}(:,redTrials(cc,:));
    end
end

plexonStructure.LFP.LFP_PerCnd = LFP_PerCnd; % sample x reptition x stim level x LED 

%% plot the data

 saveDir = [parentDirectory '\figs\' filename{:}(1:end-4) '\'];

if ~exist(saveDir)
    mkdir(saveDir);
end
% 
ColorText{1} = 'Blue LED';
ColorText{2} = 'Red LED';
trialLength = plexonStructure.stimParams.preStimTime + plexonStructure.stimParams.stimTime +plexonStructure.stimParams.postStimTime;

for ff = 1:length(spikesTrial)
    % if channel has sorted cells
    if ~isempty(plexonStructure.spikesTrial(ff))
        % run through sorted cells
        for d = 1:plexonStructure.sortedCellCount(ff+32)
            if ~isempty(spikesTrialZeroed{ff,d})
                figHandle= figure('units','normalized','outerposition',[0 0 0.5 1]);
                
                count = 1;
                for ww = 1:size(spikesTrialZeroedPerCnd{ff,d},1) % for each stim level
                    for ee = 1:size(spikesTrialZeroedPerCnd{ff,d},2) % for LED colour
                        
                        axH(count) = subplot(size(spikesTrialZeroedPerCnd{ff,d},1) , size(spikesTrialZeroedPerCnd{ff,d},2), count);
                        plotPSTH(axH(count),spikesTrialZeroedPerCnd{ff,d}{ww,ee}, [0 trialLength], [prestimTime, prestimTime+plexonStructure.stimParams.stimTime], prestimTime,500 ,5 );
                        
                        title([ColorText{ee} ' Stim Level ' num2str(ww)]);
                        count = count+1;
                    end
                end
                
                axYLim =  get(axH,'YLim');
                maxY = max([axYLim{:}]);
                set(axH,'YLim',[0 maxY]);
                
                axYLabels = get(axH,'YTickLabel');
                [~, maxIndex] = max(cellfun('size', axYLabels, 1));
                set(axH,'YTickLabel',axYLabels{maxIndex});
                
                vline2(axH,[prestimTime, prestimTime+plexonStructure.stimParams.stimTime],{'r', 'r'}, {'' ''});
                
                tightfig
                saveas(gcf, [saveDir sprintf('Chan %d Cell no. %d', ff, d ) '.png']);
                close
            end
        end
    end
end

%% Spectrograms


ColorText{1} = 'Blue LED';
ColorText{2} = 'Red LED';

for ff = 1:plexonStructure.TotalNumberOfLFPChannels
% for ff = 33:plexonStructure.TotalNumberOfLFPChannels*2
    
    figHandle= figure('units','normalized','outerposition',[0 0 0.7 1]);
    hold on 
    count = 1;
    for ww = 1:size(LFP_PerCnd{ff},3) % for each stim level
        for ee = 1:size(LFP_PerCnd{ff},4) % for LED colour
            
            params.Fs=1000; %Sampling rate
            params.tapers=[3 5]; %Taper parameters
%             params.fpass=[5 60]; %Look at 5-60 Hz band
            params.fpass=[5 140]; %Look at 5-60 Hz band
            params.trialave=1; %Average across trials
            
            
            movingwin = [0.25 0.025];
            
            axH(count) = subplot(size(LFP_PerCnd{ff},3) ,size(LFP_PerCnd{ff},4), count);
            [P2,T2,F2]=mtspecgramc(LFP_PerCnd{ff}(:,:,ww, ee),movingwin, params);
            
            imagesc(axH(count), T2-.5,F2,10*log10(P2')) %Plot power in dB
            colormap('jet');
            axis xy; xlabel('Time(s)'); ylabel('Freq (Hz)');
            colorHandle = colorbar;
            colorHandle.Label.String = 'Power (dB)';
            
            title([ColorText{ee} ' Stim Level ' num2str(ww)]);
            count = count+1;
        end
    end    
    vline2(axH,[prestimTime, prestimTime+plexonStructure.stimParams.stimTime],{'color', [1 0 1] , 'linestyle', '-', 'linewidth', 2; 'color', [1 0 1],  'linestyle', '-', 'linewidth', 2}, {'' ''});
   
    subplotEvenColorBar(figHandle)
    tightfig
    saveas(gcf, [saveDir sprintf('LFP Spect Chan %d', ff ) '.png']);
    close
end