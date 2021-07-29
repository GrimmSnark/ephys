function plotPSTH_Raster(spikes, trialXLim, stimTimes, prestimTime, xtickIncrement)


if nargin < 5 || isempty(xtickIncrement)
    xtickIncrement = 5;
end

binSize = 50; % ms

figure('Units','normalized','Position',[0 0 .3 .5]) % set up figure

%% PSTH

all = [];
for iTrial = 1:length(spikes)
    all             = [all; spikes{iTrial}'];               % Concatenate spikes of all trials
end

ax                  = subplot(2,1,1);

slength             = trialXLim(2) * 1000;                                  % Length of signal trace [ms]
nbins                = slength/binSize;                        % Bin duration in [ms]
nobins              = 1000/binSize;                            % No of bins/sec

h                   = histogram(all,nbins);
h.FaceColor         = 'k';

mVal                = max(h.Values)+round(max(h.Values)*.1);
ax.XLim             = [trialXLim(1) trialXLim(2)];
ax.YLim             = [0 mVal];
ax.XTick            = [trialXLim(1): xtickIncrement : trialXLim(2)];
ax.XTickLabels      = [trialXLim(1): xtickIncrement : trialXLim(2)- prestimTime];
ax.XLabel.String  	= 'Time [s]';
ax.YLabel.String  	= 'Spikes/Bin';

for iLab = 1:length(ax.YTickLabel)
    lab             = str2num(ax.YTickLabel{iLab});
    conv            = (lab / length(spikes)) * nobins; 	% Convert to [Hz]: avg spike count * bins/sec
    newlabel{iLab}  = num2str(round(conv));                 % Change YLabel
end
ax.YTickLabel       = newlabel;
ax.YLabel.String  	= 'Firing Rate [Hz]';

vline2(stimTimes, {'r' 'r'}, {'' ''});


%% Raster

ax = subplot(2,1,2); hold on

% For all trials...
for iTrial = 1:length(spikes)
                  
    spks            = spikes{iTrial};         % Get all spikes of respective trial    
    xspikes         = repmat(spks,3,1);         % Replicate array
    yspikes      	= nan(size(xspikes));       % NaN array
    
    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                % Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end
    
    plot(xspikes, yspikes, 'Color', 'k')
end

ax.XLim             = [trialXLim(1) trialXLim(2)];
ax.YLim             = [0 length(spikes)];
ax.XTick            = [trialXLim(1): xtickIncrement : trialXLim(2)];
ax.XTickLabels      = [trialXLim(1): xtickIncrement : trialXLim(2)- prestimTime];

ax.XLabel.String  	= 'Time [s]';
ax.YLabel.String  	= 'Trials';

vline2(stimTimes, {'r' 'r'}, {'' ''});


end