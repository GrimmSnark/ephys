function plotPSTH(axHandle, spikes, trialXLim, stimTimes, prestimTime, binSize, xtickIncrement)

if nargin < 6 || isempty(binSize)
    binSize = 50; % ms
end

if nargin < 7 || isempty(xtickIncrement)
    xtickIncrement = 5;
end


%% PSTH

all = [];
for iTrial = 1:length(spikes)
    all             = [all; spikes{iTrial}'];               % Concatenate spikes of all trials
end

ax                  = axHandle;

slength             = trialXLim(2) * 1000;                                  % Length of signal trace [ms]
nbins                = slength/binSize;                        % Bin duration in [ms]
nobins              = 1000/binSize;                            % No of bins/sec

h                   = histogram(ax,all,nbins);
h.FaceColor         = 'k';

mVal                = max(h.Values)+round(max(h.Values)*.1);
ax.XLim             = [trialXLim(1) trialXLim(2)];
ax.YLim             = [0 mVal];
ax.XTick            = [trialXLim(1): xtickIncrement : trialXLim(2)]+ prestimTime;
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

end