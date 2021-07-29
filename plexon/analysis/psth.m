function psth(spikesPerTrial,binSize)

if nargin>2 || isempty(binSize)
    binSize = 10; % default 10ms bin width
end

% find end of trial time
endOfTrial = max(spikesPerTrial{:}); % finds trial end in s

% compute plotting times
lastBin = ceil((post_time*1000)); % last bin edge in ms
edge = EdgeCalculator(bw,1,lastBin); % extract edges for psth
xmin = edge(1); % for plotting
xmax = edge(end); % for plotting


end