function plexonStructure =  loadPTPEventData(dataPTB, plexonStructure)

load(dataPTB);
codes = prairieCodes();

% load in events from PTB
PTBevents = cellfun(@str2double,stimCmpEvents(2:end,:));

paramStartPosition = find(PTBevents(:,2) ==stringEvent2Num('PARAM_START', codes));

for i =1:length(paramStartPosition)
    
    % grab cnd/ event identity from PTB codes
    block(i,:) = PTBevents(paramStartPosition(i)+1,2); % set the block (param start +1) to PTB block code
    cnd(i,:) = PTBevents(paramStartPosition(i)+2,2); % set the condition (param start +2) to PTB condition code
end

% get total number of each condition presented
cndTotal = zeros(max(cnd),1);
for i=1:max(cnd)
    cndTotal(i) = length(cnd(cnd(:)==i));
end



% transfer into structure
plexonStructure.PTBevents = PTBevents;
plexonStructure.block = block;
plexonStructure.cnd = cnd;
plexonStructure.cndTotal = cndTotal;

end