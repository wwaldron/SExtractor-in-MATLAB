function varargout = photometrytable( photoStruct, fileName, varargin )
%PHOTOMETRYTABLE Summary of this function goes here
%   Detailed explanation goes here


% Parse the inputs
prsr = inputParser;
prsr.addRequired('photoStruct', @(x) assert(isstruct(x) && all(isfield(x,{'X_WORLD','Y_WORLD'}))) );
prsr.addRequired('fileName',    @(x) assert(ischar(x)) );
prsr.addParameter('MaxSep',0.2, @(x) assert(isnumeric(x) && isscalar(x) && x > 0) );
prsr.addParameter('FilterTable',false, @(x) assert(isscalar(x) && islogical(x)));
prsr.parse(photoStruct,fileName,varargin{:});

% Find Common Sources and combine
dataStructOfTables = findcommon(photoStruct,prsr.Results.MaxSep);

% Filter the table if necessary
if prsr.Results.FilterTable
    filtInd = all(table2array(dataStructOfTables.FLUX_ISO),2);
    fieldNames = fieldnames(dataStructOfTables);
    for i = 1:length(fieldNames)
        tmpTable = dataStructOfTables.(fieldNames{i});
        tmpTable = tmpTable(filtInd,:);
        dataStructOfTables.(fieldNames{i}) = tmpTable;
    end
end

% Write Data to File
writedatatables(dataStructOfTables,fileName);

% Output desired data
if nargout == 1
    varargout{1} = dataStructOfTables;
end


end



function outputStruct = findcommon(inputStruct, maxAllowSep)

% Initial Data
nFiles     = length(inputStruct);
nPotSrcs   = 0;
% dualFields = 0;
for i = 1:nFiles
    nPotSrcs = nPotSrcs + length(inputStruct(i).X_WORLD);
    inputStruct(i).commonTo = zeros(length(inputStruct(i).X_WORLD), ...
        length(inputStruct),'uint32');
%     if isfield(inputStruct(i),'dualPhotometry')
%         dualFields = dualFields + length(inputStruct(i).dualPhotometry);
%     end
end

% First, we need to know if sources from other bands overlap. We know that,
% from the dual mode, all the sources will overlap. Therefore, we should
% only need to cross check the single mode results
for i = 1:nFiles % Over waveband
    for j = (i+1):nFiles % over other wavebands
        
        % For the two current wavebands, get a grid of the right ascensions
        % and declinations.
        [raI,raJ]   = meshgrid(inputStruct(i).X_WORLD, inputStruct(j).X_WORLD);
        [decI,decJ] = meshgrid(inputStruct(i).Y_WORLD, inputStruct(j).Y_WORLD);
        
        % Get the difference and means of the angular distances.
        dRa    = (raI - raJ)*3600;
        dDec   = (decI - decJ)*3600;
        decBar = (decI + decJ)/2;
        
        % Angular distance equation
        dTh = sqrt((dRa.*cosd(decBar)).^2 + dDec.^2);
        
        % 2-D Array that has the indices of the sources the two current
        % bands have in common.
        [commonIndJ,commonIndI] = find(dTh <= maxAllowSep);
        
        % Use bit-wise addition to indicate cross-correlation
        inputStruct(i).commonTo(commonIndI,j) = commonIndJ;
        inputStruct(j).commonTo(commonIndJ,i) = commonIndI;
        
    end
    
    if i == 1
        % create the RA/DEC arrays
        RightAscension = inputStruct(1).X_WORLD;
        Declination    = inputStruct(1).Y_WORLD;
        
        % Since this is the first to be added, each source in the first
        % filter will line up with the same index in the combined array.
        inputStruct(i).commonCombined = (1:length(inputStruct(i).X_WORLD))';
        arrLens = 0;
    else
        % Add sources that are not common to the previous bands
        tmpCommonTo = inputStruct(i).commonTo(:,1:i-1);
        common = any(tmpCommonTo, 2);
        arrLens        = [arrLens, length(RightAscension)];
        RightAscension = [RightAscension; inputStruct(i).X_WORLD(~common)]; %#ok<*AGROW>
        Declination    = [Declination;    inputStruct(i).Y_WORLD(~common)];
        
        % Add the indices of the current filter that match the combined
        % array
        nNew = sum(~common);
        inputStruct(i).commonCombined(~common) = (arrLens(i)+1):(arrLens(i)+nNew);
        inputStruct(i).commonCombined = inputStruct(i).commonCombined';    % Don't know why but this needs transposing here.
        beenAdded = ~common;
        for j = 1:i-1 % Go through the columns to add indices from previous 
            curCommon = tmpCommonTo(:,j);
            inputStruct(i).commonCombined(~beenAdded & curCommon) = ...
                curCommon(~beenAdded & curCommon) + arrLens(j);
            beenAdded = (beenAdded | curCommon);
        end
        
    end
    
end

% Clear up the memory (may not need)
clearvars commonIndI commonIndJ dDec dRa dTh decBar decI decJ raI raJ tmpCommonTo

% Create the output structure
goodFields = {'FLUX_ISO','MAG_ISO','FLUX_ISOCOR','MAG_ISOCOR','FLUX_AUTO',...
    'MAG_AUTO','FLUX_BEST','MAG_BEST'};
tmpTable = table(RightAscension,Declination);
for k = 1:length(goodFields) % Loop through the data we want
    
    for i = 1:nFiles % Loop through the single detection fields
        
        dataToWrite = zeros(size(RightAscension));
        dataToWrite(inputStruct(i).commonCombined) = inputStruct(i).(goodFields{k});
        tmpTable.(inputStruct(i).filter) = dataToWrite;
        
        for j = 1:length(inputStruct(i).dualPhotometry) % Loop through dual detection fields
            
            dataToWrite = zeros(size(RightAscension));
            dataToWrite(inputStruct(i).commonCombined) = inputStruct(i).dualPhotometry(j).(goodFields{k});
            tableFieldName = [inputStruct(i).filter,'_',inputStruct(i).dualPhotometry(j).filter];
            tmpTable.(tableFieldName) = dataToWrite;
            
        end
        
    end
    
    outputStruct.(goodFields{k}) = tmpTable;
    
end

end



% Function to write the data to a spreadsheet format
function writedatatables(structOfTables,fileName)

fieldNames = fieldnames(structOfTables);

origWarnState = warning('off');
for i = 1:length(fieldNames)
    
    tblToWrite = structOfTables.(fieldNames{i});
    writetable(tblToWrite, fileName, 'FileType', 'spreadsheet',...
        'Sheet', fieldNames{i});
    
end
warning(origWarnState);

end



