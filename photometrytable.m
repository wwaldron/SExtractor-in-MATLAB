function varargout = photometrytable( photoStruct, fileName, varargin )
%PHOTOMETRYTABLE Summary of this function goes here
%   Detailed explanation goes here


% Parse the inputs
prsr = inputParser;
prsr.addRequired('photoStruct', @(x) assert(isstruct(x) && all(isfield(x,{'X_WORLD','Y_WORLD'}))) );
prsr.addRequired('fileName',    @(x) assert(ischar(x)) );
prsr.addParameter('MaxSep',0.2, @(x) assert(isnumeric(x) && isscalar(x) && x > 0) );
prsr.parse(photoStruct,fileName,varargin{:});

% Find Common Sources and combine
dataTable = findcommon(photoStruct,prsr.Results.MaxSep);


end



function outputStruct = findcommon(inputStruct, maxAllowSep)

% Initial Data
nFiles     = length(inputStruct);
nPotSrcs   = 0;
dualFields = 0;
for i = 1:nFiles
    nPotSrcs = nPotSrcs + length(inputStruct(i).X_WORLD);
    inputStruct(i).commonTo = uint16(zeros(size(inputStruct(i).X_WORLD)));
    if isfield(inputStruct(i),'dualPhotometry')
        dualFields = dualFields + length(inputStruct(i).dualPhotometry);
    end
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
        commonInd = (dTh <= maxAllowSep);
        
        % Separate the two bands
        commonI = any(commonInd,1);
        commonJ = any(commonInd,2);
        
        % Use bit-wise addition to indicate cross-correlation
        inputStruct(i).commonTo(commonI(:)) = 2^j + inputStruct(i).commonTo(commonI(:));
        inputStruct(j).commonTo(commonJ(:)) = 2^i + inputStruct(j).commonTo(commonJ(:));
        
    end
end

% Clear up the memory
clearvars commonI commonJ commonInd dDec dRa dTh decBar decI decJ raI raJ

% Now create the RA/DEC arrays
RightAscension = inputStruct(1).X_WORLD;
Declination    = inputStruct(1).Y_WORLD;
for i = 2:nFiles
    % Initialize the common array to false (essentially making each one
    % unique
    common = false(size(inputStruct(i).X_WORLD));
    
    % Look at each preceding band in a bitwise fashion to indicate which
    % sources are no longer unique (sources that are in preceding bands)
    for j = 1:(i-1)
        common = (common | logical(bitand(inputStruct(i).commonTo,2^j)));
    end
    RightAscension = [RightAscension; inputStruct(i).X_WORLD(~common)]; %#ok<*AGROW>
    Declination    = [Declination;    inputStruct(i).Y_WORLD(~common)];
end

% Initialize the table we will use to store in our data structure
tmpTable = table(RightAscension,Declination);
for i = 1:nFiles
    tmpTable.(inputStruct(i).filter) = zeros(size(Declination));
    for j = 1:length(inputStruct(i).dualPhotometry)
        fieldName = [inputStruct(i).filter,'_',inputStruct(i).dualPhotometry(j).filter];
        tmpTable.(fieldName) = zeros(size(Declination));
    end
end

% Finally, create the output structure
goodFields = {'FLUX_ISO','MAG_ISO','FLUX_ISOCOR','MAG_ISOCOR','FLUX_AUTO',...
    'MAG_AUTO','FLUX_BEST','MAG_BEST'};
for i = 1:length(goodFields)
    
    % PLACE HOLDER
    outputStruct.(goodFields{i}) = 0;
    
end



end