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



function outputTable = findcommon(inputStruct, maxAllowSep)

nFiles     = length(inputStruct);
nPotSrcs   = 0;
dualLength = 0;
for i = 1:nFiles
    nPotSrcs = nPotSrcs + length(inputStruct(i).X_WORLD);
    if isfield(inputStruct,'dualPhotometry')
        dualLength = dualLength + length(inputStruct.dualPhotometry);
    end
end

% First, we need to know if sources from other bands overlap. We know that,
% from the dual mode, all the sources will overlap. Therefore, we should
% only need to cross check the single mode results
% pre-allocate  ... [sources, fields, bands]
fields    = fieldnames(inputStruct);
dataArray = zeros(nPotSrcs,length(fieldnames(inputStruct)),nFiles+dualLength);
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
        
        % indI(k) and indJ(k) are the two indices where we find a common
        % source between the two bands
        [indI,indJ] = find(dTh <= maxAllowSep);
        
    end
    
end



end