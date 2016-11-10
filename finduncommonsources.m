function dataStruct = finduncommonsources( dataStruct, varargin )
%FINDUNCOMMONSOURCES Uses SExtractor to find the sources in images that are
%not common to multiple wavebands
%   Often times in visible space astronomy cosmic rays strike the FPA and
%   create bright pixels that mimic true sources. This code calls
%   SExtractor (assumed to be installed already), builds/reads the created
%   catalog, and searches for sources that are not common to each waveband
%   (i.e. potential cosmic rays). The code outputs as many files are in
%   FITSFILES in an (x,y) format using Ra and Dec for reading by ds9. If
%   the config file doesn't exist, the program creates one using the
%   default configuration. Changes to several of the config options can be
%   made using Name/Value pairs. (Does not give an output FITS)


% First validate the inputs
prsr = inputParser;
prsr.StructExpand = true;
prsr.addRequired('dataStruct', ...
    @(x) assert( isstruct(x) && all(isfield(x,{'X_WORLD','Y_WORLD','catFile'})) ));
prsr.addParameter('MaxAngSepTrueSrc', 0.2,    @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0))); % ArcSecond
prsr.addParameter('RegOutputFormat','XYWorld',@(x) assert(any(strcmpi(x,{'XYWorld','XYImage','saoimage','saoimagecirc','saoworld'}))));
prsr.addParameter('CircleRegionRadius', 3,    @(x) assert(isnumeric(x) && isscalar(x) && x > 0));
prsr.addParameter('RmLargeSources', false,    @(x) assert(islogical(x) && isscalar(x)));
prsr.addParameter('LargeSourceArea',  100,    @(x) assert(isnumeric(x) && isscalar(x) && x > 0));
prsr.parse(dataStruct,varargin{:});

% Make scalars into vectors
fileListLen = length(dataStruct);
wrnSt       = warning;
warning('off');
prsr        = struct(prsr);
warning(wrnSt);
prsr.Results.MaxAngSepTrueSrc  = makevector(prsr.Results.MaxAngSepTrueSrc,fileListLen);

% Setup the star ind
nFiles = length(dataStruct);
for i = 1:nFiles
    dataStruct(i).starAndLrgInd = false(size(dataStruct(i).X_WORLD));
end

% Identify and write out the bad sources
for i = 1:nFiles % Over wavebands
    % In the wave band dimension, we know that we get an NxN array
    % corresponding to the pairing of wavebands. However, if we compare a
    % waveband to itself, the dTh will be zero. This will serve no good
    % purpose so we can ignore the diagonal. Also, the matrix is
    % symmetrical (wvb1 x wvb2 = wvb2 x wvb1) therefore we only have to
    % concern ourselves with one side of the diagonal. As you can see from
    % the loop over j, we choose the top half of the matrix. This means
    % that for N=4, instead of having to iterate 16 times, we only have to
    % iterate 6 times.
    % Note that X_WORLD is Right Ascension
    % Note that Y_WORLD is Declination
    
    if prsr.Results.MaxAngSepTrueSrc(i) > 0
        
        for j = (i+1):nFiles % over other wavebands
            
            % For the two current wavebands, get a grid of the right ascensions
            % and declinations.
            [raI,raJ]   = meshgrid(dataStruct(i).X_WORLD, dataStruct(j).X_WORLD);
            [decI,decJ] = meshgrid(dataStruct(i).Y_WORLD, dataStruct(j).Y_WORLD);
            
            % Get the difference and means of the angular distances.
            dRa    = (raI - raJ)*3600;
            dDec   = (decI - decJ)*3600;
            decBar = (decI + decJ)/2;
            
            % Angular distance equation
            dTh = sqrt((dRa.*cosd(decBar)).^2 + dDec.^2);
            
            % 2D Logical Array indicating where the angular separation is less
            % than the maximum allowed separation for a 'true' source.
            curstarAndLrgInd = (dTh <= prsr.Results.MaxAngSepTrueSrc(i));
            
            % Update the logical arrays for each waveband indicating where
            % stars are
            starI = any(curstarAndLrgInd,1);
            starJ = any(curstarAndLrgInd,2);
            dataStruct(i).starAndLrgInd = dataStruct(i).starAndLrgInd | starI(:);
            dataStruct(j).starAndLrgInd = dataStruct(j).starAndLrgInd | starJ(:);
            
        end
        
        % If the source is too large, remove it too
        if prsr.Results.RmLargeSources && ...
                all(isfield(dataStruct,{'A_IMAGE','B_IMAGE'}))
            areas = pi*dataStruct(i).A_IMAGE.*dataStruct(i).B_IMAGE;
            largeInd2Rm = areas > prsr.Results.LargeSourceArea;
            dataStruct(i).starAndLrgInd = dataStruct(i).starAndLrgInd | largeInd2Rm;
        end
        
        % Write those false sources to a file
        [catPath,catFile,~] = fileparts(dataStruct(i).catFile);
        falseRegFile = fullfile(catPath,[catFile,'_falsesource','.reg']);
        writetoregfile(falseRegFile,prsr.Results.RegOutputFormat,...
            dataStruct(i),prsr.Results.CircleRegionRadius);
        
    else
        
        % If the source is too large, remove it too
        if prsr.Results.RmLargeSources && ...
                all(isfield(dataStruct,{'A_IMAGE','B_IMAGE'}))
            areas = pi*dataStruct(i).A_IMAGE.*dataStruct(i).B_IMAGE;
            largeInd2Rm = areas > prsr.Results.LargeSourceArea;
            dataStruct(i).starAndLrgInd = dataStruct(i).starAndLrgInd & ~largeInd2Rm;
        end
        
        % Write those false sources to a file
        [catPath,catFile,~] = fileparts(dataStruct(i).catFile);
        allRegFile = fullfile(catPath,[catFile,'_allsource','.reg']);
        writetoregfile(allRegFile,prsr.Results.RegOutputFormat,...
            dataStruct(i),prsr.Results.CircleRegionRadius);
        
    end
    
end

end


function writetoregfile(file,outFormat,dataStruct,circRad)


% If desired format is SAOImage and the right fields are in place.
if strcmpi(outFormat,'saoimagecirc') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    writeInd = ~dataStruct.starAndLrgInd;
    tmpData = [dataStruct.X_IMAGE(writeInd),dataStruct.Y_IMAGE(writeInd),...
        circRad*ones(size(dataStruct.Y_IMAGE(writeInd)))]';
    fid = fopen(file,'w');
    fprintf(fid,'# format: image\n');
    fprintf(fid,'circle(%.10f,%.10f,%.10f)\n',tmpData);
    fclose(fid);
    
elseif strcmpi(outFormat,'saoimage') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE'}))
    
    writeInd = ~dataStruct.starAndLrgInd;
    tmpData = [dataStruct.X_IMAGE(writeInd),dataStruct.Y_IMAGE(writeInd),...
        dataStruct.A_IMAGE(writeInd),dataStruct.B_IMAGE(writeInd),...
        dataStruct.THETA_IMAGE(writeInd)]';
    fid = fopen(file,'w');
    fprintf(fid,'# format: image\n');
    fprintf(fid,'ellipse(%.10f,%.10f,%.10f,%.10f,%.10f)\n',tmpData);
    fclose(fid);
    
% If desired format is SAOImage but don't have ellipse parameters, but do
% have centers in image space
elseif strcmpi(outFormat,'saoimage') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    warning('User chose SAOImage, but ellipse parameters not provided. Writing out (X,Y)_IMAGE data for %s.',file);
    x = dataStruct.X_IMAGE(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
% If desired format is SAOImage but don't have ellipse parameters but do
% have ceters in world space
elseif strcmpi(outFormat,'saoimage') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    warning('User chose SAOImage, but ellipse parameters not provided. Writing out (X,Y)_WORLD data for %s.',file.');
    x = dataStruct.X_WORLD(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_WORLD(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
% If desired format is SAOWorld and have correct fields
elseif strcmpi(outFormat,'saoworld') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD','A_WORLD','B_WORLD','THETA_WORLD'}))
    
    writeInd = ~dataStruct.starAndLrgInd;
    tmpData = [dataStruct.X_WORLD(writeInd),dataStruct.Y_WORLD(writeInd),...
        dataStruct.A_WORLD(writeInd),dataStruct.B_WORLD(writeInd),...
        dataStruct.THETA_WORLD(writeInd)]';
    fid = fopen(file,'w');
    fprintf(fid,'# format: degrees (fk5)\n');
    fprintf(fid,'+ellipse(%.10f,%.10f,%.10f,%.10f,%.10f) #green\n',tmpData);
    fclose(fid);
    
% If desired format is SAOWorld but missing ellipse params but have centers
% in world coordinates
elseif strcmpi(outFormat,'saoworld') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    warning('User chose SAOWorld, but ellipse parameters not provided. Writing out (X,Y)_WORLD data for %s.',file.');
    x = dataStruct.X_WORLD(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_WORLD(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y]);
    
% If desired format is SAOWorld but missing ellipse params but have centers
% in image coordinates
elseif strcmpi(outFormat,'saoworld') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    warning('User chose SAOWorld, but ellipse parameters not provided. Writing out (X,Y)_IMAGE data for %s.',file.');
    x = dataStruct.X_IMAGE(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
% You get the picture . . .
elseif strcmpi(outFormat,'XYWorld') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    x = dataStruct.X_WORLD(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_WORLD(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
elseif strcmpi(outFormat,'XYWorld') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    warning('User specified XYWorld, but (X,Y)_WORLD was not found. Using (X,Y)_IMAGE for %s.',file.');
    x = dataStruct.X_IMAGE(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
elseif strcmpi(outFormat,'XYImage') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    x = dataStruct.X_IMAGE(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
elseif strcmpi(outFormat,'XYImage') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    warning('User specified XYImage, but (X,Y)_Image was not found. Using (X,Y)_WORLD for %s.',file.');
    x = dataStruct.X_WORLD(~dataStruct.starAndLrgInd);
    y = dataStruct.Y_WORLD(~dataStruct.starAndLrgInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
else
    
    error('Could not find (X_WORLD,Y_WORLD) or (X_IMAGE,Y_IMAGE). Make sure parameter file is properly configured and try again.');
    
end

end


