function varargout = finduncommonsources( fitsFiles, configFile, ...
    paramFile, convFile, nnFile, varargin )
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


% Before doing anything, check to see if user is running linux
if ~isunix()
    error('This program must be run on a linux distribution.');
end


% First validate the inputs
prsr = inputParser;
prsr.StructExpand = true;
prsr.addRequired('fitsFiles' ,@(x) assert(iscellstr(x) && all(cellfun(@exist,x) )));
prsr.addRequired('configFile',@(x) assert(ischar(x)));
prsr.addRequired('paramFile' ,@(x) assert(ischar(x) && exist(x,'file')));
prsr.addRequired('convFile'  ,@(x) assert(ischar(x) && exist(x,'file')));
prsr.addRequired('nnFile'    ,@(x) assert(ischar(x) && exist(x,'file')));
prsr.addParameter('SExtractorCommand','sex',@(x) assert(any(strcmp(x,{'sex','sextractor'}))));
prsr.addParameter('DetectMinArea',      5,  @(x) assert(isnumeric(x) && isvector(x) && all(x > 0))); % Pixels
prsr.addParameter('DetectThreshold',  1.5,  @(x) assert(isnumeric(x) && isvector(x) && all(x > 0))); % sigma
prsr.addParameter('AnalysisThreshold',1.5,  @(x) assert(isnumeric(x) && isvector(x) && all(x > 0))); % sigma
prsr.addParameter('MaxAngSepTrueSrc', 0.2,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0))); % ArcSecond
prsr.addParameter('Gain',               1,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0)));
prsr.addParameter('SeeingFWHM',       0.5,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0))); % ArcSecond
prsr.addParameter('PixelScale',       1.0,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0))); % ArcSecond
prsr.addParameter('RegOutputFormat','XYWorld',@(x) assert(any(strcmpi(x,{'XYWorld','XYImage','saoimage','saoworld'}))));
prsr.parse(fitsFiles,configFile,paramFile,convFile,nnFile,varargin{:});

% Make scalars into vectors
fileListLen = length(fitsFiles);
wrnSt       = warning;
warning('off');
prsr        = struct(prsr);
warning(wrnSt);
prsr.Results.DetectMinArea     = makevector(prsr.Results.DetectMinArea,fileListLen);
prsr.Results.DetectThreshold   = makevector(prsr.Results.DetectThreshold,fileListLen);
prsr.Results.AnalysisThreshold = makevector(prsr.Results.AnalysisThreshold,fileListLen);
prsr.Results.MaxAngSepTrueSrc  = makevector(prsr.Results.MaxAngSepTrueSrc,fileListLen);
prsr.Results.Gain              = makevector(prsr.Results.Gain,fileListLen);
prsr.Results.SeeingFWHM        = makevector(prsr.Results.SeeingFWHM,fileListLen);
prsr.Results.PixelScale        = makevector(prsr.Results.PixelScale,fileListLen);


% Create the config file if it does not exist
if ~logical(exist(configFile,'file'))
    % Build the command
    cmd = sprintf('%s -d > %s', prsr.Results.SExtractorCommand,configFile);
    
    % Run the command
    trysystem(cmd,configFile);
end

% Read in the text for the conf file
confStr = fileread(configFile);

% Now, change values for the optional parameters. If no inputs are given,
% the file takes precedent. That is, if the user provides no optional
% inputs for any of the above values, the values in the file will not be
% changed.
confStr = replaceline(confStr,'PARAMETERS_NAME',paramFile);
confStr = replaceline(confStr,'FILTER_NAME',convFile);
confStr = replaceline(confStr,'CHECKIMAGE_TYPE','NONE');
confStr = replaceline(confStr,'STARNNW_NAME',nnFile);

% Write the changes to the file
fidConf = fopen(configFile,'w');
fprintf(fidConf,confStr);
fclose(fidConf);

% Run SExtractor
dataStruct = runsextractor(fitsFiles,configFile,prsr);

% Correlate the Data to Find Stars
dataStruct = outputfalsesources(dataStruct,prsr.Results.MaxAngSepTrueSrc,...
    prsr.Results.RegOutputFormat);


% If the user wants an output, give the structure containing the Ra and Dec
% values for the detected cosmic rays
if nargout == 1
    varargout{1} = dataStruct;
end

end


function vec = makevector(inpt,makeLen)

if isscalar(inpt)
    vec = repmat(inpt,[1,makeLen]);
elseif length(inpt) ~= makeLen
    warning('A parameter input did not match the length of the number of input files.\nOnly using first element.');
    vec = repmat(inpt(1),[1,makeLen]);
else
    vec = inpt;
end

end


function trysystem(cmd,configFile)

% Try the command, but catch the error if one occurs. An error could be
% because the user does not have SExtractor installed, or the user does
% not have write priveliges to the directory pointed to by configFile.
[stat,msg] = system(cmd);
if logical(stat) && any(strfind(msg,'Permission'))
    fprintf(2,'SExtractor Error Message: %s\n',msg);
    error('An error occurred. It seems you don''t have write priveleges to %s\n',configFile);
elseif logical(stat) && any(strfind(msg,'directory'))
    fprintf(2,'SExtractor Error Message: %s\n',msg);
    errMsg = sprintf('An error occurred. It seemes the directory does not exist.\n');
    error(errMsg);
elseif logical(stat) && any(strfind(cmd, '-d'))
    delete(configFile)
    fprintf(2,'SExtractor Error Message: %s\n',msg);
    errMsg = sprintf('An error occurred. Perhaps SExtractor is not installed.\n\tRecommend running ''sudo apt-get install sextractor'' at the command line.\n\tNote that the newest version uses the command ''sextractor'' not ''sex''.\n');
    error(errMsg);
elseif logical(stat)
    fprintf(2,'SExtractor Error Message: %s\n',msg);
    errMsg = sprintf('An error occurred. Perhaps SExtractor is not installed.\n\tRecommend running ''sudo apt-get install sextractor'' at the command line.\n\tNote that the newest version uses the command ''sextractor'' not ''sex''.\n');
    error(errMsg); %#ok<*SPERR>
end

end


function fileStr = replaceline(fileStr,lineIdnt,newVal,varargin)

% If there are arguments in varargin, it means we should check to see if
% the parameter should be replaced. If the prsr object is using a default
% value (that is, one not passed to the parent function), then we do not
% want to replace the value in our string. However, if the user provided
% input, then we shall change the value.
if ~isempty(varargin)
    
    curOpt    = varargin{1};
    defaults  = varargin{2};
    
    % If it is not using the default value, then change the line
    if ~any(strcmpi(curOpt,defaults))
        % The item to change can be a number or string
        if isnumeric(newVal)
            repStr  = sprintf('%-17s%f',lineIdnt,newVal);
        else
            repStr  = sprintf('%-17s%s',lineIdnt,newVal);
        end
        lineIdnt = ['^.*',lineIdnt,'\s.*$'];
        fileStr = regexprep(fileStr,lineIdnt,repStr,...
            'lineanchors','dotexceptnewline','once');
    end
    
else % The string should always be replaced
    
    repStr  = sprintf('%-17s%s',lineIdnt,newVal);
    lineIdnt = ['^.*',lineIdnt,'\s.*$'];
    fileStr = regexprep(fileStr,lineIdnt,repStr,...
        'lineanchors','dotexceptnewline','once');
    
end

end


function catStruct = runsextractor(fitsFiles,confFile,prsr)
% This runs SExtractor for as many FITS files were given. This first
% modifies the conf file then runs SExtractor for the image. It outputs the
% catalog to the same directory as the input file.

for i = length(fitsFiles):-1:1
    
    % We want the catalog names and other files to match the input files.
    % Therefore, we shall get the file name and path and use that in the
    % config file
    if any(strfind(fitsFiles{i},'.fits')) % If the extension is lowercase
        catFile = strrep(fitsFiles{i},'.fits','.cat');
    elseif any(strfind(fitsFiles{i},'.FITS')) % If the extension is uppercase
        catFile = strrep(fitsFiles{i},'.FITS','.cat');
    else % Just try sticking the extension on the end
        catFile = [fitsFiles{i},'.cat'];
    end
    
    % Load in the config file
    confStr = fileread(confFile);
    
    % Modify the config file
    confStr = replaceline(confStr,'CATALOG_NAME',catFile);
    
    % Replace certain conf params
    confStr = replaceline(confStr,'DETECT_MINAREA',prsr.Results.DetectMinArea(i),...
        'DetectMinArea',prsr.UsingDefaults);
    confStr = replaceline(confStr,'DETECT_THRESH',prsr.Results.DetectThreshold(i),...
        'DetectThreshold',prsr.UsingDefaults);
    confStr = replaceline(confStr,'ANALYSIS_THRESH',prsr.Results.AnalysisThreshold(i),...
        'AnalysisThreshold',prsr.UsingDefaults);
    confStr = replaceline(confStr,'GAIN',prsr.Results.Gain(i),...
        'Gain',prsr.UsingDefaults);
    confStr = replaceline(confStr,'SEEING_FWHM',prsr.Results.SeeingFWHM(i),...
        'SeeingFWHM',prsr.UsingDefaults);
    confStr = replaceline(confStr,'PIXEL_SCALE',prsr.Results.PixelScale(i),...
        'PixelScale',prsr.UsingDefaults);
    
    % Write the changes to the file
    fidConf = fopen(confFile,'w');
    fprintf(fidConf,confStr);
    fclose(fidConf);
    
    % Build the SExtractor Command
    cmd = sprintf('%s %s -c %s',prsr.Results.SExtractorCommand,fitsFiles{i},confFile);
    
    % Try Running the Command
    trysystem(cmd,confFile);
    
    % Now, get the number of header lines
    numHdrLns = numel(strfind(fileread(catFile),'#'));
    
    % Now store the Variable names
    catFid   = fopen(catFile,'r');
    varNames = cell(1,numHdrLns);
    for j = 1:numHdrLns
        fileLine = fgetl(catFid);
        varNames{j} = strtrim(fileLine(6:29)); % The Variable column is 6:29
    end
    fclose(catFid);
    
    % Get the table and change the names
    catTable = readtable(catFile,'FileType','Text','HeaderLines',numHdrLns);
    catTable.Properties.VariableNames = varNames;
    
    % Now convert to a struct
    catStruct(i).inputFile = fitsFiles{i};
    catStruct(i).catFile   = catFile;
    for j = 1:numHdrLns
        catStruct(i).(catTable.Properties.VariableNames{j}) = catTable.(j);
    end
    
end

end


function dataStruct = outputfalsesources(dataStruct,maxSep,outFormat)

nFiles = length(dataStruct);

for i = 1:nFiles
    dataStruct(i).starInd = false(size(dataStruct(i).X_WORLD));
end

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
    
    if maxSep(i) > 0
        
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
            curStarInd = (dTh <= maxSep(i));
            
            % Update the logical arrays for each waveband indicating where
            % stars are
            starI = any(curStarInd,1);
            starJ = any(curStarInd,2);
            dataStruct(i).starInd = dataStruct(i).starInd | starI(:);
            dataStruct(j).starInd = dataStruct(j).starInd | starJ(:);
            
        end
        
        % Write those false sources to a file
        [catPath,catFile,~] = fileparts(dataStruct(i).catFile);
        falseRegFile = fullfile(catPath,[catFile,'_falsesource','.reg']);
        writetoregfile(falseRegFile,outFormat,dataStruct(i));
        
    else
        
        % Write those false sources to a file
        [catPath,catFile,~] = fileparts(dataStruct(i).catFile);
        allRegFile = fullfile(catPath,[catFile,'_allsource','.reg']);
        writetoregfile(allRegFile,outFormat,dataStruct(i));
        
    end
    
end

end


function writetoregfile(file,outFormat,dataStruct)


% If desired format is SAOImage and the right fields are in place.
if strcmpi(outFormat,'saoimage') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE'}))
    
    tmpData = [dataStruct.X_IMAGE,dataStruct.Y_IMAGE,dataStruct.A_IMAGE,...
        dataStruct.B_IMAGE,dataStruct.THETA_IMAGE];
    fid = fopen(file,'w');
    fprintf(fid,'# format: image\n');
    fprintf(fid,'ellipse(%.10f,%.10f,%.10f,%.10f,%.10f)\n',tmpData);
    fclose(fid);
    
% If desired format is SAOImage but don't have ellipse parameters, but do
% have centers in image space
elseif strcmpi(outFormat,'saoimage') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    warning('User chose SAOImage, but ellipse parameters not provided. Writing out (X,Y)_IMAGE data for %s.',file);
    x = dataStruct.X_IMAGE(~dataStruct.starInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
% If desired format is SAOImage but don't have ellipse parameters but do
% have ceters in world space
elseif strcmpi(outFormat,'saoimage') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    warning('User chose SAOImage, but ellipse parameters not provided. Writing out (X,Y)_WORLD data for %s.',file.');
    x = dataStruct.X_WORLD(~dataStruct.starInd);
    y = dataStruct.Y_WORLD(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
% If desired format is SAOWorld and have correct fields
elseif strcmpi(outFormat,'saoworld') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD','A_WORLD','B_WORLD','THETA_WORLD'}))
    
    tmpData = [dataStruct.X_WORLD,dataStruct.Y_WORLD,dataStruct.A_WORLD,...
        dataStruct.B_WORLD,dataStruct.THETA_WORLD];
    fid = fopen(file,'w');
    fprintf(fid,'# format: degrees (fk5)\n');
    fprintf(fid,'+ellipse(%.10f,%.10f,%.10f,%.10f,%.10f) #green\n',tmpData);
    fclose(fid);
    
% If desired format is SAOWorld but missing ellipse params but have centers
% in world coordinates
elseif strcmpi(outFormat,'saoworld') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    warning('User chose SAOWorld, but ellipse parameters not provided. Writing out (X,Y)_WORLD data for %s.',file.');
    x = dataStruct.X_WORLD(~dataStruct.starInd);
    y = dataStruct.Y_WORLD(~dataStruct.starInd);
    dlmwrite(file,[x,y]);
    
% If desired format is SAOWorld but missing ellipse params but have centers
% in image coordinates
elseif strcmpi(outFormat,'saoworld') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    warning('User chose SAOWorld, but ellipse parameters not provided. Writing out (X,Y)_IMAGE data for %s.',file.');
    x = dataStruct.X_IMAGE(~dataStruct.starInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
% You get the picture . . .
elseif strcmpi(outFormat,'XYWorld') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    x = dataStruct.X_WORLD(~dataStruct.starInd);
    y = dataStruct.Y_WORLD(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
elseif strcmpi(outFormat,'XYWorld') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    warning('User specified XYWorld, but (X,Y)_WORLD was not found. Using (X,Y)_IMAGE for %s.',file.');
    x = dataStruct.X_IMAGE(~dataStruct.starInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
elseif strcmpi(outFormat,'XYImage') && ...
        all(isfield(dataStruct,{'X_IMAGE','Y_IMAGE'}))
    
    x = dataStruct.X_IMAGE(~dataStruct.starInd);
    y = dataStruct.Y_IMAGE(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
elseif strcmpi(outFormat,'XYImage') && ...
        all(isfield(dataStruct,{'X_WORLD','Y_WORLD'}))
    
    warning('User specified XYImage, but (X,Y)_Image was not found. Using (X,Y)_WORLD for %s.',file.');
    x = dataStruct.X_WORLD(~dataStruct.starInd);
    y = dataStruct.Y_WORLD(~dataStruct.starInd);
    dlmwrite(file,[x,y],'delimiter','\t',...
            'precision','%.10f');
    
else
    
    error('Could not find (X_WORLD,Y_WORLD) or (X_IMAGE,Y_IMAGE). Make sure parameter file is properly configured and try again.');
    
end

end


