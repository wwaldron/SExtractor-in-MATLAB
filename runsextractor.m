function varargout = runsextractor( fitsFiles, configFile, ...
    paramFile, convFile, nnFile, varargin )
%RUNSEXTRACTOR Builds a SExtractor CONFIGFILE and runs SExtractor on the
%FITSFILES
%   Detailed explanation goes here

% Before doing anything, check to see if user is running linux
if ~isunix()
    error('This program must be run on a linux distribution.');
end


% First validate the inputs
prsr = inputParser;
prsr.StructExpand = true;
prsr.addRequired('fitsFiles' ,@(x) assert( (iscellstr(x) && all(cellfun(@exist,x))) || (ischar(x) && logical(exist(x,'file'))) ));
prsr.addRequired('configFile',@(x) assert(ischar(x)));
prsr.addRequired('paramFile' ,@(x) assert(ischar(x) && exist(x,'file')));
prsr.addRequired('convFile'  ,@(x) assert(ischar(x) && exist(x,'file')));
prsr.addRequired('nnFile'    ,@(x) assert(ischar(x) && exist(x,'file')));
prsr.addParameter('SExtractorCommand','sex',@(x) assert( any(strcmpi(x,{'sex','sextractor'})) ));
prsr.addParameter('FitsDualFiles',     {},  @(x) assert( (iscellstr(x) && all(cellfun(@exist,x))) || (ischar(x) && logical(exist(x,'file'))) ));
prsr.addParameter('DetectMinArea',      5,  @(x) assert(isnumeric(x) && isvector(x) && all(x > 0))); % Pixels
prsr.addParameter('DetectThreshold',  1.5,  @(x) assert(isnumeric(x) && isvector(x) && all(x > 0))); % sigma
prsr.addParameter('AnalysisThreshold',1.5,  @(x) assert(isnumeric(x) && isvector(x) && all(x > 0))); % sigma
prsr.addParameter('Gain',               1,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0)));
prsr.addParameter('SeeingFWHM',       0.5,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0))); % ArcSecond
prsr.addParameter('PixelScale',       1.0,  @(x) assert(isnumeric(x) && isvector(x) && all(x >= 0))); % ArcSecond
prsr.addParameter('RemoveCatFiles', false,  @(x) assert(islogical(x)) )
prsr.parse(fitsFiles,configFile,paramFile,convFile,nnFile,varargin{:});

% Make scalars into vectors
fileListLen = length(fitsFiles);
if fileListLen > 1
    wrnSt       = warning;
    warning('off');
    prsr        = struct(prsr);
    warning(wrnSt);
    prsr.Results.DetectMinArea     = makevector(prsr.Results.DetectMinArea,fileListLen);
    prsr.Results.DetectThreshold   = makevector(prsr.Results.DetectThreshold,fileListLen);
    prsr.Results.AnalysisThreshold = makevector(prsr.Results.AnalysisThreshold,fileListLen);
    prsr.Results.Gain              = makevector(prsr.Results.Gain,fileListLen);
    prsr.Results.SeeingFWHM        = makevector(prsr.Results.SeeingFWHM,fileListLen);
    prsr.Results.PixelScale        = makevector(prsr.Results.PixelScale,fileListLen);
end


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
if isempty(prsr.Results.FitsDualFiles)
    dataStruct = runandparsesextractor(fitsFiles,configFile,prsr);
else
    dataStruct = runandparsedual(fitsFiles,configFile,prsr);
end

% If the user wants an output, give the structure containing the Ra and Dec
% values for the detected cosmic rays
if nargout == 1
    varargout{1} = dataStruct;
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


function catStruct = runandparsesextractor(fitsFiles,confFile,prsr)
% This runs SExtractor for as many FITS files were given. This first
% modifies the conf file then runs SExtractor for the image. It outputs the
% catalog to the same directory as the input file.

for i = length(fitsFiles):-1:1
    
    % We want the catalog names and other files to match the input files.
    % Therefore, we shall get the file name and path and use that in the
    % config file
    [fitsPath,fitsName,~] = fileparts(fitsFiles{i});
    catFile = fullfile(fitsPath,[fitsName,'.cat']);
    
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
    
    % Now convert to a struct and delete cat file if necessary
    catStruct(i).inputFile = fitsFiles{i};
    if ~prsr.Results.RemoveCatFiles
        catStruct(i).catFile   = catFile;
    else
        delete(catFile);
        catStruct(i).catFile   = '';
    end
    for j = 1:numHdrLns
        catStruct(i).(catTable.Properties.VariableNames{j}) = catTable.(j);
    end
    
end

end


function catStruct = runandparsedual(fitsFiles,confFile,prsr)
% This sub-routine runs SExtractor in the dual mode where file1 is used to
% detect sources, but the photometry analysis is performed on file2.


for i = length(fitsFiles):-1:1
    
    fitsI = fitsFiles{i};
    
    % Modify the conf file to reflect items that should remain the same per
    % detect file
    confStr = fileread(confFile);
    
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
    [fitsPath,fitsIName,~] = fileparts(fitsFiles{i});
    
    for j = length(prsr.Results.FitsDualFiles):-1:1
        
        fitsJ = prsr.Reults.FitsDualFiles{1};
        
        if strcmpi(fitsI,fitsJ)
            continue;
        end
        
        % Create the CAT File
        [~,fitsJName,~]       = fileparts(fitsJ);
        catFile = fullfile(fitsPath,[fitsIName,'_&_',fitsJName,'.cat']);
        
        % Write the cat file to the conf string
        confStr = replaceline(confStr,'CATALOG_NAME',catFile);
        
        % Write the changes to the file
        fidConf = fopen(confFile,'w');
        fprintf(fidConf,confStr);
        fclose(fidConf);
        
        % Build the SExtractor Command
        cmd = sprintf('%s %s,%s -c %s',prsr.Results.SExtractorCommand,...
            fitsI,fitsJ,confFile);
        
        % Try Running the Command
        trysystem(cmd,confFile);
        
        % Now, get the number of header lines
        numHdrLns = numel(strfind(fileread(catFile),'#'));
        
        % Now store the Variable names
        catFid   = fopen(catFile,'r');
        varNames = cell(1,numHdrLns);
        for k = 1:numHdrLns
            fileLine = fgetl(catFid);
            varNames{k} = strtrim(fileLine(6:29)); % The Variable column is 6:29
        end
        fclose(catFid);
        
        % Get the table and change the names
        catTable = readtable(catFile,'FileType','Text','HeaderLines',numHdrLns);
        catTable.Properties.VariableNames = varNames;
        
        % Now convert to a struct and delete cat file if necessary
        catStruct(i).inputFile = fitsI;
        catStruct(i).dualPhotometry(j).photoFile = fitsJ;
        if ~prsr.Results.RemoveCatFiles
            catStruct(i).dualPhotometry(j).catFile   = catFile;
        else
            delete(catFile);
            catStruct(i).dualPhotometry(j).catFile   = '';
        end
        for k = 1:numHdrLns
            catStruct(i).dualPhotometry(j).(catTable.Properties.VariableNames{k}) = catTable.(k);
        end
        
    end
    
end


end




