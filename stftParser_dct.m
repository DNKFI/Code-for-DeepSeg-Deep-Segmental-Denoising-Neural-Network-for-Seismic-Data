function [data,opts] = stftParser_dct(fnName,x,varargin)
%STFTPARSER  is a function for parsing value-only inputs, flags, and
%   name-value pairs for the STFT and ISTFT functions. This function is for
%   internal use only. It may be removed in the future.

%   Copyright 2018-2019 The MathWorks, Inc.
%#codegen

%---------------------------------
% Validate function name
validStrings = {'stft','istft'};
fnName = validatestring(fnName,validStrings,'stftParser','fnName');
isInverse = strcmpi(fnName,'istft');

%---------------------------------
% Set defaults
nwin = 128;
defaultOverlapPercent = 0.75;
defaultWin = hann(nwin,'periodic');
noverlap = floor(nwin.*defaultOverlapPercent);
nfft = nwin;
defaultCentered = true; 
defaultTimeDimension = 'acrosscolumns';
if isInverse % If ISTFT
    % Set ISTFT defaults for ISTFT-specific parameters
    defaultMethod = 'wola';
    defaultConjSym = false;
    
else % Else STFT
    % Place holders in opts structure. Does nothing for STFT.
    method = '';
    conjsym = false;
end

%---------------------------------
% Parse varargin
if ~isempty(varargin)
    idxArgIn = 1:length(varargin);
    parseChk = true(size(varargin));
    
    % Check for value-only input
    if ~any([isa(varargin{1},'char'),isa(varargin{1},'string')])
        timeValue = varargin{1};
        parseChk(1) = false;
    else
        timeValue = [];
    end
        
    % Set indices for n-v parsing; remove indices for value-only inputs and
    % flags
    idxParse = idxArgIn(parseChk);
    
    % Check that the very first index in idxParse is a string or char
    nNameValuePairs = numel(idxParse);
    if nNameValuePairs>=1
        coder.internal.errorIf(~any([isa(varargin{idxParse(1)},'char'),...
            isa(varargin{idxParse(1)},'string')]),sprintf('signal:%s:TooManyValueOnlyInputs',fnName));
    end
    
    % Parse n-v pairs
    if coder.target('MATLAB')
        % MATLAB parser
        pstruct = inputParser;
        pstruct.FunctionName = fnName;
        pstruct.addParameter('Window',defaultWin);
        pstruct.addParameter('OverlapLength',[]);
        pstruct.addParameter('FFTLength',[]);
        pstruct.addParameter('Centered',defaultCentered);
        if isInverse
            pstruct.addParameter('Method',defaultMethod);
            pstruct.addParameter('ConjugateSymmetric',defaultConjSym);
            pstruct.addParameter('InputTimeDimension',defaultTimeDimension);
        else % Else STFT
            pstruct.addParameter('OutputTimeDimension',defaultTimeDimension);
        end
        pstruct.parse(varargin{1,idxParse});
        
        % Get window and validate
        win = pstruct.Results.Window;
        validateattributes(win,{'single','double'},...
            {'nonempty','finite','nonnan','vector'},fnName,'Window');
        win = win(:);
        nwin = length(win);
        validateattributes(nwin,{'numeric'},...
            {'scalar','integer','nonnegative','real','nonnan','nonempty',...
            'finite','>',1},fnName,'WindowLength');
        
        % Assign other n-v pairs dependent on window
        noverlap = ifEmptyUseDefault(pstruct.Results.OverlapLength,...
            floor(nwin.*defaultOverlapPercent));
        nfft = ifEmptyUseDefault(pstruct.Results.FFTLength,nwin);
              
        % Get centered value
        centered = pstruct.Results.Centered;
            
        % Get parser results for ISTFT specific n-v pairs
        if isInverse % If ISTFT
            method = pstruct.Results.Method;
            conjsym = pstruct.Results.ConjugateSymmetric;
            timeDimension = pstruct.Results.InputTimeDimension;
        else % Else STFT
            timeDimension = pstruct.Results.OutputTimeDimension;
        end
    else
        % Codegen parser
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','none', ...
            'StructExpand',false, ...
            'IgnoreNulls',true);
        if isInverse
            parms = {'Window','OverlapLength','FFTLength','Centered','Method','InputTimeDimension'};
        else
            parms = {'Window','OverlapLength','FFTLength','Centered','OutputTimeDimension'};
        end
        if ~isempty(idxParse)
            pstruct = coder.internal.parseParameterInputs(parms,poptions,varargin{1,idxParse});
            
            % Get window and validate
            win = coder.internal.getParameterValue(pstruct.Window,defaultWin,varargin{1,idxParse});
            validateattributes(win,{'single','double'},...
                {'nonempty','finite','nonnan','vector'},fnName,'Window');
            win = win(:);
            nwin = length(win);
            validateattributes(nwin,{'numeric'},...
                {'scalar','integer','nonnegative','real','nonnan','nonempty',...
                'finite','>',1},fnName,'WindowLength');
            
            % Assign other n-v pairs dependent on window
            noverlap = coder.internal.getParameterValue(pstruct.OverlapLength,floor(nwin.*defaultOverlapPercent),varargin{1,idxParse});
            nfft = coder.internal.getParameterValue(pstruct.FFTLength,nwin,varargin{1,idxParse});
            
            % Get centered value
            centered = coder.internal.getParameterValue(pstruct.Centered,defaultCentered,varargin{1,idxParse});
            
            % Get parser results for ISTFT specific n-v pairs
            if isInverse % If ISTFT
                method = coder.internal.getParameterValue(pstruct.Method,defaultMethod,varargin{1,idxParse});
                conjsym = defaultConjSym;
                timeDimension = coder.internal.getParameterValue(pstruct.InputTimeDimension,defaultTimeDimension,varargin{1,idxParse});
            else % Else STFT
                timeDimension = coder.internal.getParameterValue(pstruct.OutputTimeDimension,defaultTimeDimension,varargin{1,idxParse});
            end
        else
            win = defaultWin;
            centered = defaultCentered; 
            timeDimension = defaultTimeDimension;
            if isInverse % If ISTFT
                method = defaultMethod;
                conjsym = defaultConjSym;
            end
        end
    end
else % Use defaults
    win = defaultWin;
    centered = defaultCentered; 
    timeDimension = defaultTimeDimension;
    if isInverse % If ISTFT
        method = defaultMethod;
        conjsym = defaultConjSym;
    end
    timeValue = [];
end

%---------------------------------
% Validate attributes

% OverlapLength
validateattributes(noverlap,{'numeric'},...
    {'scalar','integer','nonnegative','real','nonnan',...
    'finite','nonempty','>=',0,'<',nwin},...
    fnName,'OverlapLength');
noverlap = double(noverlap);

% FFTLength
validateattributes(nfft,{'numeric'},...
    {'scalar','integer','nonnegative','real','nonnan',...
    'finite','nonempty','>=',nwin},fnName,'FFTLength');
nfft = double(nfft);

% Centered
validateattributes(centered,{'numeric','logical'},...
    {'real','nonnan','finite','nonempty','scalar'},fnName,'Centered');

% Validate function specific parameters
if isInverse % If ISTFT
    % Method
    validStrings = {'ola','wola'};
    method = validatestring(method,validStrings,fnName,'Method');
    
    % Conjugate Symmetry
    validateattributes(conjsym,{'numeric','logical'},...
        {'real','nonnan','finite','nonempty','scalar'},fnName,'FFTLength');
    
    % InputTimeDimension - Orientation of T-F map
    validStrings = {'acrosscolumns','downrows'};
    timeDimension = validatestring(timeDimension,validStrings,fnName,'InputTimeDimension');
else % Else STFT
    % OutputTimeDimension - Orientation of T-F map
    validStrings = {'acrosscolumns','downrows'};
    timeDimension = validatestring(timeDimension,validStrings,fnName,'OutputTimeDimension');
end

% Set opts
pVals = struct(...
    'Window',win,...
    'WindowLength',nwin,...
    'OverlapLength',noverlap,...
    'FFTLength',nfft,...
    'Centered',logical(centered),...    
    'TimeDimension',timeDimension,...
    'Method',method,...
    'ConjugateSymmetric',logical(conjsym));

%---------------------------------
% Verify data and time values
[data,opts] = verifyDataAndTime(x,timeValue,pVals,fnName);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function output = ifEmptyUseDefault(input,default)
% Check input. If it is empty, use provided default value.

if isempty(input)
    output = default;
else
    output = input;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [data,opts] = verifyDataAndTime(x,timeValue,pVals,fnName)
% Function to verify data (time signal for STFT or STFT matrix for ISTFT)
% and time values (if provided).
% Set defaults for opts structure
dataLength = 0;
timeAxisLength = 0;
timeMode = 'sp';
timeUnits = '';
initialDate = [];

%---------------------------------
% Validate data
isInverse = strcmpi(fnName,'istft');
if isInverse
    % No timetables for istft
    coder.internal.errorIf(~any(strcmpi(class(x),{'single','double'})),...
        'signal:istft:InvalidInputDataType');
    
    % Validate istft data attributes and assign to data variable
    validateattributes(x,{'single','double'},...
        {'nonsparse','finite','nonnan','3d'},'istft','S');
    
    % Set data
    [isSingle,classCast] = classToCast(pVals.Window,x);
    data = cast(x,classCast);
    if strcmpi(pVals.TimeDimension,'acrosscolumns')
        timeAxisLength = size(data,2);
    else
        timeAxisLength = size(data,1);
    end
    
    numChannels = size(data,3);
else
    if coder.target('MATLAB') && istimetable(x) % For timetables
        
        % ensure numeric input
        if ~all(varfun(@isnumeric,x,'OutputFormat','uniform')) || ~all(varfun(@ismatrix,x,'OutputFormat','uniform'))
            error(message('signal:stft:InvalidTimeTableType'));
        end
        
        % ensure data type is double or single for all the variables
        if ~all(varfun(@(x) isa(x,'double'),x,'OutputFormat','uniform'))...
                && ~all(varfun(@(x) isa(x,'single'),x,'OutputFormat','uniform'))
            error(message('signal:stft:InvalidTimeTableMixedType'));
        end

        % If table has multiple variables, each one can only contain a vector
        if size(x,2) > 1 && ~all(varfun(@(n)size(n,2),x,'OutputFormat','uniform') == 1)
            error(message('signal:stft:InvalidTimeTableType'));
        end
        
        % Set data
        [isSingle,classCast] = classToCast(pVals.Window,x{:,:});
        data = cast(x{:,:},classCast);
        
        % Assign time mode
        timeMode = 'tt';
        
        % Check times
        rowTimes = x.Properties.RowTimes;
        if isduration(rowTimes)
            ttTimeVector = seconds(rowTimes);
            timeUnits = rowTimes.Format;
        else
            d = rowTimes-rowTimes(1);
            ttTimeVector = seconds(d);
            timeUnits = rowTimes.Format;
            initialDate = rowTimes(1);
        end
    else
        % Verify class is of single or double; no timeseries
        coder.internal.errorIf(~any(strcmpi(class(x),{'single','double'})),...
            'signal:stft:InvalidInputDataType');
        
        % Validate data attributes
        validateattributes(x,{'single','double'},...
            {'nonsparse','finite','nonnan','2d','nonempty'},'stft','X');
        
        % Set data
        [isSingle,classCast] = classToCast(pVals.Window,x);
        
        % Treat row vector signal as column vecter signal for consistency
        if isrow(x)
            data = cast(x(:),classCast);
        else
            data = cast(x,classCast);
        end
    end
    [dataLength,numChannels] = size(data);
    coder.internal.errorIf(dataLength<2,'signal:stft:InvalidInputLength');
end

% Check that the number of samples is greater than the window length
if (~isInverse && pVals.WindowLength>size(data,1))
    coder.internal.error(sprintf('signal:%s:InvalidWindowLength',fnName),pVals.WindowLength);
end
if (isInverse && strcmpi(pVals.TimeDimension,'acrosscolumns') && pVals.WindowLength>size(data,1))
    coder.internal.error(sprintf('signal:%s:InvalidWindowLength',fnName),'rows',pVals.WindowLength);
end
if (isInverse && strcmpi(pVals.TimeDimension,'downrows') && pVals.WindowLength>size(data,2))
    coder.internal.error(sprintf('signal:%s:InvalidWindowLength',fnName),'columns',pVals.WindowLength);
end

% Set opts
opts = struct(...
    'DataLength',dataLength,...
    'NumChannels',numChannels,...
    'TimeAxisLength',timeAxisLength,...
    'IsSingle',isSingle,...
    'TimeMode',timeMode,...
    'TimeUnits',timeUnits,...
    'InitialDate',initialDate,...
    'EffectiveFs',2,...
    'IsNormalizedFreq',false,...
    'Window',cast(pVals.Window,classCast),...
    'WindowLength',pVals.WindowLength,...
    'OverlapLength',pVals.OverlapLength,...
    'FFTLength',pVals.FFTLength,...
    'Centered',pVals.Centered,...    
    'Method',pVals.Method,...
    'ConjugateSymmetric',pVals.ConjugateSymmetric,...
    'TimeDimension',pVals.TimeDimension);

% Determine time type
if coder.target('MATLAB') && isdatetime(timeValue)
    coder.internal.error(sprintf('signal:%s:InvalidValueOnlyInput',fnName));
end
if ~isempty(timeValue)
    if coder.target('MATLAB') && isduration(timeValue)
        opts.TimeUnits = timeValue.Format;
        if isscalar(timeValue)
            timeValue = seconds(timeValue);
            opts.TimeMode = 'ts';
        else
            coder.internal.error(sprintf('signal:%s:InvalidValueOnlyInput',fnName));
        end
    else
        if isscalar(timeValue)
            opts.TimeMode = 'fs';
        else
            coder.internal.error(sprintf('signal:%s:InvalidValueOnlyInput',fnName));
        end
    end
end

% Validate time inputs
switch opts.TimeMode
    case 'fs' % Sample rate provided
        if coder.target('MATLAB') && istimetable(x)
            coder.internal.error(sprintf('signal:%s:SampleRateAndTimetableInput',fnName));
        end
        validateattributes(timeValue, {'numeric'},{'scalar','real','finite','positive'},fnName,'sample rate')
        opts.EffectiveFs = double(timeValue);
    case 'ts' % Duration provided
        if coder.target('MATLAB') && istimetable(x)
            coder.internal.error(sprintf('signal:%s:SampleTimeAndTimetableInput',fnName));
        end
        validateattributes(timeValue, {'numeric'},{'scalar','real','finite','positive'},fnName,'sample time')
        opts.EffectiveFs = double(1/timeValue);
    case 'sp' % No time value provided; normalized frequency
        opts.EffectiveFs = 2;
        opts.IsNormalizedFreq = true;
end
% Handle timetable time vector. This is a special case called outside of
% the switch due to codegen.
if coder.target('MATLAB') && strcmpi(opts.TimeMode,'tt') % Timetable
    timeVector = ttTimeVector(:);
    validateattributes(timeVector,{'numeric'},{'vector','nonnan','finite'},fnName,'time values');
    timeVector = double(timeVector);
    opts.EffectiveFs = validateTimeValues(timeVector,fnName);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [isSingle,classCast] = classToCast(win,data)
% Function checks data and window to determine whether window and data
% should be cast to single precision

isSingle = false;
if any([isa(win,'single') isa(data,'single')])
    classCast = 'single';
    isSingle = true;
else
    classCast = 'double';
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function effectiveFs = validateTimeValues(tv,fnName)
% Check regularity of time vector intervals.

% Check values are unique
if length(tv) ~= length(unique(tv))
    coder.internal.error(sprintf('signal:%s:TimeValuesUnique',fnName));
end

% Check values are increasing
if ~issorted(tv)
    coder.internal.error(sprintf('signal:%s:TimeValuesIncreasing',fnName));
end

% Check for irregularity
needsResampling = signal.internal.utilities.isIrregular(tv);

% Calculate effective sampling frequency
if needsResampling
    coder.internal.error(sprintf('signal:%s:TimeValuesNotUniform',fnName));
else
    effectiveFs = 1/mean(diff(tv));
end
