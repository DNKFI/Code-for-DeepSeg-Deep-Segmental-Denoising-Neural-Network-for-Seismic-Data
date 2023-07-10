function varargout = stft(x,varargin)
%STFT Short-time Fourier transform.
%   S = STFT(X) returns the short-time Fourier transform (STFT) of X. X can
%   be a vector, a matrix or a timetable. If the input has multiple
%   channels, specify X as a matrix where each column corresponds to a
%   channel. If X is a timetable, it must contain finite and uniformly
%   increasing time values. For multichannel timetable input, specify X as
%   a timetable with a single variable containing a matrix or a timetable
%   with multiple variables, each containing a column vector. Precision can
%   be double or single but cannot be mixed. The output S is a matrix for
%   single-channel signals and a 3-D array for multichannel signals. Time
%   increases across the columns and frequency increases down the rows. The
%   third dimension, if present, corresponds to the input channels. If you
%   invert S using ISTFT and want the result to have the same number of
%   time samples NT as X, then (NT - NOVERLAP)/(length(WINDOW) - NOVERLAP)
%   must equal an integer.
%
%   S = STFT(X,Fs) specifies the sample rate of X in hertz as a positive
%   scalar. This parameter provides time information to the input and only
%   applies when X is a vector or a matrix.
%
%   S = STFT(X,Ts) specifies Ts as a positive scalar duration corresponding
%   to the sample time of X. This parameter provides time information to
%   the input and applies only when X is a vector or a matrix. The sample
%   rate in this case is calculated as 1/Ts.
%
%   S = STFT(...,'Window',WINDOW) specifies the window used to compute the
%   STFT. STFT divides each channel of X into segments of the same length
%   as WINDOW with an overlap between adjoining segments set to the largest
%   integer less than or equal to 75% of the window length. The function
%   then windows each segment with the vector specified in WINDOW. Use the
%   function ISCOLA to check a window for constant overlap-add (COLA)
%   compliance. COLA compliance ensures perfect reconstruction for
%   non-modified spectra. The default is a Hann window of length 128.
%
%   S = STFT(...,'OverlapLength',NOVERLAP) specifies an integer number of
%   samples of overlap between adjoining segments. NOVERLAP must be smaller
%   than the length of the window. If NOVERLAP is not specified, it is set
%   to the largest integer less than or equal to 75% of the window length.
%
%   S = STFT(...,'FFTLength',NFFT) specifies the integer number of
%   frequency points used to calculate the discrete Fourier transform. NFFT
%   must be greater than or equal to the window length. NFFT defaults to
%   the length of WINDOW. 
%
%   S = STFT(...,'Centered',CENTERED) returns a two-sided, centered
%   transform if CENTERED is true. The centered transform is computed over
%   (-pi,pi] for even NFFT and over (-pi,pi) for odd NFFT. When time
%   information is provided, the intervals become (-Fs/2,Fs/2] and
%   (-Fs/2,Fs/2), respectively. If CENTERED is false, the function computes
%   S over the interval [0,2*pi) or, if the input signal has time
%   information, over [0,Fs). CENTERED defaults to true.
%
%   S = STFT(...,'OutputTimeDimension',TIMEDIMENSION) specifies the
%   orientation of the STFT according to the location of the time
%   dimension. Set the TIMEDIMENSION value to 'downrows' if you want the
%   time dimension of S down the rows and the frequency dimension across
%   the columns. Set the TIMEDIMENSION value to 'acrosscolumns' if you want
%   the time dimension of S across the columns and the frequency dimension
%   down the rows. This argument is ignored if this function is called with
%   no output arguments. The default value is 'acrosscolumns'.
%
%   [S,F] = STFT(...) returns the frequencies at which the STFT is
%   evaluated. If the input contains time information, then F has units of
%   hertz (Hz). Otherwise it has units of rad/sample.
%
%   [S,F,T] = STFT(...) returns the times at which the STFT is evaluated.
%   If a sample rate is provided, T is a vector that contains time values
%   in seconds. If a sample time is provided, then T is a duration array
%   with the same time format as the input. If no time information is
%   provided, the output is a vector in sample numbers. S has a number of
%   rows equal to the length of the frequency vector F and a number of
%   columns equal to the length of the time vector T. 
%
%   STFT(...) with no output arguments plots the two-sided magnitude of the
%   STFT. This syntax does not support multichannel signals.
%
%    % EXAMPLE 1: 
%       % Compute and display the STFT of a chirp with normalized
%       % frequency.
%       fs = 4096;
%       t = 0:1/fs:2-1/fs;
%       x = chirp(t,250,1,500,'q');
%       stft(x)
%       
%    % EXAMPLE 2: 
%       % Compute and display the STFT of a quadratic chirp with a duration
%       % of 1 ms.
%       t = -2:1/1e3:2;
%       fo = 100;
%       f1 = 200;
%       x = chirp(t,fo,1,f1,'quadratic',[],'concave');
%       D = seconds(1e-3);
%       win = hamming(100,'periodic');
%       stft(x,D,'Window',win,'OverlapLength',98,'FFTLength',128)
% 
%    % EXAMPLE 3: 
%       % Compute and display the STFT of a complex signal.
%       fs = 3000;
%       t = 0:1/fs:1-1/fs;
%       x = exp(2j*pi*100*cos(2*pi*2*t))+randn(size(t))/100;
%       win = hann(100,'periodic');
%       stft(x,fs,'Window',win,'OverlapLength',50,'FFTLength',200)
% 
%    % EXAMPLE 4: 
%       % Compute and display the STFT of a set of intermittent sinusoid
%       % signals.
%       silence = zeros(1,1500);
%       fs = 1e3;
%       t = (0:1000-1)/fs;
%       yStep = [sin(2*pi*50*t) silence sin(2*pi*100*t) silence sin(2*pi*150*t)].';
%       t = seconds((0:length(yStep)-1)/fs).';
%       xTable = timetable(t,yStep);
%       win = blackman(300,'periodic');
%       [S,F,T] = stft(xTable,'Window',win,'OverlapLength',250);
%       smag = mag2db(abs(S));
%       pcolor(seconds(T),F,smag)
%       xlabel('Time (s)')
%       ylabel('Frequency (Hz)')
%       shading flat
%       colorbar
%       caxis(max(smag(:)) + [-60 0]) 
%
%    % EXAMPLE 5: 
%       % Compute and display the STFT of a chirp with normalized
%       % frequency over the interval 0 to 2*pi. 
%       fs = 4096;
%       t = 0:1/fs:2-1/fs;
%       x = chirp(t,250,1,500,'q');
%       stft(x,'Centered',false) 
%
%    % EXAMPLE 6:
%       % Compute the STFT of a three-channel signal
%       fs = 4096;
%       t = 0:1/fs:2-1/fs;
%       x = [chirp(t,250,1,500,'quadratic',[],'concave');
%           chirp(t,250,1,500,'quadratic',[],'convex');
%           chirp(t,250,1,500)]'; 
%       [S,F,T] = stft(x,fs);
%       % Plot the time-frequency map for the first channel
%       smag = mag2db(abs(S(:,:,1)));
%       pcolor(T,F, smag)
%       xlabel('Time (s)')
%       ylabel('Frequency (Hz)')
%       shading flat
%       colorbar
%       caxis(max(smag(:)) + [-60 0])
%       title('Input channel: 1')
%
%   See also ISTFT, ISCOLA, PSPECTRUM

% [1] Mitra, S. K. Digital Signal Processing: A Computer-Based Approach.
%     2nd Ed. New York: McGraw-Hill, 2001.

%   Copyright 2018-2019 The MathWorks, Inc.
%#codegen

%---------------------------------
% Check inputs/outputs
narginchk(1,12);
if coder.target('MATLAB') % For MATLAB
    nargoutchk(0,3);
else 
    nargoutchk(1,3);
end

%---------------------------------
% Parse inputs
[data,opts] = stftParser_dct('stft',x,varargin{:});

% No convenience plot for multichannel signals
coder.internal.errorIf(nargout == 0 && (opts.NumChannels >1),...
    'signal:stft:InvalidNumOutputMultiChannel');


%---------------------------------
% Compute STFT
[S,F,T] = computeSTFT(data,opts); 

%---------------------------------
% Set outputs 

% Format time output
if coder.target('MATLAB') && ~isempty(opts.InitialDate)
    % Set times to datetime format if time information is in datetimes
    T = seconds(T)+opts.InitialDate;
end

% Convenience plot
if nargout==0
    displaySTFT(T,F,S,opts);
end

% Set varargout
if nargout >= 1
    if strcmpi(opts.TimeDimension,'acrosscolumns')
        varargout{1} = S;
    else
        % Non-conjugate transpose if output T-F map does not have the
        % default orientation
        varargout{1} = permute(S,[2,1,3]);
    end
end
if nargout >= 2
    if opts.IsNormalizedFreq
        varargout{2} = F.*pi; % rad/sample
    else
        varargout{2} = F;
    end
end
if nargout == 3
   if coder.target('MATLAB') && isnumeric(T) && ~isempty(opts.TimeUnits)
       T = duration(0,0,T,'Format',opts.TimeUnits);
   end
   % Ensure time output is a column vector
   varargout{3} = T(:);
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [S,F,T] = computeSTFT(x,opts)
% Computes the short-time Fourier transform 
classCast = class(x); 

% Set variables 
nx = opts.DataLength;
win = opts.Window;
nwin = opts.WindowLength;
noverlap = opts.OverlapLength;
nfft = opts.FFTLength;
Fs = opts.EffectiveFs;

% Place x into columns and return the corresponding central time estimates
[xin,t] = getSTFTColumns(x,nx,nwin,noverlap,Fs);

% Apply the window to the array of offset signal segments and perform a DFT
[S,f] = computeDFT_dct(bsxfun(@times,win,xin),nfft,Fs);

% Center outputs
if opts.Centered
    [S,f] = centerOutputs(S,f);
end

% Scale frequency and time vectors in the case of normalized frequency
if opts.IsNormalizedFreq
    t = t.*opts.EffectiveFs; % samples 
end

% Set outputs
F = cast(f,classCast); 
T = cast(t,classCast); 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [S,F] = centerOutputs(S,F)
% Centers S and F

S = centerest(S);
F = centerfreq(F);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function displaySTFT(T,F,S,opts)
% Plots STFT

% Plot options
plotOpts.isFsnormalized = opts.IsNormalizedFreq;

% The function plotTFR has its own scaling. Remove scaling as necessary to
% provide plotTFR with expected inputs.
if plotOpts.isFsnormalized
    if ~opts.Centered
        timeScale = 1/(2*pi);
        freqScale = pi;
        T = T.*timeScale;
        F = F.*freqScale; 
    end
end
plotOpts.cblbl = getString(message('signal:dspdata:dspdata:MagnitudedB'));
plotOpts.title = 'Short-time Fourier Transform'; 
plotOpts.threshold = max(20*log10(abs(S(:))+eps))-60;

% Magnitude (dB) 
signalwavelet.internal.convenienceplot.plotTFR(T,F,20*log10(abs(S)+eps),plotOpts); 
