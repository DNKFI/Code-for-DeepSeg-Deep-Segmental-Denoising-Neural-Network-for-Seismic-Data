function [X,T] = istft(S,varargin)
%ISTFT Inverse short-time Fourier transform.
%   X = ISTFT(S) returns the inverse short-time Fourier transform (ISTFT)
%   of S. For single-channel signals, specify S as a matrix with time
%   increasing across the columns and frequency increasing down the rows.
%   For multichannel signals, specify S as a 3-D array with the third
%   dimension corresponding to the channels. S is expected to be two-sided
%   and centered.
%
%   X = ISTFT(S,Fs) specifies the sample rate of X in hertz as a positive
%   scalar.
%
%   X = ISTFT(S,Ts) specifies Ts as a positive scalar duration
%   corresponding to the sample time of X. The sample rate in this case is
%   calculated as 1/Ts.
%
%   X = ISTFT(...,'Window',WINDOW) specifies the window used in calculating
%   the ISTFT. Perfect time-domain reconstruction requires the ISTFT window
%   to match the window used to generate the STFT. Use the function ISCOLA
%   to check a window/overlap combination for constant overlap-add (COLA)
%   compliance. COLA compliance is a requirement for perfect reconstruction
%   for non-modified spectra. The default is a Hann window of length 128.
%
%   X = ISTFT(...,'OverlapLength',NOVERLAP) specifies an integer number of
%   samples of overlap between adjoining segments. NOVERLAP must be smaller
%   than the length of the window. Perfect time-domain reconstruction
%   requires the ISTFT NOVERLAP to match the NOVERLAP used to generate the
%   STFT. Use the function ISCOLA to check a window/overlap combination for
%   constant overlap-add (COLA) compliance. COLA compliance is a
%   requirement for perfect reconstruction for non-modified spectra. If
%   NOVERLAP is not specified, it is set to the largest integer less than
%   or equal to 75% of the window length.
%
%   X = ISTFT(...,'FFTLength',NFFT) specifies the integer number of
%   frequency points used to calculate the discrete Fourier transform.
%   Perfect time-domain reconstruction requires the ISTFT NFFT to match the
%   NFFT used to generate the STFT. NFFT defaults to the length of WINDOW.
%
%   X = ISTFT(...,'Method',METHOD) specifies the inversion method to use:
%       'ola'  - Overlap-Add 
%       'wola' - Weighted Overlap-Add
%   If the method is set to 'wola', a second window is applied after the
%   inverse DFT and prior to the final overlap-add stage, and a Griffin-Lim
%   normalization is performed. Typically, the analysis window used to
%   generate the STFT is the same as the synthesis window used during the
%   ISTFT. METHOD defaults to 'wola'.
% 
%   X = ISTFT(...,'ConjugateSymmetric',CONJUGATESYMMETRIC) is specified as
%   a logical true if S is symmetric or a logical false if S is
%   nonsymmetric. When S (the STFT matrix input to the ISTFT function) is
%   not exactly conjugate symmetric due to round-off error,
%   CONJUGATESYMMETRIC set to true ensures S is treated as if it were
%   conjugate symmetric. If S is conjugate symmetric, then the inverse
%   transform computation is faster, and the output is real. This
%   name-value pair is not supported for code generation. Defaults to
%   false.
%
%   X = ISTFT(...,'Centered',CENTERED) treats the input S as a two-sided,
%   centered transform if CENTERED is true. The function rearranges S so
%   that the zero-frequency component is the first row of the array. If
%   CENTERED is false, S is not rearranged. CENTERED defaults to true.
%
%   X = ISTFT(...,'InputTimeDimension',TIMEDIMENSION) specifies the
%   orientation of input S according to the location of the time dimension.
%   If TIMEDIMENSION is set to 'downrows', ISTFT assumes that the time
%   dimension of S is down the rows and the frequency dimension is across
%   the columns. If TIMEDIMENSION is set to 'acrosscolumns', ISTFT assumes
%   that the time dimension of S is across the columns and the frequency
%   dimension is down the rows. The default value is 'acrosscolumns'.
%
%   [X,T] = ISTFT(...) returns time vector T. If a sample rate is provided,
%   T is a vector of time values in seconds. If a sample time is provided,
%   then T is a duration array with the same time format as the input. If
%   no time information is provided, the output is a vector of sample
%   numbers. X contains the reconstructed time-domain signal for all the
%   input channels.
%
%    % EXAMPLE 1: 
%       % Compute the ISTFT of a real signal using the overlap-add method.
%       fs = 10240;
%       t = 0:1/fs:0.5-1/fs;
%       x = 5*sin(2*pi*t*10);
%       D = duration(0,0,1/fs);
%       win = hamming(512,'periodic');
%       S = stft(x,D,'Window',win,'OverlapLength',384,'FFTLength',1024);
%       [X,T] = istft(S,D,'Window',win,'OverlapLength',384,...
%           'FFTLength',1024,'Method','ola','ConjugateSymmetric',true); 
%       
%       % Plot original and resynthesized signals.  
%       plot(t,x,seconds(T),X,'-.')
%       axis tight
%       xlabel('Time (s)')
%       ylabel('Amplitude (V)')
%       title('Original and Reconstructed Signal')
%       legend('Original','Reconstructed')
% 
%    % EXAMPLE 2: 
%       % Compute the ISTFT of a complex signal using weighted overlap-add
%       % method.
%       fs = 3000;
%       t = 0:1/fs:1-1/fs;
%       x = exp(2j*pi*100*cos(2*pi*2*t))+randn(size(t))/100;
%       nwin = 100; 
%       win = hann(nwin,'periodic');
%       xZero = [zeros(1,nwin) x  zeros(1,nwin)]; % Zero pad to fix edges 
%       S = stft(xZero,fs,'Window',win,'OverlapLength',75);
%       [X,T] = istft(S,fs,'Window',win,'OverlapLength',75);
%       
%       % Remove zeros.
%       X(1:nwin) = []; 
%       X(end-nwin+1:end) = []; 
%       T = T(1:end-2*nwin); 
%       
%       % Plot original and resynthesized signals. 
%       plot(t,abs(x),T,abs(X),'-.')
%       xlabel('Time (s)')
%       ylabel('Amplitude (V)')
%       title('Original and Reconstructed Signal')
%       legend('Original','Reconstructed')
%
%    % EXAMPLE 3:
%       % Compute the ISTFT of a multichannel signal
%       fs = 4096;
%       t = 0:1/fs:2-1/fs;
%       x = [chirp(t,250,1,500,'quadratic',[],'concave');
%           chirp(t,250,1,500,'quadratic',[],'convex');
%           chirp(t,250,1,500)]';
%       win = hamming(256,'periodic');   
%       S = stft(x,fs,'Window',win,'OverlapLength',128,'FFTLength',1024);
%       [X,T] = istft(S,fs,'Window',win,'OverlapLength',128,...
%           'FFTLength',1024,'Method','ola');
% 
%       %Plot the first 0.1 seconds of the first channel of the original signal. Overlay the first channel of the reconstructed signal.
%       plot(t,x(:,1),T,X(:,1),'-.')
%       xlabel('Time (s)')
%       ylabel('Amplitude (V)')
%       title('Original and Reconstructed Signal')
%       legend('Original','Reconstructed')
%       xlim([0, 0.1])
%       title('Input channel: 1')
% 
%   See also STFT, ISCOLA, PSPECTRUM, IFFT

% [1] Crochiere, R. E. "A Weighted Overlap-Add Method of Short-Time
%     Fourier Analysis/Synthesis." IEEE Transactions on Acoustics, Speech,
%     and Signal Processing. Vol. ASSP-28, Feb. 1980, pp. 99-102.
% [2] Griffin, D. W. and J. S. Lim. "Signal Estimation from Modified
%     Short-Time Fourier Transform." IEEE Transactions on Acoustics,
%     Speech, and Signal Processing. Vol. ASSP-32, No. 2, April 1984.
% [3] Portnoff, M. R. "Time-Frequency Representation of Digital Signals
%     and Systems Based on Short-Time Fourier Analysis." IEEE Transactions
%     on Acoustics, Speech, and Signal Processing. Vol. ASSP-28, Feb. 1980,
%     pp. 55-69.
        
%   Copyright 2018-2019 The MathWorks, Inc.
%#codegen

%---------------------------------
% Check inputs/outputs
narginchk(1,16);
if coder.target('MATLAB') % For MATLAB
    nargoutchk(0,2);
else % For everything else 
    nargoutchk(1,2);
end

%---------------------------------
% Parse inputs
[data,opts] = stftParser_dct('istft',S,varargin{:});

%---------------------------------
% Non-conjugate transpose if T-F map input does not have the default
% orientation
if strcmpi(opts.TimeDimension,'acrosscolumns')
    inputData = data;
else
    inputData = permute(data,[2,1,3]);
end

%---------------------------------
% Compute ISTFT
[X,T] = computeISTFT(inputData,opts);

%---------------------------------
% Update time vector if a duration was provided
if coder.target('MATLAB') && strcmpi(opts.TimeMode,'ts')
    T = seconds(T);
    T.Format = opts.TimeUnits;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [X,T] = computeISTFT(s,opts)
% Computes inverse short-time Fourier transform 
classCast = class(s); 

% Set variables 
win = opts.Window;
numCh = opts.NumChannels;
nwin = opts.WindowLength;
noverlap = opts.OverlapLength;
nfft = opts.FFTLength;
hop = nwin-noverlap;
nseg = opts.TimeAxisLength;
xlen = nwin + (nseg-1)*hop;
Fs = opts.EffectiveFs;

% Uncenter
if opts.Centered
    s = uncenter(s);
end

% IDFT
if coder.target('MATLAB') && opts.ConjugateSymmetric
    xifft = dct(s,nfft,'Type',4);
else
    xifft = dct(s,nfft,'Type',4);
end

xifft = xifft(1:min(nwin,size(xifft,1)),1:nseg,1:numCh);

% Initialize time-domain signal
if isreal(xifft)
    x = zeros(xlen,1,numCh,classCast);
else
    x = complex(zeros(xlen,1,numCh,classCast));
end

% Set method
if strcmpi(opts.Method,'ola') 
    a = 0;
else % Else WOLA
    a = 1;
end

% Initialize normalization value
normVal = zeros(xlen,numCh);
winNominator = repmat(win.^a,1,1,numCh);
winDeNominator = repmat(win.^(a+1),1,numCh);

% Overlap-add
for ii = 1:nseg
    x(((ii-1)*hop+1):((ii-1)*hop+nwin),1,:) = x(((ii-1)*hop+1):((ii-1)*hop+nwin),1,:) ...
        + xifft(:,ii,:).*winNominator;
    normVal(((ii-1)*hop+1):((ii-1)*hop+nwin),:) = normVal(((ii-1)*hop+1):((ii-1)*hop+nwin),:)+winDeNominator;
end

% Normalize
normVal(normVal<(nseg*eps)) = 1; % Don't normalize by small values
X = squeeze(x)./normVal;

% Time vector
T = cast(((0:size(x,1)-1).')./Fs,classCast);

% Scale time vector in the case of normalized frequency
if opts.IsNormalizedFreq
    T = T.*opts.EffectiveFs; % sample
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function s = uncenter(s)
% Uncenter s

n = size(s,1);
if iseven(n)
  % even (nyquist is at end of spectrum)
  s = circshift(s,-(n/2-1));
else
  % odd
  s = ifftshift(s,1);
end