function y = centerest(y)
%CENTEREST centers the spectral estimates so that DC term is in the center
%   This function is for internal use only. It may be removed in the future. 
%
%   See also CENTERFREQ.
%
%   Copyright 2018 The MathWorks, Inc. 
%#codegen

n = size(y,1);
if iseven(n) 
  %even (nyquist is at end of spectrum)
  y = circshift(y,n/2-1);
else
  % odd
  y = fftshift(y,1);
end
