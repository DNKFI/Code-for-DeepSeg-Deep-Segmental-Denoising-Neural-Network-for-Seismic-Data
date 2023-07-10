function f = centerfreq(f)
%CENTEREST centers the frequency vector so that DC term is in the center
%   This function is for internal use only. It may be removed in the future. 
%   
%   See also CENTEREST.
%
%   Copyright 2018 The MathWorks, Inc. 
%#codegen

n = numel(f);
if iseven(n)
  %even (nyquist is at end of spectrum)
  f = f - f(n/2);
else
  % odd
  f = f - f((n+1)/2);
end
