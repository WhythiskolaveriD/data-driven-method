function [ ps] = Phase_calc(fts,ffts)%
%UNTITLED3 Summary of this function goes here
% calculates the phase difference between two time series based on the
% hilbert transform method explained by Phillip et al.

% Phillips, T., R. S. Nerem, B. Fox-Kemper, J. S. Famiglietti, and B. Rajagopalan (2012),
% The influence of ENSO on global terrestrial water storage using GRACE, Geophysical
% Research Letters, 39 (16), L16,705, doi:10.1029/2012GL052495.
%--------------------------------------------------------------------------------
% written by Bramha Dutt Vishwakarma, Institute of Geodesy, University of
% Stuttgart.      ----- 12 July 2015 ----
%--------------------------------------------------------------------------------%%

%% declaration of the variables.
[~,c] = size(fts);

ps(1:c) = 0;
%time = (1:r)';
%rr = isnan(fts);

%time(rr(:,1)==1) = NaN; % time vector
%time = time(all(~isnan(time),2),:); % extract values and leave NaN
fts = fts(all(~isnan(fts),2),:); % extract values and leave NaN
ffts = ffts(all(~isnan(ffts),2),:); % extract values and leave NaN

[rn,~] = size(fts);

 for i=1:c
A = [ones(rn,1) real(hilbert(ffts(:,i))) imag(hilbert(ffts(:,i)))]; % design matrix

abc = (A'*A)\(A'*fts(:,i));

ps(i) = atan2(abc(3,1),abc(2,1))*(180/pi); % phase

 end
end



