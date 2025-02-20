function [xalign, tkl_pts, signal_notalign, cormax] = ...
    alignmentwrt1(signal, startindex, windowlength_samples, ...
    rate, frequencyrange_bins)
%===================================================================
% we align signals wrt to sensor 1
% with structured delays wrt the sensor locations
% That means that from the delays we extract the predicted delays
% with all combinations of sensors
%===================================================================
% [xalign, tkl_pts] = 
%    alignmentwrt1(signal,startindex,windowlength_samples, ...
%      rate, frequencyrange_bins)
%===================================================================
% Inputs:
%   - signal: array size (N,M)
%             N the observation number
%             M sensor number
%   - startindex : beginning index
%   - windowlength_samples :
%   - rate : oersampling rate, tyical 4
%   - frequencyrange_bins: if absent no filtering
%==
% Outputs:
%   - xalign : time-aligned signals
%   - tkl_pts : delays wrt 1 in samples
%
% Used functions
%   - delayedsignalF.m
%====================================================================
% last modified : 8/11
%

signal = signal ./ (ones(size(signal,1),1)*std(signal,[],1));
if nargin == 5
    [forward,  reverse] = butter(2,2*frequencyrange_bins);
    signal              = filter(forward,reverse,signal);
end
signal = signal ./ (ones(size(signal,1),1)*std(signal,[],1));
signal_notalign = signal;

signal    = resample(signal,rate, 1);
windowlength_samples = rate*windowlength_samples;
startindex = (startindex-1)*rate+1;

M         = size(signal,2);
tkl       = zeros(M,1);
cormax    = zeros(M-1,1);

xe        = signal(startindex:startindex+windowlength_samples-1,:);
xt1       = xe(:,1);
for ks              = 1:M-1
    xt2             = xe(:,ks+1);
    corkl           = xcorr(xt1,xt2, 'coeff');
    [bid, indmaxkl] = max(corkl);
    cormax(ks) = bid;
    taukl           = windowlength_samples - indmaxkl;
    % Parabolic interpolation
    yykl            = corkl(indmaxkl-1:indmaxkl+1);
    xxkl            = (-1:1)';
    alphakl         = [ones(3,1) xxkl xxkl.^2]\yykl;
    % we perform index in 1:Combi
    tkl(ks+1)        = (taukl + alphakl(2)/2/alphakl(3));
end
%==
xalign                  = zeros(windowlength_samples,M);
xalign(:,1)             = xe(:,1);
for ss = 2:M
    xalign(:,ss)        = delayedsignalF(xe(:,ss),-tkl(ss));
end
tkl_pts = tkl/rate;
xalign = resample(xalign,1,rate);
%====================================================================
