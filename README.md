# Neural-Data-Analysis-Crap
function y = filterSignal(x, Fs)
% Filter raw signal
%   y = filterSignal(x, Fs) filters the signal x. Each column in x is one
%   recording channel. Fs is the sampling frequency. The filter delay is
%   compensated in the output y.
Lpass=300;
Hstop=3000;
fBand=[Lpass/(Fs/2),Hstop/(Fs/2)];
[b,a]=cheby2(6,50,fBand); % create the chebyshev type II bandpass filter
% h=fvtool(b,a) visualize the filter if needed
for i=1:4
    y(:,i)=filtfilt(b,a,x(:,i)); % zero-phase delay correction
end

%% BUTTERWORTH FILTER
Lpass = 300;
Hpass = 3000;
%normalize the frequencies
nfBand=[Lpass/(samplingRate/2),Hpass/(samplingRate/2)];
%%
[b,a] = butter(6, nfBand, 'bandpass')
h=fvtool(b,a)


