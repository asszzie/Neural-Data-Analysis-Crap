load NDA_rawdata
%select a subset of the data called xx for each channel separatedly
x  = gain * double(x);
xx = x([100:200000],1);
xx2 = x([100:200000],2);
xx3 = x([100:200000],3);
xx4 = x([100:200000],4);

% define a bandpass filter (300-3000 Hz) 
Fs = samplingRate;
Lpass=300;
Hstop=3000;
fBand=[Lpass/(Fs/2),Hstop/(Fs/2)];% normalization of frequency limits
[b,a]=cheby2(7,50,fBand,'bandpass'); % create the chebyshev type II bandpass filter
%h=fvtool(b,a) %visualize the filter if needed

% apply the filter described by a and b with zero-phase delay correction
y=filtfilt(b,a,xx); 
y2=filtfilt(b,a,xx2); 
y3=filtfilt(b,a,xx3); 
y4=filtfilt(b,a,xx4); 

%-------------------figure 1-------------------------
% create all subplots of raw and filtered data subset 
subplot(4,2,1);
plot(xx)
title('Channel 1 raw', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,2);
plot(y)
title('Channel 1 filtered', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,3);
plot(xx2)
title('Channel 2 raw', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,4);
plot(y2)
title('Channel 2 filtered', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,5);
plot(xx3)
title('Channel 3 raw', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,6);
plot(y3)
title('Channel 3 filtered', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,7);
plot(xx4)
title('Channel 4 raw', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -400 400])

subplot(4,2,8);
plot(y4)
title('Channel 4 filtered', 'FontSize', 16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize', 11)
axis([0 2*10^5 -400 400])
%saveas(3,'Figure 1 - week 1-version2.png')

%----------------------------figure 2-------------------------------
%Determine Threshold and detect spikes for subset data in 'seg' 
seg = [y, y2, y3, y4]; % stores segment of data to plot with 4 channels

[s,t,spikes] = detectSpikes2(seg, Fs); 

[p,l] = findpeaks(-s); % p is not important, l is collecting the first position of each column
l = [1;l]; % add 1 to include the first channel
s1 = s([l(1):l(2)-1]);    % first channel with sample location of spikes
s2 = s([l(2):l(3)-1]);     % second channel with sample location of spikes etc.
s3 = s([l(3):l(4)-1]);
s4 = s([l(4):length(s)]);
close all
%Make figure 2 for channel 1
sigma=median(abs(seg)./0.6745); % sigma is the estimated standard deviation 
Threshold=-5*sigma;

figure;
plot(y)
hold on
plot(get(gca,'xlim'), [Threshold(:,1) Threshold(:,1)], 'r'); %draw threshold line using the size x-axis

Spike_ch1 =spikes(:,1);
figure(1)
hold on
plot(find(Spike_ch1<0), Spike_ch1(Spike_ch1<0), '^r', 'MarkerSize', 4)
title('Spikes below Threshold \newline          Channel 1', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -140 100])
box off

%idem dito for channel 2 to 4
figure;
plot(y2)
hold on
plot(get(gca,'xlim'), [Threshold(:,2) Threshold(:,2)], 'r'); %draw threshold line using the size x-axis
Spike_ch2 =spikes(:,2);
figure(2)
hold on
plot(find(Spike_ch2<0), Spike_ch2(Spike_ch2<0), '^r', 'MarkerSize', 4)
title('Spikes below Threshold \newline          Channel 2', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -140 100])
box off
hold off

figure;
plot(y3)
hold on
plot(get(gca,'xlim'), [Threshold(:,3) Threshold(:,3)], 'r'); %draw threshold line using the size x-axis
Spike_ch3 =spikes(:,3);
figure(3)
hold on
plot(find(Spike_ch3<0), Spike_ch3(Spike_ch3<0), '^r', 'MarkerSize', 4)
title('Spikes below Threshold \newline          Channel 3', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -140 100])
box off
hold off

figure;
plot(y4)
hold on
plot(get(gca,'xlim'), [Threshold(:,4) Threshold(:,4)], 'r'); %draw threshold line using the size x-axis
Spike_ch4 =spikes(:,4);
figure(4)
hold on
plot(find(Spike_ch4<0), Spike_ch4(Spike_ch4<0), '^r', 'MarkerSize', 4)
title('Spikes below Threshold \newline          Channel 4', 'FontSize',16)
xlabel('Time (ms)','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
axis([0 2*10^5 -140 100])
box off
hold off

%************************WATCH OUT! FOR FIGURE 3 WE USE ALL RECORDED
%DATA!!************************************************************
%%
load NDA_rawdata
x  = gain * double(x);
% run code
y = filterSignal(x,samplingRate);
[s, t, spikes] = detectSpikes2(y,samplingRate);

%Waveform extraction
Waveforms_Spikes = nan(30,4,length(s)); %create a matrix that has a 30x4xSpikes appearance 
Spike_values = [s y(s)];% Make a Matrix with index number and values of spikes at the time of crossing

for i = 1:length(s)%go through all spikes and put them in the mentioned matrix

    Waveforms_Spikes(:,:,i) = y(s(i)-15:s(i)+14,:);
    
end
%-------------figure 3a: the 100 first spikes--------------
figure
hold on
wave=reshape(permute(w,[2,1,3]),30,55909*4);
for j = 1:100 %plot the first 100 spikes.
   plot(wave(:,j))
   hold on
   drawnow update
end
title('Waveform of the Hunderd First Spikes (x-axis corresponds to 1 ms) \newlineAll channels', 'FontSize',16)
xlabel('Samples','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
hold off

%--------------figure 3b: the 100 biggest spikes-------------
all_spikes = reshape(spikes,1,length(spikes)*4); % make row vector to sort
biggest = sort(all_spikes);
biggest_100 = biggest([1:100]);

%recollect the indices of the 100 biggest spikes in original spike vector
index_100=[];
for i=1:100
    index_100(i) = find(spikes==biggest_100(i));
end

%draw 100 biggest spikes (waveform)
figure;
Waveform_100=[];
for i = 1:100 %go through 
    Waveforms_100(:,:,i) = y(index_100(i)-15:index_100(i)+14);
    plot(Waveforms_Spikes(:,:,i)) 
    hold on
    drawnow update
end
title('Waveform of the Hunderd Biggest Spikes (x-axis corresponds to 1 ms) \newlineAll channels', 'FontSize',16)
xlabel('Samples','FontSize',11)
ylabel('Voltage (\muV)', 'FontSize',11)
hold off
