clear;
close all;
PLOT = 1;
% data given in assignment
N=9;data = [1 0 0 1 0 1 0 1 0];
% for random data generation
% N = 10000;data = randi([0 1], 1, N);

% set symbol time / symbol duration
ts = 1;
% set carrier details
fc = 1;
Ac = 2;
% set sampling frequency
fs = 100; % should be > 1/ts
% set SNR
SNRdB = 5;
% SNRdB = 10;
% SNRdB = 15;
% SNRdB = 20;
% SNRdB = 25;

SNR = 10^(SNRdB/10); %for voltage
BER_th = (1/2)*erfc(sqrt(SNR));

sample_per_bit = ts*fs;

t = (1/fs:1/fs:N*ts);
message = repelem(2*data - 1, sample_per_bit); %input sequence as a square wave
carrier = Ac * sin(2*pi*fc*t); %carrier wave
mod = message .* carrier; %modulation

rx = awgn(mod, 10*log10(SNR/(.5*(sample_per_bit))),'measured'); % AWGN noise addition
rx_bef = rx; %transmitted wave before filter.
rx = lowpass(rx, 2, fs); %transmitted wave after filter.

demod_bef = rx_bef.*carrier; %demodulated wave before filter.
demod = rx .* carrier; %demodulated wave after filter.
    
data_rx = zeros(1, N); %output data after filter
corr = zeros(1, sample_per_bit*N);
for i=1:N
    corr((i-1)*sample_per_bit+1:i*sample_per_bit) = trapz(t((i-1)*sample_per_bit+1:i*sample_per_bit), demod((i-1)*sample_per_bit+1:i*sample_per_bit));
    data_rx(i) = sum(corr((i-1)*sample_per_bit+1:i*sample_per_bit)) > 0;
end

data_rx_bef = zeros(1, N); %output data before filter
corr_bef = zeros(1, sample_per_bit*N);
for i=1:N
    corr_bef((i-1)*sample_per_bit+1:i*sample_per_bit) = trapz(t((i-1)*sample_per_bit+1:i*sample_per_bit), demod_bef((i-1)*sample_per_bit+1:i*sample_per_bit));
    data_rx_bef(i) = sum(corr_bef((i-1)*sample_per_bit+1:i*sample_per_bit)) > 0;
end

%Plotting of graphs
if PLOT
    tiledlayout(5, 1);
    ax = nexttile;
    plot(ax, t, mod);
    title(ax, 'Modulated signal');
    xlabel('Time (seconds)-->')
    ylabel('Amplitude-->')
    
    ax = nexttile;
    plot(ax, t, rx_bef);
    title(ax, 'Received signal (before filter)');
    xlabel('Time (seconds)-->')
    ylabel('Amplitude-->')
    
    ax = nexttile;
    plot(ax, t, rx);
    title(ax, 'Received signal (after filter)');
    xlabel('Time (seconds)-->')
    ylabel('Amplitude-->')
    
    ax = nexttile;
    plot(ax, t, demod_bef);
    hold on;
    plot(ax, t, corr_bef);
    title(ax, 'Demodulated signal (before filter)');
    xlabel('Time (seconds)-->')
    ylabel('Amplitude-->')
    hold off;   
    
    ax = nexttile;
    plot(ax, t, demod);
    hold on;
    plot(ax, t, corr);
    title(ax, 'Demodulated signal (after filter)');
    xlabel('Time (seconds)-->')
    ylabel('Amplitude-->')
    hold off;
    
%     figure(2);
%     plot(t, message);
%     title('Input Sequence');
%     xlabel('Time (seconds)-->')
%     ylabel('Amplitude-->')
%     hold on;
%     plot(t, corr);
%     legend('sent data', 'decoded data');
%     hold off;
end

%BER calculation
[nerr, ber] = biterr(data_rx,data);
ber