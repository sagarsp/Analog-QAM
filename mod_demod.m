%% Project
%
% Author : Sagar Patel
% 
% Aim : To implement the analog Quadrature amlitude modultaion scheme on stereo.
% Here we approch, modulating two signals i.e. left and right channel signal of stereo
% signal, with same frequency of sinusoid. But one of the signal will get
% modulated by sine and other with cosine. As the sine and cosine are of
% pi/2 phase, the corelation (i.e. integreation over a period) is zero. So
% to use this phenomena, both the modulated signal is added with
% each-other and then we will transmit them. Here to give the glimps of
% analog signal we are using audio in wav file format, which is loss-less
% format among the all audio compression standard. Though, now-a-days the
% audio becomes digital. So one more thing to do is upsample the signal, so
% that it will look like upsampled(20-20kHz upsampled by 10. so 2-2kHz) 
% analog singnal. And then we will modulate both the singal by 16kHz
% sinusiod wave. Then add it. And transmitting the modulated singal. 

% On the other hand, we receive the signal which is sent through the
% channel, maybe wireless or wired depeding upon the choice w.r.t economy or
% environment(beheaviour). So there will be some noise received at the
% receiver end. But for simplicity, to understand the modulation and
% demodulation process, assuming the nosie on the receiver end is zero i.e.
% ideal condition. Now, To demodulate the signal, we pass the received
% signal in two parallel circuits. In which, one mixes the sine wave and
% other mixes the cosine wave by mixer. Here also assumed that the carreir
% at the receiver side is in sync. with the received signal carreir. We can
% extract carrier from the received signal also. Now the next thing is to
% get the left and right channel signal back by passing the output of
% mixer into low-pass filter at the receiver end. Here we are using 
% ButterWorth Low-Pass Filter. Then we will downsample the signals and we
% will get both the singal. If we are modulating stereo, then at the other 
% end we need to get the stereo sound back by formating the two signals at
% the received end.  

%% Initilizing and Reading the data
clear all;
clc;
close all;
clf;
[signal fs]=wavread('sample.wav');
Time=length(signal)/fs;  % In Seconds
%% Left and Right Audio Signals
s_l=signal(:,1)';
s_r=signal(:,2)';
%% Plotting the Data
figure(1);
t=0:1/fs:Time-(1/fs);
subplot(321);
plot(t',signal);
title('Stereo Signal');
xlabel('Time (seconds)');ylabel('Amplitude');
subplot(322);
plot((50000:51000)/fs,signal(50000:51000,:));
title('Stereo Signal (ZOOM UP) <50000-51000>');
xlabel('Time (seconds)');ylabel('Amplitude');
subplot(323);
plot(t',s_l);
title('Left Signal');
xlabel('Time (seconds)');ylabel('Amplitude');
subplot(324);
plot((50000:51000)/fs,s_l(50000:51000));
title('Left Signal (ZOOM UP) <50000-51000>');
xlabel('Time (seconds)');ylabel('Amplitude');
subplot(325);
plot(t',s_r);
title('Right channel Signal');
xlabel('Time (seconds)');ylabel('Amplitude');
subplot(326);
plot((50000:51000)/fs,s_r(50000:51000));
title('Right Signal (ZOOM UP) <50000-51000>');
xlabel('Time (seconds)');ylabel('Amplitude');
%% Listening the Sound
display('Lets listen to the stereo Sound <press any key>');
pause;
sound(signal,fs);
display('Stereo Sound');
pause;
display('Lets listen to left channel signal <press any key>');
pause;
sound([s_l',zeros(length(s_l),1)],fs);
display('Left sound');
pause;
display('Lets listen to right channel signal <press any key>');
pause;
sound([zeros(length(s_r),1),s_r'],fs);
display('Right sound');
pause;
%% Up Sampling to Modulate the Signal
L = 10;         % Up-Sampling Co-Efficient
x1 = s_l;
x2 = s_r;
% Generating the interpolated output sequence
y1 = interp(x1,L);
y2 = interp(x2,L);
% Ploting the input (left) and the output (upsampled left) sequences
figure(2);
subplot(2,1,1);
n=10000:10100;
stem(n,x1(n(1):n(length(n))));
title('Left Sequence <ZOOM UP>');
xlabel('Time index n'); ylabel('Amplitude');
subplot(2,1,2);
m = L*n(1):(n(length(n))*L)-1;
stem(m,y1(m(1):m(length(m))));
title('Up Sampled Left Sequence <ZOOM UP>');
xlabel('Time index n'); ylabel('Amplitude');
s_l_up = y1;
s_r_up = y2;
%% Generating a carrier
fc = 16000;
A = 1;
t = 0:1/fs:(length(s_l_up)-1)/fs;
carrier = A*cos(2*pi*fc*t);
carrier2 = A*sin(2*pi*fc*t);
figure(3);
subplot(221);
plot(t,carrier);
title('Carrier Wave - Cosine (16 kHz)');
xlabel('Time (seconds)');
ylabel('Amplitude (seconds)');
subplot(222);
plot(t(100:150),carrier(100:150));
title('Carrier Wave - Cosine <ZOOM UP>');
xlabel('Time (seconds)');
ylabel('Amplitude (seconds)');
subplot(223);
plot(t,carrier2);
title('Carrier Wave - Sine (16 kHz)');
xlabel('Time (seconds)');
ylabel('Amplitude (seconds)');
subplot(224);
plot(t(100:150),carrier2(100:150));
title('Carrier Wave - Sine <ZOOM UP>');
xlabel('Time (seconds)');
ylabel('Amplitude (seconds)');

figure(4);
freqz(carrier,1,512,fs);
title('Frequency Response of Carrier - Cosine');
figure(5);
freqz(carrier2,1,512,fs);
title('Frequency Response of Carrier - Sine');
%% Modulation the Snare signal with Carrier
mod1 = s_l_up.* carrier;
mod2 = s_r_up.* carrier2;
mod = mod1 + mod2;
%% Plotting in Time domain
figure(6);
subplot(311);
plot(s_l_up(104000:105000));
hold on;
plot(s_r_up(104000:105000),'-g');
title('Left(B) and Right(G) in time domain');
subplot(312);
plot(carrier(104000:105000));
title('Carrier in time domain');
subplot(313);
plot(mod1(104000:105000));
hold on;
plot(mod2(104000:105000),'-g');
title('Modulated - Left(B) Right(G) - in time domain');

figure(7);
subplot(311);
plot(mod1(104000:105000));
title('Left Modulated in time domain');
subplot(312);
plot(mod2(104000:105000));
title('Right Modulated in time domain');
subplot(313);
plot(mod(104000:105000));
title('Modulated signal in time domain');
%% Plotting the FFT of Signals
figure(8);
NFFT = 65536*2*2;
wd = [-pi:2*pi/(NFFT-1):pi]; %digital frequency in rad/s
fd = [-pi:2*pi/(NFFT-1):pi]/(2*pi); %digital frequency in hz
Fa = fd * fs; %analog freq in hz

subplot(321);
plot(Fa,abs(fftshift(fft(s_l_up,NFFT))));
title('FFT of Up-Sampled of Left signal');
subplot(322);
plot(Fa,abs(fftshift(fft(s_r_up,NFFT))));
title('FFT of Up-Sampled of Right signal');
subplot(323);
plot(Fa,abs(fftshift(fft(carrier,NFFT))));
title('FFT of Carrier Signal - Cosine');
subplot(324);
plot(Fa,abs(fftshift(fft(carrier2,NFFT))));
title('FFT of Carrier Signal - Sine');
subplot(325);
plot(Fa,abs(fftshift(fft(mod1,NFFT))));
title('FFT of modulated signal');
subplot(326);
plot(Fa,abs(fftshift(fft(mod2,NFFT))));
title('FFT of modulated signal');
%% Demodulation
rec=mod;
%Here we can generate carrier locally or by extracting from signal
%So using the transmitted carreir to co-relate
rec1 = rec.*carrier;
rec2 = rec.*carrier2;
%% Plotting the FFt of Demodulated rec1 and rec2
figure(9);
subplot(311);
plot(Fa,abs(fftshift(fft(rec,NFFT))));
title('Received Signal Spectrum');
subplot(312);
plot(Fa,abs(fftshift(fft(rec1,NFFT))));
title('rec1 Signal Spectrum');
subplot(313);
plot(Fa,abs(fftshift(fft(rec2,NFFT))));
title('rec2 Signal Spectrum');

%% Generating a Low Pass Filter
F_Pass= 4000;    %PassBand Freq = 4kHz
F_Stop= 10000;   %StopBand Freq = 10kHz

F_s=fs;

fd_p = F_Pass / F_s;
fd_s = F_Stop / F_s;

Rp = 0.5;
Rs = 30;
[Order, Cut_off] = BUTTORD(fd_p, fd_s, Rp, Rs);
[num,den] = BUTTER(Order,Cut_off);

%% Demodulating the signal by extracting the carrier
clear mod1 mod2;
rec_1_Filtrd = filter(num,den,rec1);
rec_2_Filtrd = filter(num,den,rec2);
figure(10);
[H,F] = FREQZ(num,den,512,F_s*2);
plot(F,abs(H));
title('Frequency Respose of LowPass Filter');
%% Plotting the FFT
figure(11);
subplot(311);
plot(Fa,abs(fftshift(fft(rec,NFFT))));
title('Received Signal Spectrum');
subplot(312);
plot(Fa,abs(fftshift(fft(rec1,NFFT))));
hold on;
f_vary = [-flipud(F); F];
h_vary = [flipud(abs(H)); abs(H)];
plot(f_vary,2000*h_vary,'r');
title('Received 1 Signal Spectrum and LowPass Response');
subplot(313);
plot(Fa,abs(fftshift(fft(rec_1_Filtrd,NFFT))));
title('Low Pass Filtered Signal Spectrum');
%% Downsampling Data
lft_dwn = rec_1_Filtrd(1:L:end);
rht_dwn = rec_2_Filtrd(1:L:end);
%% Generating Stereo
final = [lft_dwn' rht_dwn'];
%% Plotting Received Downsampled signals in Time domain
figure(12);
t=0:1/fs:Time-(1/fs);
subplot(311);
plot(t,lft_dwn);
title('Left Signal in Time Domain');
xlabel('Time (Seconds)');ylabel('Amplitude');
subplot(312);
plot(t,rht_dwn);
title('Right Signal in Time Domain');
xlabel('Time (Seconds)');ylabel('Amplitude');
subplot(313);
plot(t',final);
title('Stereo Signal (Left = B and Right = G) in Time Domain');
xlabel('Time (Seconds)');ylabel('Amplitude');
%% Plotting Error Signal (We have not included AWGN. So Error = 0)
E_l = s_l - lft_dwn;
E_r = s_r - rht_dwn;
E_str = signal - final;
figure(13);
t=0:1/fs:Time-(1/fs);
subplot(311);
plot(t,E_l);
title('Error in Left Signal in Time Domain');
xlabel('Time (Seconds)');ylabel('Amplitude');
subplot(312);
plot(t,E_r);
title('Error in Right Signal in Time Domain');
xlabel('Time (Seconds)');ylabel('Amplitude');
subplot(313);
plot(t',E_str);
title('Error in Stereo Signal (Left = B and Right = G) in Time Domain');
xlabel('Time (Seconds)');ylabel('Amplitude');
%% Listening to the Saperate channels
display('Left Received and Down-sampled Signal <press any key>');
pause;
sound([lft_dwn' zeros(length(lft_dwn),1)],fs);
display('Left Signal');
pause;
display('Right Received and Down-sampled Signal <press any key>');
pause;
sound([zeros(length(rht_dwn),1) rht_dwn'],fs);
display('Right Signal');
pause;
%% Playing a stereo Signal
display('Received Stereo Signal <press any key>');
pause;
sound(final,fs);
display('Received Signal');
pause;
