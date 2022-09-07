%% Clearing Everything Before Initialization
clear;
close all;
clc;
Trange = [0 5 -2 2];

% Here the user shall be given the choice to
B=5*10^3;

radio_file = menu('Choose the radio signal','Signal_1','Signal_2','Signal_3');
if radio_file == 1
[y, fs] = audioread('Signal_1.wav');
info = audioinfo('Signal_1.wav');
end
if radio_file == 2
[y, fs] = audioread('Signal_2.wav');
info = audioinfo('Signal_2.wav');
end
if radio_file == 3
[y, fs] = audioread('Signal_3.wav');
info = audioinfo('Signal_3.wav');
end

ts=1/fs;                                             %sampling time              
t=[0 : ts : info.Duration]; t=t(1:end-1);            %time period
m = transpose(y(:,1));
Lfft=length(t);
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
figure(1);
subplot(3,2,1);
axis(Trange);
plot(t,m);
title('message signal (t-domain)');


M_s=fftshift(fft(m, Lfft));
Frange=[-20*10^3 20*10^3 -1*10^3 1*10^3];
subplot(3,2,2);
fd1 = plot(freqs, abs(M_s));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');


%% Carrier
%In broadcast applications, for instance, a transmit antenna can radiate only a narrow band without
% distortion. This means that to avoid distortion caused by the transmit antenna, we must have
% fc/B>>1. Thus broadcast band AM radio, with B = 5 kHz and the band of 550:1600 kHz for fc
fc=550*10^3;                  %carrier frequency
c=cos(2*pi*fc*t);            %carrier function


%% Modulation
s_dsb=m.*(cos(2*pi*fc*t));

subplot(3,2,3);

plot(t,s_dsb);
axis(Trange);
title('modulated signal (t-domain)');

S_dsb=fftshift(fft(s_dsb,Lfft));

Frange=[-20*10^3 20*10^3 -1*10^3 1*10^3];
subplot(3,2,4);
fd1 = plot(freqs, abs(S_dsb));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');
%% Demoduation
s_dem = s_dsb .* (2*cos(2*pi*fc*t));
% low pass filter
cutoff_freq = 8000;
cutoff_coefficient = cutoff_freq/(fs/2);
[b,a]=butter(5, cutoff_coefficient);
s_dem=filter(b,a,s_dem);
S_dem = fftshift(fft(s_dem, Lfft));


subplot(3,2,5);
axis(Trange);
plot(s_dem);
title('De-modulated with receiver signal (t-domain)');

DEM_s=fftshift(fft(S_dem, Lfft));
Frange=[-20*10^3 20*10^3 -1*10^3 1*10^3];
subplot(3,2,6);
fd1 = plot(freqs, abs(S_dem));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('De-modulated signal spectrum (f-domain)');


%% Superheterodyne receiver
% if >= B
% fc >= if + B
f_if = 445*10^3;
Fs = 200*10^4;

% plotting
figure(2);
% plot message
subplot(6,2,1);
plot(t,m);
axis(Trange);
title('message signal (t-domain)');


M_s=fftshift(fft(m, Lfft));
Frange=[-2*10^3 2*10^3 -1*10^3 1*10^3];
subplot(6,2,2);
fd1 = plot(freqs, abs(M_s));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum (f-domain)');

% plot modulated
subplot(6,2,3);
plot(t,s_dsb);
axis(Trange);
title('modulated signal (t-domain)');

Frange=[-20*10^3 20*10^3 -1*10^3 1*10^3];
subplot(6,2,4);
fd1 = plot(freqs, abs(S_dsb));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('modulated signal spectrum (f-domain)');

%% RF filter
%if bandpass, then, F_up < fc+2If, F_lower > fc-B_band
L_cutoff_freq=(fc-f_if);
U_cuttoff_freq=(fc+f_if);
lrf_limit_coefficient= L_cutoff_freq/(Fs/2);
urf_limit_coefficient= U_cuttoff_freq/(Fs/2);
[brf,arf]=butter(5,[lrf_limit_coefficient urf_limit_coefficient]);
rff=3.*filter(brf,arf,s_dsb);

subplot(6,2,5);
axis(Trange);
plot(rff);
title('RF Filtered signal (t-domain)');

RFF_s=fftshift(fft(rff, Lfft));
Frange=[-20*10^3 20*10^3 -1*10^3 1*10^3];
subplot(6,2,6);
fd1 = plot(freqs, abs(RFF_s));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('RF Filtered signal spectrum (f-domain)');

%% mixer
%       Local Oscillator Parameters
%-----------------------------------------------

f_lo = fc+f_if;
f_gen=cos(2*pi*f_lo*t);
mixed=rff.*f_gen;

subplot(6,2,7);
plot(t, mixed);
axis(Trange);
title('Mixed signal (t-domain)');

MIXED_s = fftshift(fft(mixed,Lfft));

subplot(6,2,8)
title('Response after mixing signals ');
xlabel('Frequency in Hz');
ylabel('Amplitude');
axis([-10^6 10^6 -1000 1000]);
plot(freqs, abs(MIXED_s));
title('Mixed signal spectrum (f-domain)');

%% IF bandpass filter
lif_cutoff_cof= (f_if-B/2)/(Fs/2); %removing lower adjacent station
uif_cutoff_cof= (f_if+B/2)/(Fs/2); %removing upper adjacent station
[bif, aif] = butter(6,[lif_cutoff_cof, uif_cutoff_cof]);
iff = 3*filter(bif,aif,mixed); %Remove the neighboring stations

subplot(6,2,9);
axis(Trange);
plot(t, (iff));
title('IF Filtered signal (t-domain)');

%cheking in frequency domain
IFF_s = fftshift(fft(iff, Lfft)); %Freq Response of Filtered Stations

subplot(6,2,10)
plot(freqs, abs(IFF_s));
title('IF Filter Response (f-domain) ')
axis([-20*10^3 20*10^3 -1000 1000]);
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% receiver demodulation
s_recdem = iff .* (2*cos(2*pi*f_if*t));
% low pass filter
cutoff_freq = 8000;
cutoff_coefficient = cutoff_freq/(fs/2);
[b,a]=butter(5, cutoff_coefficient);
s_recdem=filter(b,a,s_dem);

S_recdem = fftshift(fft(s_recdem, Lfft));

sound(s_recdem, fs);

subplot(6,2,11);
plot(t, (s_recdem));
axis(Trange);
title('Demodulated signal with receiver signal (t-domain)');

%cheking in frequency domain, verifying fft
subplot(6,2,12);
plot(freqs, abs(S_recdem));
title('Demodulated signal with receiver signal (f-domain)');
axis([-20*10^3 20*10^3 -1000 1000]);
xlabel('Frequency (Hz;)');
ylabel('Magnitude');