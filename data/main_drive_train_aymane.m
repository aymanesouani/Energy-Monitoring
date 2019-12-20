
%% 
close all
clear all 
clc

%% Load data
load('pwm.mat');

%% Spectrogram

%% Spectrogram courant
[si support_temps support_frequence]= spectrogram(ipwm, []);
figure(1)
spectrogram(ipwm,'yaxis')
figure(2)
subplot(1,2,1)
imagesc(abs(si))
subplot(1,2,2)
imagesc(angle(si))
 %% Spectrogram son
ss = spectrogram(spwm);
figure(1)
spectrogram(spwm,'yaxis')
figure(2)
subplot(1,2,1)
imagesc(abs(ss))
subplot(1,2,2)
imagesc(angle(ss))

 %% Spectrogram vibrations
sv = spectrogram(vpwm);
figure(1)
spectrogram(spwm,'yaxis')
figure(2)
subplot(1,2,1)
imagesc(abs(sv))
subplot(1,2,2)
imagesc(angle(sv))

%% Spectral analysis : Power spoectral density 

% Courant
% FFT
N = length(ipwm);
idft = fft(ipwm);
idft = idft(1:N/2+1);
psdi =  abs(idft).^2;
psdi(2:end-1) = 2*psdi(2:end-1);

figure(1)
plot(10*log10(psdi))
grid on
title('Periodogram Using FFT ipwm')
xlabel('Freq')
ylabel('Power dB')


taille_fenetre = N/100;
window = hanning(taille_fenetre);
figure(2)
subplot(2,1,1)
plot(ipwm)
title('Representation temporelle du courant')
subplot(2,1,2)
pwelch(ipwm, window)

% Vibrations
% FFT
N = length(vpwm);
vdft = fft(vpwm);
vdft = vdft(1:N/2+1);
psdv =  abs(vdft).^2;
psdv(2:end-1) = 2*psdv(2:end-1);




figure(3)
plot(10*log10(psdv))
grid on
title('Periodogram Using FFT vpwm')
xlabel('Freq')
ylabel('Power dB')


taille_fenetre = N/100;
window = hanning(taille_fenetre);
figure(4)
subplot(2,1,1)
plot(vpwm)
title('Representation temporelle des vibrations')
subplot(2,1,2)
pwelch(vpwm, window)


% Vibrations
% FFT
N = length(spwm);
sdft = fft(spwm);
sdft = sdft(1:N/2+1);
psds =  abs(sdft).^2;
psds(2:end-1) = 2*psds(2:end-1);




figure(5)
plot(10*log10(psds))
grid on
title('Periodogram Using FFT spwm')
xlabel('Freq')
ylabel('Power dB')


taille_fenetre = N/50;
window = hanning(taille_fenetre);
figure(6)
subplot(2,1,1)
plot(spwm)
title('Representation temporelle du son')
subplot(2,1,2)
pwelch(spwm, window)

%% Simu de coherence spectrale

t = 0:1/1e3:2;
signal_entree = chirp(t,0,1,250);

figure(1)
plot(signal_entree)
title('Signal entrée : Chirp')

windowSize = 2; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
x = filter(b,a,signal_entree);

figure(2)
plot(x)
noise = wgn(,power);
signal_bruite =  signal_entree + noise ;
figure(3)
plot(signal_bruite)




