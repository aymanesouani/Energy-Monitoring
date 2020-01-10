
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

ipwm_carre = ipwm.^2;
figure(1)
subplot(3,1,1)
coherence_icarre_vibration = mscohere(ipwm_carre, vpwm, window, [],[], fs);
plot(coherence_icarre_vibration)
hold on
plot(mscohere(ipwm, vpwm, window, [],[], fs))
title('Coherence spectrale i^2 /vibrations')


subplot(3,1,2)
coherence_i_son = mscohere(ipwm, spwm, window, [],[], fs);
plot(coherence_i_son)
title('Coherence spectrale i/son')

subplot(3,1,3)
coherence_son_vibration = mscohere(spwm, vpwm, window, [],[], fs);
plot(coherence_son_vibration)
title('Coherence spectrale son/vibrations')



%% Estimation de vibrations avec i_carre
Ndiv = 2560;
fenetre = hanning(Ndiv);
Cii = cpsd(ipwm_carre, ipwm_carre,fenetre , [], [], fs);
Civ = cpsd(vpwm, ipwm_carre, fenetre, [], [], fs);
% Filtre
H_iv = Civ ./ Cii ; 
figure(1)
subplot(1,2,1)
plot(abs(H_iv))
title("Module fréquentiel du filtre")
subplot(1,2,2)
plot(180*angle(H_iv)/pi)
title("Argument(H(f)) ")

% Reconstruction signal temporel:
ifft_H_iv = ifft(H_iv);
figure(2)
plot(abs(ifft_H_iv))

% Estimation : 

Vpwm_estime = conv(ipwm_carre, real(ifft_H_iv));
figure(3)
subplot(2,1,1)
plot(vpwm)
title("Vibrations réelles")
subplot(2,1,2)
plot(Vpwm_estime)
title("Vibrations estimées")


Cvv = cpsd(Vpwm_estime, Vpwm_estime,fenetre , [], [], fs);
Csv = cpsd(spwm, Vpwm_estime(1:256000), fenetre, [], [], fs);
% Filtre
H_sv = Csv ./ Cvv ; 
figure(4)
subplot(1,2,1)
plot(abs(H_sv))
title("Module fréquentiel du filtre")
subplot(1,2,2)
plot(180*angle(H_sv)/pi)
title("Argument(H(f)) ")

% Reconstruction signal temporel:
ifft_H_sv = ifft(H_sv);
figure(5)
plot(abs(ifft_H_sv))

% Estimation : 

Spwm_estime = conv(Vpwm_estime, real(ifft_H_sv));
figure(6)
subplot(2,1,1)
plot(spwm)
title("Son réel")
subplot(2,1,2)
plot(Spwm_estime)
title("Son estimé")


% Plot des DSP :
figure(7)
subplot(2,1,1)
plot(pwelch(vpwm, fenetre))
hold on
plot(pwelch(Vpwm_estime, fenetre))

subplot(2,1,2)
plot(pwelch(spwm, fenetre))
hold on
plot(pwelch(Spwm_estime, fenetre))

figure(8)
subplot(2,1,1)
plot(pwelch(vpwm, fenetre))
hold on
plot(pwelch(Vpwm_estime, fenetre))

subplot(2,1,2)
plot(pwelch(spwm, fenetre))
hold on
plot(pwelch(Spwm_estime, fenetre))





