%% Estimation de V à partir de I^2
clear all
close all
clc

load('pwm.mat')
Ts = 1/fs*(0:(length(ipwm)-1));
Fs = 0:floor(fs/2);

Nx = length(ipwm);
nsc = floor(Nx/150);
nov = floor(nsc/2);

x = ipwm.*ipwm;
y = spwm;

[pxx,fxx] = cpsd(x,x,hanning(nsc),[],[],fs,'twosided');
pxx = abs(pxx);
[pyy,fyy] = cpsd(y,y,hanning(nsc),[],[],fs);
pyy = abs(pyy);
[pyx,fyx] = cpsd(y,x,hanning(nsc),[],[],fs,'twosided');

Tvs = pyx./pxx;
tvs = fftshift(ifft(Tvs,'symmetric'));

% Affichage réponse impulsionnelle
figure
plot(tvs)
title('Réponse impulsionnelle')
yest = conv(x,tvs,'same');

% Affichage des DSP et cohérences spectrales
pyyest = pwelch(yest,hanning(nsc),[],[],fs);
figure
subplot(211)
plot(fyy,mag2db(abs(pyy)))
hold on
plot(fyy,mag2db(abs(pyyest)))
legend('DSP du son','DSP du son estimé')
xlim([0 fyy(end)])
ylim([min(mag2db(abs(pyyest))) max(mag2db(abs(pyyest)))])
hold off

% subplot(212)
% [Cyyest,f] = mscohere(y,yest,hanning(nsc),[],[],fs);
% plot(f,abs(Cyyest))
% hold on
% yyaxis right
% plot(f,20*log10(abs(Cyyest)))
% hold off
% xlim([0 f(end)])

subplot(212)
[Cyyest,f] = mscohere(y,yest,hanning(nsc),[],[],fs);
[Cyx,f] = mscohere(y,x,hanning(nsc),[],[],fs);
plot(f,(abs(Cyx)))
hold on
plot(f,(abs(Cyyest)))
hold off
xlim([0 f(end)])
legend('Cyx','Cyyest')


% Affichage des données temporelles
figure
window = 10000:10500;
plot(Ts(window),y(window))
hold on
plot(Ts(window),yest(window))
xlim([-inf inf])
legend('Son à estimer','Son estimé')
hold off