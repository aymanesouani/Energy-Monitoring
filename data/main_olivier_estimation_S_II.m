%% Estimation de V � partir de I^2
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

% Affichage r�ponse impulsionnelle
figure
plot(tvs)
title('R�ponse impulsionnelle')
yest = conv(x,tvs,'same');

% Affichage des DSP et coh�rences spectrales
pyyest = pwelch(yest,hanning(nsc),[],[],fs);
figure
subplot(211)
plot(fyy,mag2db(abs(pyy)))
hold on
plot(fyy,mag2db(abs(pyyest)))
legend('DSP du son','DSP du son estim�')
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
plot(f,(abs(Cyyest)))
hold on
plot(f,(abs(Cyx)))
hold off
xlim([0 f(end)])
legend('Cyyest','Cyx')


% Affichage des donn�es temporelles
figure
window = 10000:10500;
stem(Ts(window),y(window))
hold on
stem(Ts(window),yest(window))
xlim([-inf inf])
legend('Son � estimer','Son estim�')
hold off