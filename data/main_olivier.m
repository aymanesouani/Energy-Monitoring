clear all
close all
clc

load('pwm.mat')
signal = [ipwm,vpwm,spwm];
Ts = 1/fs*(0:(length(ipwm)-1));
Fs = 0:floor(fs/2);

%% Stationnaire ou non ?

Nx = length(ipwm);
nsc = floor(Nx/20);
nov = floor(nsc/2);
figure
for i=1:3
    t = spectrogram(signal(:,i),hanning(nsc),nov);
    subplot(3,1,i)
    imagesc(abs(t))
end

%% DSP

Nx = length(ipwm);
nsc = floor(Nx/10);
nov = floor(nsc/2);
figure
for i=1:3
    pxx = fft(signal(:,i));
    subplot(3,1,i)
    plot(abs(pxx))
end


%% pwelch

Nx = length(ipwm);
nsc = floor(Nx/100);
nov = floor(nsc/2);
figure
for i=1:3
    [pxx,f] = pwelch(signal(:,i),hanning(nsc),[],[],fs);
%     subplot(2,3,i)
%     plot(Ts,ipwm)
    subplot(3,1,i)
    plot(f,mag2db(abs(pxx)))
    xlim([0 f(end)])
end

%% Mise en évidence d'un lien linéaire

Nx = length(ipwm);
nsc = floor(Nx/100);
nov = floor(nsc/2);
figure
for i=1:2
    [Cxy,f] = mscohere(signal(:,i),signal(:,i+1),hanning(nsc),[],[],fs);
    subplot(3,1,i)
    plot(f,abs(Cxy))
    hold on
    yyaxis right
    plot(f,20*log10(abs(Cxy)))
    hold off
    xlim([0 f(end)])

end

%% Estimation de S à partir de V
Nx = length(ipwm);
nsc = floor(Nx/100);
nov = floor(nsc/2);

[pyy,f] = pwelch(signal(:,3),hanning(nsc),[],[],fs);
figure
plot(f,mag2db(abs(pyy)))

% Estimation de S
[Tvs,f] = tfestimate(vpwm,spwm,hanning(nsc),[],[],fs,'Estimator','H2');
tvs = ifft(Tvs,'symmetric');
yest = conv(vpwm,tvs,'same');

pyyest = pwelch(yest,hanning(nsc),[],[],fs);
hold on
plot(f,mag2db(abs(pyyest)))
legend('DSP du son','DSP du son estimée')
xlim([0 f(end)])
ylim([min(mag2db(abs(pyyest))) max(mag2db(abs(pyyest)))])
hold off

figure
window = 1:500;
stem(Ts(window),spwm(window))
hold on
stem(Ts(window),yest(window))
legend('son à estimer','son estimé')
hold off




