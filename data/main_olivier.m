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

%% Mise en �vidence d'un lien lin�aire

Nx = length(ipwm);
nsc = floor(Nx/100);
nov = floor(nsc/2);
figure
for i=1:2
    [Cxy,f] = mscohere(signal(:,i),signal(:,i+1),hanning(nsc),[],[],fs);
    subplot(2,1,i)
    plot(f,abs(Cxy))
    hold on
    yyaxis right
    plot(f,20*log10(abs(Cxy)))
    hold off
    xlim([0 f(end)])

end

%% Estimation de S � partir de V
clear all
close all
clc

load('pwm.mat')
signal = [ipwm,vpwm,spwm];
Ts = 1/fs*(0:(length(ipwm)-1));
Fs = 0:floor(fs/2);

Nx = length(ipwm);
nsc = floor(Nx/150);
nov = floor(nsc/2);

% x = signal(:,1).*signal(:,1);
x = signal(:,2);
y = signal(:,3);

[pxx,fxx] = cpsd(x,x,hanning(nsc),[],[],fs,'twosided');
pxx = abs(pxx);
[pyy,fyy] = cpsd(y,y,hanning(nsc),[],[],fs);
pyy = abs(pyy);
[pyx,fyx] = cpsd(y,x,hanning(nsc),[],[],fs,'twosided');

Tvs = pyx./pxx;
tvs = fftshift(ifft(Tvs,'symmetric'));

% tvs = flip(fftshift(ifft(ifftshift(Tvs),'symmetric')));
% tvs = tvs(1:floor(length(tvs)/2)+1);
figure
plot(tvs)
title('R�ponse impulsionnelle')
% yest = conv(tvs,vpwm);
yest = conv(vpwm,tvs,'same');


pyyest = pwelch(yest,hanning(nsc),[],[],fs);
figure
subplot(211)
plot(fyy,mag2db(abs(pyy)))
hold on
plot(fyy,mag2db(abs(pyyest)))
legend('DSP du son','DSP du son estim�e')
xlim([0 fyy(end)])
ylim([min(mag2db(abs(pyyest))) max(mag2db(abs(pyyest)))])
hold off

subplot(212)
[Cyyest,f] = mscohere(y,yest,hanning(nsc),[],[],fs);
plot(f,abs(Cyyest))
hold on
yyaxis right
plot(f,20*log10(abs(Cyyest)))
hold off
xlim([0 f(end)])

figure
window = 10000:20000;
stem(Ts(window),spwm(window))
hold on
stem(Ts(window),yest(window))
xlim([-inf inf])
legend('son � estimer','son estim�')
hold off

%% Mod�lisation TF et TFI de H
clear all
close all
clc


load('pwm.mat')
signal = [ipwm,vpwm,spwm];
Ts = 1/fs*(0:(length(ipwm)-1));
Fs = 0:floor(fs/2);

[b,a] = ellip(5,0.7,50,0.1);
h=impz(b,a,50);
h=h(1:end);

x = vpwm;

y=conv(x,h,'same');

Nx = length(x);
nsc = floor(Nx/200);
nov = floor(nsc/2);

[pxx,fxx] = cpsd(x,x,hanning(nsc),[],[],fs,'twosided');
pxx = abs(pxx);
[pyy,fyy] = cpsd(y,y,hanning(nsc),[],[],fs);
pyy = abs(pyy);
[pyx,fyx] = cpsd(y,x,hanning(nsc),[],[],fs,'twosided');
Tvs = pyx./pxx;

tvs = fftshift(ifft((Tvs),'symmetric'));
% tvs = fftshift(real(ifft((Tvs))));

% Tvs = [conj(flip(Tvs(2:end)));Tvs];
% tvs = ifft(ifftshift(Tvs),'symmetric');
% tvs = flip(fftshift(ifft(ifftshift(Tvs),'symmetric')));
% tvs = tvs(floor(length(tvs)/2):end);
figure
plot(tvs)
hold on
plot(h)
legend('estim�e','vraie')
title('R�ponse impulsionnelle')

yest = conv(x,tvs,'same');


pyyest = pwelch(yest,hanning(nsc),[],[],fs);
figure
plot(fyy,mag2db(abs(pyy)))
hold on
plot(fyy,mag2db(abs(pyyest)))
legend('DSP du son','DSP du son estim�e')
xlim([0 fyy(end)])
ylim([min(mag2db(abs(pyyest))) max(mag2db(abs(pyyest)))])
hold off

figure
window = 50000:50500;
stem(y(window))
hold on
stem(yest(window))
legend('son � estimer','son estim�')
hold off
