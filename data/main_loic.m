clear all
close all
clc

load('pwm.mat')


%Periodogram de welch power spectral density estimate
Ns = length(ipwm);
Ndiv=floor(Ns/100);
%[p2,f2]=pwelsh(signal,hanning(Ns),max(256,2^nextpow2(Ns/10)),fs);
[pi fi]=pwelch(ipwm,hanning(Ndiv),[],[],fs);
[ps fS]=pwelch(spwm,hanning(Ndiv),[],[],fs);
[pv fv]=pwelch(vpwm,hanning(Ndiv),[],[],fs);

%Coherence spectrale
civ = mscohere(ipwm,vpwm,hanning(Ndiv),[],[],fs);
cvs = mscohere(vpwm,spwm,hanning(Ndiv),[],[],fs);
cis = mscohere(ipwm,spwm,hanning(Ndiv),[],[],fs);

%Filtrage wiener approx entre vibration et son
pv_norm = pv./pv;
H=cvs./pv_norm;
h=real(ifft(H,'symmetric'));


%Reconstruction
test=vpwm(1:2049);
s_pred=conv(vpwm,h,'same');
mse = mse(spwm,s_pred)

figure()
subplot(311)
plot(ipwm)
title("courant en t")
subplot(312)

plot(spwm)
title("son en t")
subplot(313)
plot(vpwm)
title("vibration en t")

figure()
subplot(311)
plot(fi,20*log(abs(pi)))
title("DSP courant")
subplot(312)
plot(fS,20*log(abs(ps)))
title("DSP son")
subplot(313)
plot(fv,20*log(abs(pv)))
title("DSP vibration")

figure()
subplot(311)
plot(fv,civ)
title("CS courant/vibration")
subplot(312)
plot(fS,cvs)
title("CS vibration/son")
subplot(313)
plot(fv,cis)
title("CS courant/son")

figure()
plot(fS,20*log(abs(H)))
title("fonction de transfert H")

figure()
subplot(121)
plot(s_pred)
title("son prédit avec v")
subplot(122)
plot(s_pred)
title("son mesuré")


%% Modélisation TF et TFI de H
clear all
close all
clc

h=zeros(500,1);
h(1:250)=1;
h(251:500)=-1;

x = zeros(2000,1);
x(1000)=1;

y=conv(x,h,'same');
[Syy, fi]=pwelch(y,hanning(100));
[Sxx, fi]=pwelch(x,hanning(100));

Cxy= mscohere(x,y,hanning(100));
Sxx = Sx

H = Cxy./Sxx;
fa=5000;
fs = 100000;
t = 0:1/fs:0.05;
ref = sin(2*pi*fa*t);
bruit = 0.4*sin(2*pi*8*fa*t);
signal = ref + bruit;


figure()
plot(h)
title("h en tempo")
figure()
plot(y)
title("conv")

% figure()
% plot(signal)