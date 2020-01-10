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
H=cvs.*ps;
%h=real(ifft(H,'symmetric'));
H = [conj(flip(H(2:end)));H];
h_p = ifft(ifftshift(H),'symmetric');


%Reconstruction
test=vpwm(1:length(h_p));
s_pred=conv(vpwm,h_p,'same');
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
plot(20*log(abs(H)))
title("fonction de transfert H")

figure()
subplot(121)
plot(s_pred)
title("son prédit avec v")
subplot(122)
plot(spwm)
title("son mesuré")


%% Modélisation TF et TFI de H ARCHI FAUX C4EST PAS LINEAIRE
clear all
close all
clc

fa=5000;
fs = 100000;
t = 1/fs:1/fs:0.5;
ref = sin(2*pi*fa*t);
bruit = 0.4*sin(2*pi*8*fa*t);
signal = ref + bruit;

h = rand(1000,1);
test=fft(h);
x=signal;

y=conv(x,h,'valid');
y=[y,zeros(1,length(h)-1)];

[Syy, fi]=pwelch(y,hanning(100), [], length(t)/2,fs);
[Sxx, fi]=pwelch(x,hanning(100), [], length(t)/2,fs);

Cxy= mscohere(x,y,hanning(100), [], length(t)/2,fs);

H = Cxy.*Syy;
% H = [conj(flip(H(2:end)));H];
% h_p = ifft(ifftshift(H),'symmetric');




% figure()
% plot(h)
% title("h en tempo")
% figure()
% plot(y)
% title("conv")
% figure()
% plot(h_p)
% title("rep. imp.")
% figure()
% plot(Cxy)

%% Reconstruction d'un signal temporelle à partir de la moitier d'un spectre
clear all
close all
clc

fa=5000;
fs = 100000;
t = 0:1/fs:0.5;
ref = sin(2*pi*fa*t);
bruit = 0.4*sin(2*pi*8*fa*t);
signal = ref + bruit;

S=(fft(ref));
S=S(1:25000);
S = [conj(flip(S(2:end))),S];
s = real(ifft(ifftshift(S),'symmetric'));
mse(s,ref)
figure()
plot(abs(S))
figure()
subplot(121)
plot(s)
title("reconstruit")
subplot(122)
plot(ref)
title("ref")


%% test olivier

clear all
close all
clc

load('pwm.mat')

%Periodogram de welch power spectral density estimate
Ns = length(ipwm);
Ndiv=floor(Ns/100);

[b,a] = ellip(5,0.7,50,0.1);
h=impz(b,a,200);
h=h(15:end);

x = vpwm;

y=conv(x,h,'same');

%perio de welch
[Sxx f1]=pwelch(x,hanning(Ndiv),[],[],fs);
[Syy f2]=pwelch(y,hanning(Ndiv),[],[],fs);


%Coherence spectrale
%Cxy = mscohere(x,y,hanning(Ndiv),[],[],fs);
Cxy = cpsd(x,y,hanning(Ndiv),[],[],fs);
Cyy = cpsd(y,y,hanning(Ndiv),[],[],fs);
Cxx = cpsd(x,x,hanning(Ndiv),[],[],fs);
%Filtrage wiener approx entre vibration et son
%H=Cxy.*Syy;
H=Cxy./Cxx;
H = [conj(flip(H(1:end)));H];
h_p = ifft(ifftshift(H),'symmetric');
%h_p = ifft(H,'symmetric');

%Reconstruction
s_pred=filter(h_p,1,vpwm);
mse = mse(spwm,s_pred)


figure()
plot(f1,abs(Cxy))
title("CSxy")

figure()
subplot(121)
plot(h_p)
title("h predit")
xlim([0 200])
subplot(122)
plot(h)
title("h ref")

figure()
subplot(121)
plot(s_pred)
title("son prédit avec v")
subplot(122)
plot(y)
title("son mesuré")

s_pred=s_pred/max(s_pred);

%soundsc(s_pred,fs)

