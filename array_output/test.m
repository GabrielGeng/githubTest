
n = (1:48000).';
fs = 16000;
s1 = 0.999.^n.*ones(48000,1);
s1 = s1.*exp(1j*2*pi*680/fs*n);
nfft = length(s1);
f = fs*(0:nfft/2)/nfft;
S_1 = fft(s1,nfft);
figure
plot(f,abs(S_1(1:nfft/2+1)));
%soundsc(real(s1),fs)
s2 = zeros(48000,1);
s2(1:1000) = 1;
s2 = s2.*exp(1j*2*pi*680/fs*n);
S_2 = fft(s2,nfft);
figure
plot(f,abs(S_2(1:nfft/2+1)));
soundsc(real(s2),fs)


%duration = 3;   % 3 seconds
[HanMeili,fs]= audioread([pwd,'\meili.wav']);
L = length(HanMeili);
nfft = L;
t = (0:L-1)/fs;
f = fs*(0:nfft/2)/nfft;
Y = fft(HanMeili/nfft,nfft);
Y = 2*Y(1:nfft/2+1);
figure
subplot(2,1,1)
plot(t,HanMeili);
title('The speech of HanMeiLi');
xlabel('Times in [s]');
subplot(2,1,2)
plot(f,abs(Y));
title('The frequency of HanMeiLi');
xlabel('Frequency in [Hz]');
soundsc(real(HanMeili),fs);









