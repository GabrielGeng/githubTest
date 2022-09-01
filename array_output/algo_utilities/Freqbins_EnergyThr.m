function [ pressureFreqbinsEngThr ] = Freqbins_EnergyThr( Noise,frLen,nfft,ovrLapFac, scaling )

%% Collecting frames
noFreq = nfft/2-1;
noFr = 50;
noiseframe = NaN(frLen,noFr,4);
noisefrfft = NaN(nfft,noFr,4);

for frIdx = 1:noFr
    noiseframe(:,frIdx,:) = chopsignal(frIdx,Noise,frLen,ovrLapFac,noFr);
    win = hamming(frLen);
    for senIdx = 1:4
        noiseframe(:,frIdx,senIdx)= noiseframe(:,frIdx,senIdx).*win;  % windowing
    end
    noisefrfft(:,frIdx,:) = fft(noiseframe(:,frIdx,:),nfft);
end

pressureFreqbinsEngThr = zeros(noFreq,1);
for fftIdx = 1:noFreq
    pressureFreqbinsEngThr(fftIdx) = scaling*median( norm(noisefrfft(fftIdx,:,1)) );
end

end

