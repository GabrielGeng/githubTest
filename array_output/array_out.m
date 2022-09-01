clear all; close all; clc;
addpath(genpath(pwd));

%% Room setup parameters
nSensors = 4;
senDist = 0.05;
roomDim  = [8 6 4];             % Room dimensions [x y z] (m)      
arrayCenter = [roomDim(1:2)/2 1.5];
arrayHeight = arrayCenter(3);
micPoints = generateUlaCoords(arrayCenter,nSensors,senDist,0,arrayHeight);
micPos = micPoints; % Receiver position [x y z] (m)

dist       = [1.5;  1.5;  1.5];      % distance of source from sensors [m]
if ~exist('azimuth','var')&&~exist('elevation','var')
     %azimuth    = [0; 30];      % azimuth angle [degree]
     %elevation  = [0; 0];       % elevation angle [degree]
     azimuth = 40;
     elevation = 0;
end

if ~exist('T60','var')&&~exist('SNR','var')
    T60    = 0;                  % Reverberation time (s)
    SNR    = [20 20*ones(1,nSensors-1)];     % SNR (dB)
end

c      = 340;                    % Sound velocity (m/s)
fs     = 22050;                  % Sample frequency (samples/s)
%frLen  = 1024;                   % Frame length
%nfft   = frLen;

srcNum     = length(azimuth);
dist       = dist(1:srcNum);
azimuth    = azimuth(1:srcNum);
elevation  = elevation(1:srcNum);

srcPos = repmat(arrayCenter,srcNum,1) + [dist.*cosd(azimuth+90).*cosd(elevation), dist.*sind(azimuth+90).*cosd(elevation), (dist.*sind(elevation))];
visualroom(srcPos, micPos, roomDim);

%duration = 19;
%signal_size  = fs*duration;

%srcFiles = cell(3,2);
%srcFiles(1, :) = {'ConcatenatedByREJU_M1.wav', 'ConcatenatedByREJU_M2.wav'};
%srcFiles(2, :) = {'ConcatenatedByREJU_F1.wav', 'ConcatenatedByREJU_F2.wav'};
%srcFiles(3, :) = {'ConcatenatedByREJU_M7.wav', 'ConcatenatedByREJU_M8.wav'};

%srcSig = [];
%for srcIdx = 1:3
%    srcsig_col = [];
%    for fileIdx = 1:2
%        srcsigtemp = audioread([pwd,'\Timit_Sound\', srcFiles{srcIdx, fileIdx}],[1 signal_size]);
%        srcsig_col = [srcsig_col; srcsigtemp/std(srcsigtemp)/15];
%    end
%    srcSig = [srcSig srcsig_col];
%end

%srcSig = srcSig(:,1:srcNum);
%signal_size = size(srcFiles, 2)*signal_size;
clear srcsig_col srcsigtemp;
srcSig = audioread([pwd,'\data_project\9mm.wav']);
srcSig = srcSig(:,1:srcNum);
signal_size = size(srcSig,1);
%srcSig = srcSig(1:4000,:);
%signal_size = size(srcSig,1);
%samp = linspace(1,40000,40000).';
%srcSig = cos(2*pi*(fs/2)*samp);
%signal_size = size(srcSig,1);

%% Room impulse response
if T60~=0
    order = -1;                 % -1 equals maximum reflection order
    rir_len = 2048; %T60*fs;
    mtype = {'omnidirectional', 'omnidirectional', 'omnidirectional','omnidirectional'};    % Type of microphone
    orientation = {[0 0], [0 0], [0 0], [0 0]};     % Microphone orientation (rad)
    hp_filter = 0;              % Disable high-pass filter
    RIR = zeros(rir_len,srcNum,4);
    for srcIdx = 1:srcNum
        for senIdx = 1:4
            RIR(:,srcIdx,senIdx) = rir_generator(c, fs, micPos(senIdx,:), srcPos(srcIdx,:), roomDim, T60, rir_len, mtype{senIdx}, order, length(roomDim), orientation{senIdx}, hp_filter);
        end
    end
%elseif T60==0 % free-space RIR
%    rir_len = 1024;
%    RIR = zeros(rir_len,srcNum,4);
%    for srcIdx = 1:srcNum
%        dist =  norm(srcPos(srcIdx,:) - micPos);
%        sampdelay = round(dist*fs/c);
%        STER = [1; cosd(azimuth(srcIdx))*cosd(elevation(srcIdx)); sind(azimuth(srcIdx))*cosd(elevation(srcIdx)); sind(elevation(srcIdx))];
%        RIR(sampdelay,srcIdx,:) = STER;
%    end
elseif T60==0
    T60=0.2;
    rir_len = 2048;
    mtype = 'omnidirectional';
    order = 0;
    RIR = zeros(rir_len,srcNum,4);
    for srcIdx = 1:srcNum
        for senIdx = 1:4
             RIR(:,srcIdx,senIdx)= rir_generator(c, fs, micPos(senIdx,:),srcPos(srcIdx,:),roomDim,...
                        T60,rir_len,mtype,order,length(roomDim),[],1);
        end
    end
end
clear order rir_len mtype orientation hp_filter;

%% Plot RIR
figure;
for senIdx = 1:4
    subplot(4,1,senIdx);
    plot( (1:size(RIR,1))*1000/fs , RIR(:,1,senIdx));
end
xlabel('Time (ms)');
%
% figure;
% for senIdx = 1:4
%     subplot(4,1,senIdx);
%     spec = fft(RIR(:,1,senIdx));
%     specLen = length(spec);
%     plot(1:specLen/2,20*log10(abs(spec(1:specLen/2))));
% end

%% Mic received signals
AuData = zeros(signal_size,4);
Noise  = zeros(signal_size,4);

% Fs = 48000;
% [Temp, ~] = audioread([pwd,'\MeetingRoomD_28_July_2015\MeetingRoomD_050.wav'],[1, Fs*duration*3]); % noise only
% for senIdx = 1:4
%     RecNoise(:,senIdx) = resample(double(Temp(:,senIdx)),fs,Fs);
% end
% clear Temp
    
for senIdx = 1:4
    for srcIdx = 1:srcNum
        AuData(:,senIdx) = filter(RIR(:,srcIdx,senIdx),1,srcSig(:,srcIdx)) + AuData(:,senIdx);
    end
    if SNR(senIdx)~=Inf
        noiseType = 'WGN';
        [AuData(:,senIdx),Noise(:,senIdx)] = gennoise(AuData(:,senIdx), randn(signal_size,1), SNR(senIdx)); % WGN
        
%         noiseType = 'Recorded';
%         [AuData(:,senIdx),Noise(:,senIdx)] = gennoise(AuData(:,senIdx), RecNoise(:, senIdx), SNR(senIdx));
    end
end

disp(['T60 = ',num2str(T60*1000),' ms']);
disp(['SNR for WGN = ',num2str(SNR),' dB']);

% clear RecNoise

%% Save and Listening
% soundsc(AuData(:,1),fs);
save('AVS_sim_data');

%% Plot waveform
plotAVSwave(AuData,srcSig,1,signal_size,fs);

%% Plot spectrogram
% figure;
% spectrogram(srcSig(:,1),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('source signal 1');
% figure;
% spectrogram(srcSig(:,2),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('source signal 2');
% figure;
% spectrogram(AuData(:,1),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('received signal o');
% figure;
% spectrogram(AuData(:,2),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('received signal x');
% figure;
% spectrogram(AuData(:,3),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('received signal y');
% figure;
% spectrogram(AuData(:,4),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('received signal z');
% figure;
% spectrogram(Noise(:,1),frLen,0.75*frLen,nfft,fs,'yaxis');
% title('noise signal o');
%% DOA Estimate
X = AuData.';
[J,Frame] = size(X);
f_c = 250;
R_x = X*X'./Frame; 
theta = -180:0.5:180;
Index = linspace(0,J-1,J);
v = [sin(theta*pi/180);-cos(theta*pi/180)];       % size 2*721
p = [(-(J-1)/2 + Index.')*senDist,zeros(J,1)];             % size 6*2
a_theta = exp(-1j*2*pi*f_c*(p*v)./c);             % steer vector(match weights) 6*721
P_pm = diag(a_theta'*R_x*a_theta)./J;             % spatial power size 721*1
n_source = 1;                                       % the number of sources
[U ,eigval] = eig(R_x);                             % the columns of U is the corresponding eigenvectors
                                                    % and eigval is a diagonal matrix of eigenvalues
[eigval,index]  = sort(diag(eigval),1,'descend');   % eigval size is 6*1, the index in descend order
U   = U(:,index);                                   % descend the eigenvectors
Us  = U(:,1:n_source);                              % signal subspace
Un  = U(:,n_source+1:J);                            % noise subspace 6*4
P_sm = 1./diag(a_theta'*(Un*Un')*a_theta);          % pseudo power 721*1
linspec = {'b-','LineWidth',2};
figure
plot(theta, 10*log10(abs(P_pm)), linspec{:});
figure
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');