function[AuData,DatafileName] = avsReceivedSignal(azimuth,elevation,T60,SNR)
%% Room setup parameters
roomDim  = [8 6 4];             % Room dimensions [x y z] (m)
micPos   = [4 3 1.5];           % Receiver position [x y z] (m)
dist       = [1.5;  1.5;  1.5];      % distance of source from sensors [m]
c      = 345;                    % Sound velocity (m/s)
fs     = 16000;                  % Sample frequency (samples/s)
frLen  = 1024;                   % Frame length
nfft   = frLen;

srcNum     = length(azimuth);
dist       = dist(1:srcNum);
azimuth    = azimuth(1:srcNum);
elevation  = elevation(1:srcNum);

u_real     = [ones(srcNum,1) cosd(azimuth).*cosd(elevation) sind(azimuth).*cosd(elevation) sind(elevation)];  % true doa vector
% doa_dif_1 = acosd(u_real(1,2:4)*u_real(2,2:4)')
% doa_dif_2 = acosd(u_real(2,2:4)*u_real(3,2:4)')
senNum = length(u_real);

srcPos = repmat(micPos,srcNum,1) + [dist.*cosd(azimuth).*cosd(elevation), dist.*sind(azimuth).*cosd(elevation), (dist.*sind(elevation))];
visualroom(srcPos, micPos, roomDim);

duration = 19;
signal_size  = fs*duration;

srcFiles = cell(3,2);
srcFiles(1, :) = {'ConcatenatedByREJU_M1.wav', 'ConcatenatedByREJU_M2.wav'};
%srcFiles(2, :) = {'ConcatenatedByREJU_M1.wav', 'ConcatenatedByREJU_M2.wav'};
srcFiles(2, :) = {'ConcatenatedByREJU_F1.wav', 'ConcatenatedByREJU_F2.wav'};
srcFiles(3, :) = {'ConcatenatedByREJU_M7.wav', 'ConcatenatedByREJU_M8.wav'};

srcSig = [];
for srcIdx = 1:3
    srcsig_col = [];
    for fileIdx = 1:2
        srcsigtemp = audioread([pwd,'\Timit_Sound\', srcFiles{srcIdx, fileIdx}],[1 signal_size]);
        srcsig_col = [srcsig_col; srcsigtemp/std(srcsigtemp)/15];
    end
    srcSig = [srcSig srcsig_col];
end

srcSig = srcSig(:,1:srcNum);
signal_size = size(srcFiles, 2)*signal_size;
clear srcsig_col srcsigtemp;

%% Room impulse response
if T60~=0
    order = -1;                 % -1 equals maximum reflection order
    rir_len = 2048; %T60*fs;
    mtype = {'omnidirectional', 'bidirectional', 'bidirectional', 'bidirectional'};    % Type of microphone
    orientation = {[0 0], [0 0], [pi/2, 0], [0, pi/2]};     % Microphone orientation (rad)
    hp_filter = 0;              % Disable high-pass filter
    RIR = zeros(rir_len,srcNum,senNum);
    for srcIdx = 1:srcNum
        for senIdx = 1:senNum
            RIR(:,srcIdx,senIdx) = rir_generator(c, fs, micPos, srcPos(srcIdx,:), roomDim, T60, rir_len, mtype{senIdx}, order, length(roomDim), orientation{senIdx}, hp_filter);
        end
    end
elseif T60==0 % free-space RIR
    rir_len = 1024;
    RIR = zeros(rir_len,srcNum,senNum);
    for srcIdx = 1:srcNum
        dist =  norm(srcPos(srcIdx,:) - micPos);
        sampdelay = round(dist*fs/c);
        STER = [1; cosd(azimuth(srcIdx))*cosd(elevation(srcIdx)); sind(azimuth(srcIdx))*cosd(elevation(srcIdx)); sind(elevation(srcIdx))];
        RIR(sampdelay,srcIdx,:) = STER;
    end
end
clear order rir_len mtype orientation hp_filter;

%% Plot RIR
figure;
for senIdx = 1:senNum
    subplot(senNum,1,senIdx);
    plot( (1:size(RIR,1))*1000/fs , RIR(:,1,senIdx));
end
xlabel('Time (ms)');

%% Obtain direct sound (free-space)
if 1
    order = 0;                 % -1 equals maximum reflection order
    rir_len = 2048; %T60*fs;
    mtype = {'omnidirectional', 'bidirectional', 'bidirectional', 'bidirectional'};    % Type of microphone
    orientation = {[0 0], [0 0], [pi/2, 0], [0, pi/2]};     % Microphone orientation (rad)
    hp_filter = 0;              % Disable high-pass filter
    RIR_direct = zeros(rir_len,srcNum,senNum);
    for srcIdx = 1:srcNum
        for senIdx = 1:senNum
            RIR_direct(:,srcIdx,senIdx) = rir_generator(c, fs, micPos, srcPos(srcIdx,:), roomDim, T60, rir_len, mtype{senIdx}, order, length(roomDim), orientation{senIdx}, hp_filter);
        end
    end
end

figure;
for senIdx = 1:senNum
    subplot(senNum,1,senIdx);
    plot( (1:size(RIR_direct,1))*1000/fs , RIR_direct(:,1,senIdx));
end
xlabel('Time (ms)');


%% ALL reflections
RIR_reflect = RIR - RIR_direct;

figure;
for senIdx = 1:senNum
    subplot(senNum,1,senIdx);
    plot( (1:size(RIR_reflect,1))*1000/fs , RIR_reflect(:,1,senIdx));
end
xlabel('Time (ms)');

%% Mic received signals
AuData = zeros(signal_size,senNum);
Noise  = zeros(signal_size,senNum);
load('Randn_Noise.mat','Randn_Noise');  
    
for senIdx = 1:senNum
    for srcIdx = 1:srcNum
        AuData(:,senIdx) = filter(RIR(:,srcIdx,senIdx),1,srcSig(:,srcIdx)) + AuData(:,senIdx);
    end
    if SNR(senIdx)~=Inf
        noiseType = 'WGN';
%         [AuData(:,senIdx),Noise(:,senIdx)] = gennoise(AuData(:,senIdx), randn(signal_size,1), SNR(senIdx)); % WGN
        [AuData(:,senIdx),Noise(:,senIdx)] = gennoise(AuData(:,senIdx), Randn_Noise(:,senIdx), SNR(senIdx));
        
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
DatafileName = 'AVS_sim_data';

%% Plot waveform
plotAVSwave(AuData,srcSig,1,signal_size,fs);
end

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