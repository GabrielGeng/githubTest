clear all;close all;
%% Abtain freq information
duration = 3;   % 3 seconds
[temp,fs]= audioread([pwd,'\Timit_Sound\ConcatenatedByREJU_F1.wav']);
frame = duration * fs;
temp = temp(1:frame,1);
[temp_1,fs]= audioread([pwd,'\Timit_Sound\concatenatedByREJU_M7.wav']);
temp_1 = temp_1(1:frame,1);
L = frame;
nfft = L;
t = (0:L-1)/fs;
f = fs*(0:nfft/2)/nfft;

Y = fft(temp/nfft,nfft);
Y = 2*Y(1:nfft/2+1);

Y_1 = fft(temp_1/nfft,nfft);
Y_1 = 2*Y_1(1:nfft/2+1);

figure
subplot(2,2,1)
plot(t,temp);
subplot(2,2,2)
plot(f,abs(Y));
subplot(2,2,3)
plot(t,temp_1);
subplot(2,2,4)
plot(f,abs(Y_1));


%soundsc(temp_1,fs);
%% Array output
c = 340;                                % the sound speed 
nSensors = 4;
senDist = 0.025;
dx = senDist;
dy = 0; 
J = nSensors;
Index = linspace(0,J-1,J);
p = [(-(J-1)/2 + Index.')*senDist,zeros(J,1)];             % size 6*2


direction = [-75,-70];                         % row vector containing angles for which constraints hold
n_source = length(direction);                                       % the number of sources
src_dist = 1;
src_position = src_dist*[cosd(90+direction);sind(90+direction)]; %each column a source
src_sen_dist = NaN(J,n_source);
for srcIdx = 1:n_source
    src_sen_dist(:,srcIdx) = vecnorm(src_position(:,srcIdx).' - p,2,2);
end
sampledelay = fs*src_sen_dist/c;
filterLen = 1024;
h = 1:filterLen;
h_filter = NaN(J,n_source,filterLen);
for senIdx = 1:J
    for srcIdx = 1:n_source
        h_filter(senIdx,srcIdx,:) = sinc(h - sampledelay(senIdx,srcIdx));
    end
end
X = NaN(L,J);
for senIdx = 1:J
    X(:,senIdx) = filter(squeeze(h_filter(senIdx,1,:)),1,temp) + filter(squeeze(h_filter(senIdx,2,:)),1,temp_1);
end

linspec = {'rx','MarkerSize',12,'LineWidth',2};
src_linspec = {'bo','MarkerSize',12,'LineWidth',2};
figure;
plot(p(:,1),p(:,2),linspec{:});  
hold on;
plot(src_position(1,:),src_position(2,:),'bo','MarkerSize',12);
title('Sensor and source positions');
xlabel('x position in meters');
ylabel('y position in meters');
disp('The four microphones are ready !');
xlim([-1.5,1.5]);
ylim([-0.2,1.7]);

% SNR = 40;
% randnoise = randn(L,1);
% for senIdx = 1:J
%     [X(:,senIdx),Noise] = generate_noise(X(:,senIdx),randnoise,SNR);
% end
% soundsc(real(X(:,2)),fs);
save('Observation_wb.mat','X','fs');

%% Array Free-Space Model


%% STFT
%frLen = 2048;
frLen = 1024;
nostft = 2*frLen;
overLapFac = 0.75;
%overLapFac = 0;

noFrame = floor((size(X,1)-frLen)./((1-overLapFac)*frLen));
xframe = zeros(frLen,J,noFrame);
Xstft = zeros(nostft,J,noFrame);
shift = (1-overLapFac)*frLen;
for frIdx =1:noFrame
    if frIdx == 1
        xframe(:,:,frIdx) = X(1:frLen,:);%gengjh extract the first frame,index 1 to 128
        st = 1;
        ed = frLen;
    else
        st = (frIdx-1)*shift+1;
        ed = st+frLen-1;%gengjh ed - st = 127
        xframe(:,:,frIdx) = X(st:ed,:);%gengjh extract every 128 points
    end
    win = hamming(frLen);
    for senIdx = 1:J
        xframe(:,senIdx,frIdx)= xframe(:,senIdx,frIdx).*win;  % windowing %gengjh we must make sure the dimension between xframe and hamming is the same
    end
    Xstft(:,:,frIdx) = fft(xframe(:,:,frIdx),nostft);    
end
%% DOA
noFreq = nostft/2-1;
f_c = fs*(1:nostft/2)/nostft;

Fre_X = Xstft(1:noFreq,:,:);
P_music = 0;
% P_pm = 0;
theta = -90:1:90;
Index = linspace(0,J-1,J);
v = [sin(theta*pi/180);-cos(theta*pi/180)];       % size 2*721
p = [(-(J-1)/2 + Index.')*senDist,zeros(J,1)]; 

for freqIdx = 1:noFreq
R_x = squeeze(Fre_X(freqIdx,:,:))*squeeze(Fre_X(freqIdx,:,:))'./noFrame; 

      
a_theta = exp(-1j*2*pi*f_c(freqIdx)*(p*v)./c);             % steer vector(match weights) 6*721

[U ,eigval] = eig(R_x);                             % the columns of U is the corresponding eigenvectors                                                    % and eigval is a diagonal matrix of eigenvalues
[eigval,index]  = sort(diag(eigval),1,'descend');   % eigval size is 6*1, the index in descend order
U   = U(:,index);                                   % descend the eigenvectors
Us  = U(:,1:n_source);                              % signal subspace
Un  = U(:,n_source+1:J);                            % noise subspace 6*4

P_music = P_music + diag(a_theta'*(Un*Un')*a_theta);
% P_pm = P_pm + diag(a_theta'*R_x*a_theta)./J; 
end
P_sm = noFrame./P_music;          % pseudo power 721*1
figure
linspec = {'b-','LineWidth',2};
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');
xlim([-90,90]);
           
% figure
% linspec = {'b-','LineWidth',2};
% plot(theta, 10*log10(abs(P_pm)), linspec{:});
% xlim([-90,90]);

%% Find the local and global maximum
P_middle = abs(P_sm(2:end-1));
P_front = abs(P_sm(1:end-2));
P_back = abs(P_sm(3:end));
logic_front = (P_middle - P_front)>0;
logic_back = (P_middle - P_back)>0;
logic = logic_front & logic_back;
P_middle(~logic) = min(P_middle);
P_local = [abs(P_sm(1));P_middle;abs(P_sm(end))];
[~,doa_Idx] = maxk(P_local,n_source);
doa = theta(doa_Idx);
[~,minIdx] = min(abs(doa));
doa_source = doa(minIdx);
[~,maxIdx] = max(abs(doa));
interfer = doa(maxIdx);
disp(['The desired source DOA with MUSIC is: ',num2str(doa_source),' deg']);
disp(['The interfering DOA with MUSIC is: ',num2str(interfer),' deg']);
%% Noise
 function[corruptsignal,noise] = generate_noise(clearsignal,randnoise,SNR)
signalpower = (clearsignal'*clearsignal)./length(clearsignal);
randnoisepower = (randnoise'*randnoise)./length(randnoise);
noisepower = signalpower./(10^(0.1*SNR));
noise = sqrt(noisepower./randnoisepower)*randnoise;
corruptsignal = noise + clearsignal;
end
 