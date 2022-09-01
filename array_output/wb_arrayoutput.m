clear all;close all;
%% Abtain freq information
duration = 3;   % 3 seconds
%[temp,fs]= audioread([pwd,'\data_project\M16.wav']);
[temp,fs]= audioread([pwd,'\Timit_Sound\ConcatenatedByREJU_F1.wav']);
frame = duration * fs;
temp = temp(1:frame,1);
%[temp_1,fs]= audioread([pwd,'\data_project\9mm.wav']);
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
nSensors = 4;
senDist = 0.25;
dx = senDist;
dy = 0; 
J = nSensors;
my_array = arrays.ULA(J,dx,dy);
disp('The four microphones are ready !');
figure
my_array.plot();  
c = 340;                                % the sound speed 

direction = [-56,41];                         % row vector containing angles for which constraints hold
n_source = length(direction);                                       % the number of sources
%target = [0,1];                         % row vector containing target values for the beampattern

b = beamformer;
set(b, 'array',         my_array);
set(b, 'n_sensor',      J);
set(b, 'sound_speed',   c);

F1 = fs*(1:nfft/2)/nfft;
F2 = flip(F1(1:nfft/2-1));
F= [0,F1,F2];
%F = fs*(-nfft/2:nfft/2-1)/nfft;
Y1 = fft(temp,nfft);
Y2 = fft(temp_1,nfft);
steermatrix = zeros(J,n_source,nfft);
for freqIdx = 1:nfft
    steermatrix(:,:,freqIdx) =  array_response_vector(b,direction, F(freqIdx));
end

Fre1 = zeros(nfft,J);
Fre2 = zeros(nfft,J);
for senIdx = 1:J
    Fre1(:,senIdx) = Y1.*squeeze(steermatrix(senIdx,1,:));
end
X1 = ifft(Fre1);
%soundsc(real(X1(:,2)),fs);
for senIdx = 1:J
    Fre2(:,senIdx) = Y2.*squeeze(steermatrix(senIdx,n_source,:));
end
X2 = ifft(Fre2);
%soundsc(real(X2(:,2)),fs);
X = X1+X2;
SNR = 40;
randnoise = randn(L,1);
for senIdx = 1:J
    [X(:,senIdx),Noise] = generate_noise(X(:,senIdx),randnoise,SNR);
end
%soundsc(real(X(:,2)),fs);
save('Observation_wb.mat','X','fs');


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

%% Null-steering
nfft = size(X,1);
Freq_X = fft(X,nfft);
F1 = fs*(1:nfft/2)/nfft;
F2 = flip(F1(1:nfft/2-1));
F= [0,F1,F2];
v = [sin(interfer*pi/180);-cos(interfer*pi/180)]; 
compensate_matrix = zeros(nfft,J);
for freqIdx = 1:nfft
    compensate_matrix(freqIdx,:) = exp(1j*2*pi*F(freqIdx)*(p*v)./c); 
end
comp_Freq_X = Freq_X.*compensate_matrix;
compensate_X = ifft(comp_Freq_X);
Block_vector = zeros(J,1);
Block_vector(1:2) = [1,-1];
Align_out = compensate_X*Block_vector;


%Align_out= compensate_X(:,1)-compensate_X(:,2);
soundsc(real(Align_out(:,1)),fs);
%% Noise
 function[corruptsignal,noise] = generate_noise(clearsignal,randnoise,SNR)
signalpower = (clearsignal'*clearsignal)./length(clearsignal);
randnoisepower = (randnoise'*randnoise)./length(randnoise);
noisepower = signalpower./(10^(0.1*SNR));
noise = sqrt(noisepower./randnoisepower)*randnoise;
corruptsignal = noise + clearsignal;
end
 