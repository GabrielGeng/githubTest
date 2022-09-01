clear all;close all;
%% Load data
[temp1,fs]= audioread([pwd,'\audio\-20_70-01.wav']);
temp2 = audioread([pwd,'\audio\-20_70-02.wav']);
temp3 = audioread([pwd,'\audio\-20_70-03.wav']);
temp4 = audioread([pwd,'\audio\-20_70-04.wav']);

nfft = length(temp1);
f = fs*(0:round(nfft/2))/nfft;
S_1 = fft(temp1/nfft,nfft);
figure
plot(f,abs(S_1(1:round(nfft/2)+1)));

soundsc(temp1(:,1),fs);
% A = load('Observation_nb.mat');
% X = A.X;
% fs =A.fs;
X = [temp1,temp2,temp3,temp4];
nSensors = size(X,2);
J = nSensors;
senDist = 0.025;
c = 340;
n_source = 2;
Index = linspace(0,J-1,J);
p = [(-(J-1)/2 + Index.')*senDist,zeros(J,1)];             % size 6*2
linspec = {'rx','MarkerSize',12,'LineWidth',2};
figure
plot(p(:,1),p(:,2),linspec{:});  
title('Sensor positions');
xlabel('x position in meters');
ylabel('y position in meters');
disp('The four microphones are ready !');
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
% theta = -180:1:180;
v = [sin(theta*pi/180);-cos(theta*pi/180)];       % size 2*721
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






