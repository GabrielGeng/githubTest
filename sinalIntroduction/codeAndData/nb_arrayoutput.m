clear all;close all;
%% generate narrowband signal
n = (1:48000).';
fs = 16000;
L = length(n);
temp = 0.999.^n.*ones(48000,1);
temp = temp.*exp(1j*2*pi*680/fs*n);
temp = exp(1j*2*pi*680/fs*n);

nfft = length(temp);
f = fs*(0:nfft/2)/nfft;
S_1 = fft(temp/nfft,nfft);
figure
plot(f,abs(S_1(1:nfft/2+1)));
%soundsc(real(s1),fs)
temp_1 = zeros(48000,1);
temp_1(1:1000) = 1;
temp_1 = temp_1.*exp(1j*2*pi*680/fs*n);
S_2 = fft(temp_1/nfft,nfft);
figure
plot(f,abs(S_2(1:nfft/2+1)));
%soundsc(real(s2),fs)

%% Array output
nSensors = 4;
senDist = 0.25;
dx = senDist;
dy = 0; 
J = nSensors;
c = 340;                                % the sound speed 
direction = [-6,52];                         % row vector containing angles for which constraints hold
nb_f = 680;

senIdx = linspace(0,J-1,J);
p = [(-(J-1)/2 + senIdx.')*senDist,zeros(J,1)]; 
v = [sin(direction*pi/180);-cos(direction*pi/180)];  % size 2*361

linspec = {'rx','MarkerSize',12,'LineWidth',2};

figure
plot(p(:,1),p(:,2),linspec{:});  
title('Sensor positions');
xlabel('x position in meters');
ylabel('y position in meters');
disp('The four microphones are ready !');


steervector = exp(-1j*2*pi*nb_f*(p*v)./c);   % size J*361
SNR = 40;
X = zeros(L,J);
randnoise = randn(L,1);
for senIdx = 1:J
    X(:,senIdx)= temp.*steervector(senIdx,1)+temp_1.*steervector(senIdx,2);
    [X(:,senIdx),Noise] = generate_noise(X(:,senIdx),randnoise,SNR);
end

%soundsc(real(X(1,:)),fs);
save('Observation_nb.mat','X','fs');
%% DOA
n_source = 2;                                       % the number of sources
X = X.';
[J,Frame] = size(X);
f_c = 680;
R_x = X*X'./Frame; 
theta = -90:0.5:90;
Index = linspace(0,J-1,J);
v = [sin(theta*pi/180);-cos(theta*pi/180)];       % size 2*721
p = [(-(J-1)/2 + Index.')*senDist,zeros(J,1)];             % size 6*2
a_theta = exp(-1j*2*pi*f_c*(p*v)./c);             % steer vector(match weights) 6*721

[U ,eigval] = eig(R_x);                             % the columns of U is the corresponding eigenvectors                                                    % and eigval is a diagonal matrix of eigenvalues
[eigval,index]  = sort(diag(eigval),1,'descend');   % eigval size is 6*1, the index in descend order
U   = U(:,index);                                   % descend the eigenvectors
Us  = U(:,1:n_source);                              % signal subspace
Un  = U(:,n_source+1:J);                            % noise subspace 6*4
P_sm = 1./diag(a_theta'*(Un*Un')*a_theta);          % pseudo power 721*1
figure
linspec = {'b-','LineWidth',2};
plot(theta, 10*log10(abs(P_sm)), linspec{:});
title('MUSIC pseudo power spectrum')
xlabel('Angle in [degrees]');
ylabel('Power spectrum in [dB]');
xlim([-90,90]);

%% Find the Global maximum and visualization
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
comp_vector = steervector(:,1)';
compensate = X.*(comp_vector.');
Block_vector = [1,-1,0,0];
align_out= Block_vector*compensate;
%soundsc(real(align_out(1,:)),fs);
%% NOISE
 function[corruptsignal,noise] = generate_noise(clearsignal,randnoise,SNR)
signalpower = (clearsignal'*clearsignal)./length(clearsignal);
randnoisepower = (randnoise'*randnoise)./length(randnoise);
noisepower = signalpower./(10^(0.1*SNR));
noise = sqrt(noisepower./randnoisepower)*randnoise;
corruptsignal = noise + clearsignal;
end





