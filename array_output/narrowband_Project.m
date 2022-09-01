clear all;close all;
%% Load data
A = load('Observation_nb.mat');
X = A.X;
fs =A.fs;
%soundsc(real(X(:,1)),fs);
[Frame,nSensors] = size(X);
J = nSensors;
senDist = 0.25;
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
%% DOA
f_c = 680;
Trans_X = X.';
R_x = Trans_X*Trans_X'./Frame; 
theta = -90:0.5:90;
v = [sin(theta*pi/180);-cos(theta*pi/180)];       % size 2*721
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
v = [sin(interfer*pi/180);-cos(interfer*pi/180)]; 
comp_vector = exp(1j*2*pi*f_c*(p*v)./c);
compensate = X.*(comp_vector).';
Block_vector = [1;-1;0;0];
align_out= compensate*Block_vector;
%soundsc(real(align_out(:,1)),fs);

