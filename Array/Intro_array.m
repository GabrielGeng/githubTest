%% evaluate beampattern
J =10;
angles = -180:1:180;
theta = (angles*pi/180).';
my_pattern = 1/J^2 * abs(sin(sin(theta)*J*pi/2)./sin(sin(theta)*pi/2)).^2;  % size 361*1
my_pattern(isnan(my_pattern))=1;
%B = 10*log10(my_pattern);
B = my_pattern;
linspec = {'b-','LineWidth',2};
figure(1)
plot(angles, B, linspec{:});
xlim([-180 180]);
xlabel('Angle in [degrees]');
ylabel('Beampattern');
figure(2)
polarplot(theta,B,linspec{:});
rlim([0 1]);
%% steering (Delay and Sum)
J =10;
theta_0 = 120*pi/180;
angles = -180:1:180;
theta = (angles*pi/180).';
sen_Idx = linspace(0,J-1,J).';
my_pattern = 1/J^2 * abs(sum(exp(-1j*pi*sen_Idx*(sin(theta)-sin(theta_0)).'))).^2;  % size 361*1
my_pattern(isnan(my_pattern))=1;
%B = 10*log10(my_pattern);
B = my_pattern;
linspec = {'b-','LineWidth',2};
figure(3)
plot(angles, B, linspec{:});
xlim([-180 180]);
xlabel('Angle in [degrees]');
ylabel('Beampattern');
figure(4)
polarplot(theta,B,linspec{:});
rlim([0 1]);