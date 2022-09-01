function plotAVSwave( AuData,srcSig,st,ed,fs,fig_handle )

srcNum = size(srcSig,2);

if nargin == 6
    figure(fig_handle);
    clf;
else figure;
end

for srcIdx = 1:srcNum
    subplot(srcNum+1,1,srcIdx);
    plotframeloc(srcSig(:,srcIdx),st,ed,fs); title(['Source signal ',num2str(srcIdx)]);
    ylim([-1, 1])
end
subplot(srcNum+1,1,srcNum+1);
plotframeloc(AuData(:,1,1),st,ed,fs); title('Omnidirectional mic received signal')


end


function plotframeloc(signal,st,ed,fs)

t = (1:length(signal))/fs;

plot(t,signal);


hold on;
v = axis;
plot(st/fs*ones(10,1),linspace(v(3),v(4),10),'r')
plot(ed/fs*ones(10,1),linspace(v(3),v(4),10),'r')
xlabel('Time (s)');
ylabel('Magnitude');

end