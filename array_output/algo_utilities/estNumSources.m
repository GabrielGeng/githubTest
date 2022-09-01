function [srcNum_est] = estNumSources(azi, ele)

doas = radtodeg([azi, ele]);
doas(:,1) = wrapTo360(doas(:,1));

BinCenter{1} = 0:5:360;
BinCenter{2} = -90:5:90;
Count = hist3(doas, BinCenter);
lkDir_azi = BinCenter{1};
lkDir_ele = BinCenter{2};

h = fspecial('gaussian', [9, 9], 1.5);
SmoothedCount = filter2(h, Count);
SmoothedCount = filter2(h, SmoothedCount);

spec = SmoothedCount.';
spec = spec./max(max(spec));


% for plotting
if 0
    azimuth    = [110; 170];      % azimuth angle [degree]
    elevation  = [-10; 15];       % elevation angle [degree]
    figure(7); clf; hold on;
    myhist3_normalized(doas, BinCenter);
    xlim([0, 360]); ylim([-90, 90]);
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    for srcIdx = 1:length(azimuth)
        plot3(azimuth(srcIdx),elevation(srcIdx),1.2,'ro','markersize',6,'markerfacecolor','w');
    end
    
    figure(8); clf; hold on;
    surf(lkDir_azi, lkDir_ele, spec);
    xlim([0, 360]); ylim([-90, 90]);
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    for srcIdx = 1:length(azimuth)
        plot3(azimuth(srcIdx),elevation(srcIdx),1.2,'ro','markersize',6,'markerfacecolor','w');
    end
end

% estimate No. of peaks
[ind_ele,ind_azi] = find_2d_peaks(spec, 4);

peaks = NaN(1, 4);
thr = 0.2;
peaks(1) = spec(ind_ele(1), ind_azi(1));
srcNum_est = 1;
for ii = 2:length(ind_ele)
    peaks(ii) = spec(ind_ele(ii), ind_azi(ii));
    if peaks(ii) >= thr*peaks(1)
        srcNum_est = srcNum_est + 1;
    end
end
% azi_est = lkDir_azi(ind_azi);
% ele_est = lkDir_ele(ind_ele);

end

