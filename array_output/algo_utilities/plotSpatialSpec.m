function plotSpatialSpec( Spec,lkDir_azi,lkDir_ele,azimuth,elevation,type,fig_handle)

Spec_max = max(max(Spec)); % normalize

figure(fig_handle);
clf;
switch lower(type)
    case 'log'
        mesh(lkDir_azi,lkDir_ele,log(Spec/Spec_max));
    case 'rel'
        Spec = Spec - min(min(Spec));
        mesh(lkDir_azi,lkDir_ele,Spec/max(max(Spec)));
end

hold on;
v = axis;
view(0,90);
for srcIdx = 1:length(azimuth)
    p1 = plot3(azimuth(srcIdx),elevation(srcIdx),v(6),'ko','markersize',6,'markerfacecolor','k');
end
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');zlabel('Power');
axis([0 360 -90 90]);

end

