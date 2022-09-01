function [ RMSAE, AZI, ELE, AE ] = calculateRMSAE( azimuth, elevation, AZI, ELE )
% tentative code;

AZI = AZI(:, ~isnan(AZI(1,:)));
ELE = ELE(:, ~isnan(ELE(1,:)));

srcNum      = length(azimuth);
block_No    = size(AZI,2);

%% Sorting
if length(azimuth)>=2
    if azimuth(1)<=azimuth(2)
        AZI = sort(AZI,1,'ascend');
    elseif azimuth(1)>azimuth(2)
        AZI = sort(AZI,1,'descend');
    end
    if elevation(1)<=elevation(2)
        ELE = sort(ELE,1,'ascend');
    elseif elevation(1)>elevation(2)
        ELE = sort(ELE,1,'descend');
    end
end
%% Error
AE = zeros(block_No,1);
for blockIndex = 1:block_No
    AZI_blk = AZI(:,blockIndex);
    ELE_blk = ELE(:,blockIndex);
    srcNum_est = sum(~isnan(AZI_blk));
    
    if srcNum==2 && srcNum_est==1
        AZI_blk(2) = AZI_blk(1);
        ELE_blk(2) = ELE_blk(1);
        AE_blk = angular_error(azimuth,elevation,AZI_blk,ELE_blk);
    elseif srcNum==3 && srcNum_est==1
        AZI_blk(2) = AZI_blk(1);
        ELE_blk(2) = ELE_blk(1);
        AZI_blk(3) = AZI_blk(1);
        ELE_blk(3) = ELE_blk(1);
        AE_blk = angular_error(azimuth,elevation,AZI_blk,ELE_blk);
    elseif srcNum==3 && srcNum_est==2
        AZI_tmp([1 3]) = AZI_blk(1:2);
        ELE_tmp([1 3]) = ELE_blk(1:2);
        AZI_tmp(2) = AZI_blk(1);
        ELE_tmp(2) = ELE_blk(1);
        AE_blk1 = angular_error(azimuth,elevation,AZI_tmp,ELE_tmp);
        AZI_tmp(2) = AZI_blk(2);
        ELE_tmp(2) = ELE_blk(2);
        AE_blk2 = angular_error(azimuth,elevation,AZI_tmp,ELE_tmp);
        if AE_blk2 > AE_blk1
            AZI_tmp(2) = AZI_blk(1);
            ELE_tmp(2) = ELE_blk(1);
        end
        AZI_blk = AZI_tmp;
        ELE_blk = ELE_tmp;
        AE_blk = min([AE_blk1, AE_blk2]);
    else
        AE_blk = angular_error(azimuth,elevation,AZI_blk,ELE_blk);
    end
    AZI(:,blockIndex) = AZI_blk;
    ELE(:,blockIndex) = ELE_blk;
    AE(blockIndex)    = AE_blk;
end

RMSAE = sqrt( sum(AE.^2)/block_No );

end

function [AE] = angular_error(azimuth,elevation,Azi_est,Ele_est)
srcNum  = length(azimuth);
Error   = NaN(srcNum,1);
for srcIdx = 1:srcNum
    u_real = [cosd(azimuth(srcIdx))*cosd(elevation(srcIdx))  sind(azimuth(srcIdx)).*cosd(elevation(srcIdx))  sind(elevation(srcIdx))];
    u_est  = [cosd(Azi_est(srcIdx))*cosd(Ele_est(srcIdx))  sind(Azi_est(srcIdx)).*cosd(Ele_est(srcIdx))  sind(Ele_est(srcIdx))];
    dif_u = norm(u_real - u_est);
    Error(srcIdx) = 2*asind(dif_u/2);
end
AE = mean(Error);
end

