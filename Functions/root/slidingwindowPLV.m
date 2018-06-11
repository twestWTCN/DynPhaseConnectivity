function [PLV PLV_tvec] = slidingwindowPLV(period,phi,overlap,type)
[slide_dphi_12,sind] = slideWindow(diff(phi,1,2), floor(period), floor(period*overlap));
slide_dphi_12 = wrapToPi(slide_dphi_12);

if isequal(type,'WPLV')
    win = hanning(size(slide_dphi_12,1));
    for i = 1:size(slide_dphi_12,2)
        PLV(i) = abs(wmean(exp(1i*slide_dphi_12(:,i)),win,1));
    end
else
    PLV = abs(mean(exp(1i*slide_dphi_12),1));
    disp('USING PLV!')
end
% PLV_tvec = (round(median(sind,1)));
% PLV = movmean(PLV,8);

% PLVi = abs(sum(exp(1i*slide_dphi_12),1));
% PLV_tvec = (round(median(sind,1)));

% PLI = abs(mean(sign(slide_dphi_12),1));  % PLI
% PLV= abs(mean(abs(slide_dphi_12).*sign(slide_dphi_12),1))./ mean(abs(slide_dphi_12),1);  %wPLI
PLV_tvec = (round(median(sind,1)));