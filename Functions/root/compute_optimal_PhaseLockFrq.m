function [maxfrq maxPLV] = compute_optimal_PhaseLockFrq(R,Xdata,band)
frqlist = R.PA.frqrange{band};
for frqn = 1:numel(frqlist)
    [amp phi dphi_12 , ~, ~] = comp_instant_angle_phase(Xdata,frqlist(frqn),R.PA.bwid,band);
    %Epoched
    WinSize = R.PA.slidingwindow*Xdata.fsample;
    [dphi_12,sind] = slideWindow(dphi_12, floor(WinSize*5), 0);
    dphi_12 = wrapToPi(dphi_12);
    
    win = hanning(size(dphi_12,1));
    for i = 1:size(dphi_12,2)
        WPLV(i) = abs(wmean(exp(1i*dphi_12(:,i)),win,1));
    end    
    PLV(frqn) = mean(abs(mean(exp(1i*dphi_12),1)));
    PLI(frqn) = mean(abs(mean(sign(dphi_12),1)));
end
crit = eval(R.PA.optimalPLFrqMeth);
[~,loc] = findpeaks(crit);
if isempty(loc); loc = find(crit==max(crit)); end
[~,i] = min(abs(loc - median(1:numel(frqlist)))); % Find the peak closest to the centre of the band
maxfrq = frqlist(i);
maxPLV = crit(i);
