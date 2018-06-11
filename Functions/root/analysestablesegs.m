function [phi_dist amp_dist seg_ddt segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,refseries,period,mwid,fsamp,SNR_eps,minsep,zerodel,Ampeps,snrw)
if isempty(minsep)
    minsep = 0;
end

amp = refseries;
consecSegs = SplitVec(qstable,'consecutive');
% consecSegs = fixGap(consecSegs,minsep);
% lengths
segL_ddt = cellfun('length',consecSegs);
seglist = find(segL_ddt>(period));

seg_ddt = [consecSegs{segL_ddt>(period)}];
segL_ddt = segL_ddt(segL_ddt>(period))/fsamp;
% phase angles of segments
phi_dist = [];
amp_dist = [];
H = [];
for j = 1:numel(seglist)
    if sum(mean(snrw(consecSegs{seglist(j)},1:2))>SNR_eps) == 2 % && abs(circ_mean(wrapToPi(tseries(consecSegs{seglist(j)}))))>(pi/12)%%&& ((median(amp(consecSegs{seglist(j)},3))-median(amp(:,3)))/median(amp(:,3))*100)>50
        phi_dist(j) = circ_mean(wrapToPi(tseries(consecSegs{seglist(j)})));
        %         phi_2_dist(j) = circ_mean(wrapTo2Pi(tseries(consecSegs{j},2)));
        %         phi_1_dist(j) = circ_mean(wrapTo2Pi(tseries(consecSegs{j},1)));
        amp_dist(1,j) = (median(amp(consecSegs{seglist(j)},1))-Ampeps(1))/(Ampeps(1))*100;
        amp_dist(2,j) = (median(amp(consecSegs{seglist(j)},2))-Ampeps(2))/(Ampeps(2))*100;
        amp_dist(3,j) = (median(amp(consecSegs{seglist(j)},3))-Ampeps(3))/(Ampeps(3))*100;
        % Amp correlations
        x1 = amp(consecSegs{seglist(j)},1); x2 = amp(consecSegs{seglist(j)},2);
%         [r,i] = my_xcorr(x1,x2,-10:10); %[r,i] = xcorr(x1,x2,fsamp,'unbiased');
r = 1; i =0;
        [r,ii] = max(r);
        %         subplot(1,2,1); plot(1:numel(x1),x1,1:numel(x2),x2);
        %         subplot(1,2,2); plot(x1+x2);
        H(1,j) = r;
        H(2,j) = i(ii)/fsamp;
        %         H(3,j) = mean(amp(consecSegs{seglist(j)},3))/mean(diff(amp(consecSegs{seglist(j)},3)));
        %         H(1,j) = mean(amp(consecSegs{seglist(j)},1))/mean(diff(amp(consecSegs{seglist(j)},1)));
        %         H(2,j) = mean(amp(consecSegs{seglist(j)},2))/mean(diff(amp(consecSegs{seglist(j)},2)));
        %         H(3,j) = mean(amp(consecSegs{seglist(j)},3))/mean(diff(amp(consecSegs{seglist(j)},3)));
        % ####        amp_dist(3,j) = max(amp(consecSegs{j},3))/mean(amp(:,3));
    else
        phi_dist(j) = NaN;
        amp_dist(1,j) = NaN;
        amp_dist(2,j) = NaN;
        amp_dist(3,j) = NaN;
        H(1,j) = NaN;
        H(2,j) =NaN;
        
    end
end
% a= 1;
% phi_dist(isnan(phi_dist)) = [];
% amp_dist(isnan(amp_dist)) = [];
