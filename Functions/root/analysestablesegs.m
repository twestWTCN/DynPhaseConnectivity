function [phi_dist amp_dist seg_ddt segL_ddt consecSegs H] = analysestablesegs(qstable,tseries,refseries,period,mwid,fsamp,SNR_eps,minsep,zerodel,Ampeps,snrw)
% if isempty(minsep)
%     minsep = 0;
% end

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
    win = tukeywin(size(consecSegs{seglist(j)},2),0.25);
    for sig = 1:4
        phi_dist(j,1) = circ_mean(wrapToPi(tseries(consecSegs{seglist(j)})),win);
        wamp = amp(consecSegs{seglist(j)},sig); %.*win;
        amp_dist(sig,j,1) = median(wamp); %( (median(wamp)-median(amp(:,sig)))/median(amp(:,sig)) )*100;
    end
end
a= 1;
% phi_dist(isnan(phi_dist)) = [];
% amp_dist(isnan(amp_dist)) = [];
