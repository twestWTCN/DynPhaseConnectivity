function [phi_dist amp_dist segL_ddt consecSegs] = analysestablesegs_PLV_vs_Phi_Amp(PLV_stable,RP,Amp,period,fsamp)
% if isempty(minsep)
%     minsep = 0;
% end
consecSegs = SplitVec(PLV_stable,'consecutive');
segL_ddt = cellfun('length',consecSegs);
seglist = find(segL_ddt>(period));

seg_ddt = [consecSegs{segL_ddt>(period)}];
segL_ddt = segL_ddt(segL_ddt>(period))/fsamp;
% phase angles of segments
phi_dist = [];
amp_dist = [];

for j = 1:numel(seglist)
    win = tukeywin(size(consecSegs{seglist(j)},2),0.25);
    for sig = 1:4
        phi_dist(j,1) = circ_mean(wrapToPi(RP(consecSegs{seglist(j)})),win);
        wamp = Amp(consecSegs{seglist(j)},sig); %.*win;
        amp_dist(sig,j,1) = median(wamp); %( (median(wamp)-median(amp(:,sig)))/median(amp(:,sig)) )*100;
    end
end
a= 1;
% phi_dist(isnan(phi_dist)) = [];
% amp_dist(isnan(amp_dist)) = [];
