function [phishift phipeak] = findAmpPhi(R,amp,phi,phiBin)
% phiBin = linspace(-pi,pi,N);
phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
if strncmp(R.PA.plotting.realignMeth,'WghtedPrctleAmp',15)
    prc = str2double(R.PA.plotting.realignMeth(end-1:end));
else
    prc = 85;
end
for i = 1:length(phiBin)-1
    binDat = amp(:,phi>=phiBin(i) & phi<=phiBin(i+1))';
    ampBinNumel(i) = numel(binDat);
    ampBinMu(:,i) = nanmean(binDat);
    WampBinMean(:,i) = mean(binDat)*(numel(binDat)/numel(amp)); %%     Weighted mean
    WampBinPrct(:,i) = prctile(binDat,prc)*(numel(binDat)/numel(amp)); %%     Weighted mean
    RPBinMean(:,i) = circ_mean(phi(:,phi>=phiBin(i) & phi<=phiBin(i+1))');
    binN(:,i) = numel(binDat);
    ampBinSEM(:,i) = nanstd(binDat)/size(binDat,1);
end
WLrgampBinMu = ampBinMu.*(ampBinNumel./max(ampBinNumel));
switch R.PA.plotting.realignMeth(isstrprop(R.PA.plotting.realignMeth,'alpha'))
    case 'PhiHist'
        phishift = phiBinMid(binN==max(binN));
        phishift = wrapToPi(phi - phishift(1)); % Use the first maximum
    case 'MaxLBAmp'
        phishift = phiBinMid(ampBinMu==max(ampBinMu));
        phishift = wrapToPi(phi - phishift(1)); % Use the first maximum
    case 'WghtedPrctleAmp'
        phishift = phiBinMid(WampBinPrct==max(WampBinPrct));
        phipeak = RPBinMean(WampBinPrct==max(WampBinPrct));
        if isempty(phishift); phishift = 0; disp('Couldnt Shift Phi!!!'); end
        phishift = wrapToPi(phi - phishift(1)); % Use the first maximum
    case 'WghtedMaxBin'
        phishift = phiBinMid(WLrgampBinMu==max(WLrgampBinMu));
        phishift = wrapToPi(phi - phishift(1)); % Use the first maximum
    case 'WghtedMeanAmp'
        phishift = phiBinMid(WampBinMean==max(WampBinMean));
        phishift = wrapToPi(phi - phishift(1)); % Use the first maximum
    case 'PhiMean'
        phishift = wrapToPi(phi - circ_mean(phi(~isnan(phi))'));
    case 'noshift'
        phishift = phi;
    otherwise 
        disp('No Valid Phase Centring Method Specified')
end
if isempty(phipeak); phipeak = phiBin(3); disp('Couldnt Shift Phi!!!'); end
phipeak = phipeak(1);

% phishift = phiBinMid(ampBinMu==max(ampBinMu));
% % wrapToPi(phiBinMid - phishift)
% phishift = wrapToPi(phi - phishift(1));