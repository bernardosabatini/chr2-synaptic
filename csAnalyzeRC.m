function [rm, rs, cm] = csAnalyzeRC(dData, pulseSize, acqRate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	nPnts=length(dData);
	xx=1/acqRate*[0:(nPnts-1)];
	tt=floor(2*nPnts/3);
	bl=mean(dData(tt:end)); % assume last 1/3 is in steady state
	[downPeak, downPeakInd]=min(dData);
	nnPnts=length(dData)-downPeakInd+1;
	
	ff=fit(xx(1:(end-downPeakInd+1))', dData(downPeakInd:end)'-bl, 'exp1', 'StartPoint', [downPeak, -.1] );
%	figure; plot(ff, xx, dData-bl)
	tau=-1/ff.b;
	peak0=min(ff.a*exp(-(downPeakInd-1)*ff.b/acqRate), downPeak);
	rs=1000*pulseSize/peak0;
	rmprs=1000*pulseSize/bl;
	rm=rmprs-rs;
	ppp=rs*rm/(rs+rm);
	cm=tau/ppp;
	
	rmF=-ff.a/pulseSize*1000;
	rmE=bl/pulseSize*1000;
	cm=1000*tau/rmE;
end

