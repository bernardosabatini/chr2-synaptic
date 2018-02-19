function [rPeak, rEnd, tau ] = csCurrentClampPulseAnalysis(dData, acqRate, pulseSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	rPeak=0;
	rEnd=0;
	tau=0;
	
	nPnts=length(dData);
	xx=1/acqRate*[0:(nPnts-1)];
	
	if pulseSize<0
		[vPeak, iPeak]=min(dData);
	else
		[vPeak, iPeak]=max(dData);
	end
	
	vEnd=mean(dData(end-round(nPnts/5):end));
	
	try
		rPeak=1000*vPeak/pulseSize;
		rEnd=1000*vEnd/pulseSize;
		
		ff=fit(xx(1:iPeak)', dData(1:iPeak)'-vPeak, 'exp1', 'StartPoint', [vPeak, -10] );
%		figure; plot(ff, xx, dData-vPeak)
 		tau=-1/ff.b;
	catch
		disp('******** CC PULSE FIT did not work');
	end	
	
end

