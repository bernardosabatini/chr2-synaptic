function [ output_args ] = csRunAnalysis_AMPANMDA_SP153( cellList )
% csRunAnalysis_SecondHalfFig4
%	has the check pulse in a different place than the standard

	csRunAnalysis_Flex(cellList, ...
		'pulseStart', 1000, ...
		'maxRest', 1000, ...
		'minRest', -300, ...
		'maxRestSD', 25, ...
		'minRm', 40, ...
		'maxRs', 50, ...		
		'chargeWindow', 50, ...
		'fakeShift', -100,...
		'checkPulseStart', 2800, ...		
		'checkPulseEnd', 2850, ...		
		'autoPosNegPeak', 1	);
