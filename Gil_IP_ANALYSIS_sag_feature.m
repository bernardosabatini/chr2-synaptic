
M=figure;
hold on
C=figure;
hold on
L=figure;
hold on

avgTrace.sag.all=[];
avgTrace.sag.M=[];
avgTrace.sag.C=[];
avgTrace.sag.L=[];

allCounter=[];
allCounter.all=0;
allCounter.M=0;
allCounter.C=0;
allCounter.L=0;

for counter=1:length(csAllCells)
	newCell=csAllCells(counter);
	disp([newCell.mouseID ' ' num2str(newCell.cellID)])
	aps=find((newCell.pulseI==-100) & isnumeric(newCell.pulseV) & ...
		(newCell.restMean<-50) & (newCell.restMean>-80) & (newCell.restSD<5)...
		& (newCell.checkPulseRpeak>100));
	
	if ~isempty(aps) && ~isempty(newCell.Injection)
		aps=aps(1);
		
		zone=getZone(newCell);
		
		allCounter.all=allCounter.all+1;
		if isempty(avgTrace.sag.all)
			avgTrace.sag.all=newCell.acq{aps}.data;
		else
			avgTrace.sag.all=...
				avgTrace.sag.all+newCell.acq{aps}.data;
		end
		
		allCounter.(zone)=...
			allCounter.(zone)+1;
		if isempty(avgTrace.sag.(zone))
			avgTrace.sag.(zone)=newCell.acq{aps}.data;
		else
			avgTrace.sag.(zone)=...
				avgTrace.sag.(zone)+newCell.acq{aps}.data;
		end
		
		eval(['figure(' zone ')']);
		plot(newCell.acq{aps}.data);
	end
end

for c={'all', 'M', 'C', 'L'}
	cc=c{1};
	%	length(apResults.(cc).nAP)
	avgTrace.sag.(cc)=...
		avgTrace.sag.(cc) / ...
		allCounter.(cc);
end

