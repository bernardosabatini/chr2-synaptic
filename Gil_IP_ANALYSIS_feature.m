
M=figure;
hold on
C=figure;
hold on
L=figure;
hold on

avgTrace.all=[];
avgTrace.M=[];
avgTrace.C=[];
avgTrace.L=[];


otherFields={'ML', 'DV', 'restMode', 'restMean', 'restMedian', 'restSD', ...
	'pulseI', 'pulseV', 'pulseRm', 'sagV', 'checkPulseRpeakMean', 'checkPulseRendMean'...
	};

apResults=[];
apResults.all=[];
apResults.M=[];
apResults.C=[];
apResults.L=[];

allCounter=0;
for counter=1:length(csAllCells)
	newCell=csAllCells(counter);
	disp([newCell.mouseID ' ' num2str(newCell.cellID)])
	aps=find(isnumeric(newCell.nAP) & (newCell.nAP>0) ...
		& (newCell.restMean<-50) & (newCell.restSD<5)...
		& (newCell.checkPulseRpeak>100));
	
	if ~isempty(aps) && ~isempty(newCell.Injection)
		aps=aps(1);
		allFields=fieldnames(newCell.pulseAP{aps});
		
		zone=getZone(newCell);
		
		if isempty(avgTrace.all)
			avgTrace.all=newCell.acq{aps}.data;
		else
			avgTrace.all=...
				avgTrace.all+newCell.acq{aps}.data;
		end
		
		if isempty(avgTrace.(zone))
			avgTrace.(zone)=newCell.acq{aps}.data;
		else
			avgTrace.(zone)=...
				avgTrace.(zone)+newCell.acq{aps}.data;
		end
		
		eval(['figure(' zone ')']);
		plot(newCell.acq{aps}.data);
		
		for fc=1:length(allFields)
			fns=allFields{fc};
			value=newCell.pulseAP{aps}.(fns)(1);
			if isfield(apResults.all, fns)
				apResults.all.(fns)(end+1)=value;
			else
				apResults.all.(fns)=value;
			end
			
			if isfield(apResults.(zone), fns)
				apResults.(zone).(fns)(end+1)=value;
			else
				apResults.(zone).(fns)=value;
			end
		end
		
		for fc=1:length(otherFields)
			fns=otherFields{fc};
			value=newCell.(fns);
			if length(value)>1
				value=value(aps);
			end
			if isfield(apResults.all, fns)
				apResults.all.(fns)(end+1)=value;
			else
				apResults.all.(fns)=value;
			end
			
			if isfield(apResults.(zone), fns)
				apResults.(zone).(fns)(end+1)=value;
			else
				apResults.(zone).(fns)=value;
			end
		end
	end
end

for c={'all', 'M', 'C', 'L'}
	cc=c{1};
	%	length(apResults.(cc).nAP)
	avgTrace.(cc)=...
		avgTrace.(cc) / ...
		length(apResults.(cc).nAP);
end

