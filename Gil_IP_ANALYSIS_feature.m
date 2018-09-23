
M=figure;
hold on
C=figure;
hold on
L=figure;
hold on


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
	aps=find(isnumeric(newCell.nAP) & newCell.nAP>0);
	
	if ~isempty(aps) && ~isempty(newCell.Injection)
		aps=aps(1);
		allFields=fieldnames(newCell.pulseAP{aps});
		
		for fc=1:length(allFields) 
			fns=allFields{fc};
			value=newCell.pulseAP{aps}.(fns)(1);
			if isfield(apResults.all, fns)
				apResults.all.(fns)(end+1)=value;
			else
				apResults.all.(fns)=value;
			end
			
			if isfield(apResults.(newCell.Injection), fns)
				apResults.(newCell.Injection).(fns)(end+1)=value;
			else
				apResults.(newCell.Injection).(fns)=value;
			end
			
			eval(['figure(' newCell.Injection ')']);
			plot(newCell.acq{aps}.data);
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
			
			if isfield(apResults.(newCell.Injection), fns)
				apResults.(newCell.Injection).(fns)(end+1)=value;
			else
				apResults.(newCell.Injection).(fns)=value;
			end
		end
	end
end

