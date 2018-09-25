% SPfiles=dir('SP*.mat');
% 
% allCells=cell(1,length(SPfiles));
% for counter=1:length(SPfiles)
% 	cc=SPfiles(counter);
% 	load(cc.name)
% 	allCells{counter}=newCell;
% end
% 
% apFields={'AP_AHP_V',...
%        'AP_thresh_V',...
%     'AP_thresh_time',...
%            'AP_HW_V',...
%              'AP_HW',...
%        'AP_max_dVdT'};
%    
% apFieldsN=length(apFields);
	   
bExact=[-100 -50 0 20 40 50 60 80 100 150 200 250];
bLow=[-100 0 19 49 100 150 200 250 ];
bHigh=[-100 0 31 61 100 150 200 250 ];
pCount=length(bExact);

avgV=zeros(4,pCount);
avgI=zeros(4,pCount);
avgA=zeros(4,pCount);
avgN=zeros(4,pCount);

stdV=zeros(4,pCount);
stdI=zeros(4,pCount);
stdA=zeros(4,pCount);
stdN=zeros(4,pCount);

allCounter=0;
for runThru=1:2 % first run calcualtes the averages.  Second does the Standard Deviation
	for counter=1:length(csAllCells)
		newCell=csAllCells(counter);
		disp([newCell.mouseID ' ' newCell.cellID])
		aps=newCell.nAP;
		aps(isnan(aps))=0;

		if sum(aps)>10
			for bCounter=1:pCount
				ff1=[];
				ff=find((newCell.pulseI==bExact(bCounter)));
				if ~isempty(ff)
					ff1=ff(1);
	% 			else
	% 				ff=find((newCell.pulseI>bLow(bCounter)) & (newCell.pulseI<bHigh(bCounter)));
	% 				if ~isempty(ff)
	% 					ff1=ff(1);
	% 				end
				end
				if ~isempty(ff1) && ~isnan(newCell.pulseV(ff1)) && ...
					(newCell.restMean(ff1)<-50) && (newCell.restMean(ff1)>-80) && (newCell.restSD(ff1)<5)...
					&& (newCell.checkPulseRpeak(ff1)>100)

					zone=getZone(newCell);
					if zone=='M'
						indZone=2;
					elseif zone=='C'
						indZone=3;
					elseif zone=='L'
						indZone=4;
					end

					for ind=[1 indZone]  % run through it once to put in the ALL pile and then in zone specific pile
						if runThru==1
							avgV(ind, bCounter)=avgV(ind, bCounter)+newCell.pulseV(ff1);
							avgI(ind, bCounter)=avgI(ind, bCounter)+newCell.pulseI(ff1);
							avgA(ind, bCounter)=avgA(ind, bCounter)+newCell.nAP(ff1);
							avgN(ind, bCounter)=avgN(ind, bCounter)+1;
						else
							stdV(ind, bCounter)=stdV(ind, bCounter)+(newCell.pulseV(ff1)-avgV(ind, bCounter))^2;
							stdI(ind, bCounter)=stdI(ind, bCounter)+(newCell.pulseI(ff1)-avgI(ind, bCounter))^2;
							stdA(ind, bCounter)=stdA(ind, bCounter)+(newCell.nAP(ff1)-avgA(ind, bCounter))^2;
							stdN(ind, bCounter)=stdN(ind, bCounter)+1;
						end
					end
				end
			end
		end
	end
	if runThru==1
		avgI=avgI./avgN;
		avgV=avgV./avgN;
		avgA=avgA./avgN;
	else
		stdI=(stdI).^(0.5)./stdN;
		stdV=(stdV).^(0.5)./stdN;
		stdA=(stdA).^(0.5)./stdN;
	end
end

figure;
for c=2:4
	vv=avgV(c,:);
	ii=avgI(c,:);
	vv(isnan(vv))=[];
	ii(isnan(ii))=[];
	plot(ii, vv)
	hold on
	vvs=stdV(c,:);
	iis=stdI(c,:);
	vvs(isnan(vvs))=[];
	iis(isnan(iis))=[];
	errorbar(ii, vv, vvs, vvs)
end

figure;
for c=2:4
	aa=avgA(c,:);
	ii=avgI(c,:);
	aa(isnan(aa))=[];
	ii(isnan(ii))=[];
	
	plot(ii, aa)
	hold on
	aas=stdA(c,:);
	iis=stdI(c,:);
	aas(isnan(aas))=[];
	iis(isnan(iis))=[];
	errorbar(ii, aa, aas, aas)	
end

figure;
for c=2:4
	aa=avgA(c,:);
	vv=avgV(c,:);
	aa(isnan(aa))=[];
	vv(isnan(vv))=[];
	
	plot(vv, aa)
	hold on
	aas=stdA(c,:);
	vvs=stdV(c,:);
	aas(isnan(aas))=[];
	vvs(isnan(vvs))=[];
	errorbar(vv, aa, aas, aas)	
end