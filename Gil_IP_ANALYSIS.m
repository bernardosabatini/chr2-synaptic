% SPfiles=dir('SP*.mat');
% 
% allCells=cell(1,length(SPfiles));
% for counter=1:length(SPfiles)
% 	cc=SPfiles(counter);
% 	load(cc.name)
% 	allCells{counter}=newCell;
% end

apFields={'AP_AHP_V',...
       'AP_thresh_V',...
    'AP_thresh_time',...
           'AP_HW_V',...
             'AP_HW',...
       'AP_max_dVdT'};
   
apFieldsN=length(apFields);
	   
bExact=[-100 -50 0 20 40 50 60 80 100 150 200 250];
bLow=[-100 0 19 49 100 150 200 250 ];
bHigh=[-100 0 31 61 100 150 200 250 ];
pCount=length(bExact);

avgV=zeros(4,pCount);
avgI=zeros(4,pCount);
avgA=zeros(4,pCount);
avgN=zeros(4,pCount);

allCounter=0;
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
				(newCell.restMean(ff1)<-50) && (newCell.restSD(ff1)<5)...
				&& (newCell.checkPulseRpeak(ff1)>100)
			
				avgV(1, bCounter)=avgV(1, bCounter)+newCell.pulseV(ff1);
				avgI(1, bCounter)=avgI(1, bCounter)+newCell.pulseI(ff1);
				avgA(1, bCounter)=avgA(1, bCounter)+newCell.nAP(ff1);
				avgN(1, bCounter)=avgN(1, bCounter)+1;

				zone=getZone(newCell);
				
				if zone=='M'
					ind=2;
				elseif zone=='C'
					ind=3;
				elseif zone=='L'
					ind=4;
				end
				
				avgV(ind, bCounter)=avgV(ind, bCounter)+newCell.pulseV(ff1);
				avgI(ind, bCounter)=avgI(ind, bCounter)+newCell.pulseI(ff1);
				avgA(ind, bCounter)=avgA(ind, bCounter)+newCell.nAP(ff1);
				avgN(ind, bCounter)=avgN(ind, bCounter)+1;
			end
		end
	end
end

avgI=avgI./avgN;
avgV=avgV./avgN;
avgA=avgA./avgN;

figure;
for c=1:4
	vv=avgV(c,:);
	ii=avgI(c,:);
	vv(isnan(vv))=[];
	ii(isnan(ii))=[];
	
	plot(ii, vv)
	hold on
end

figure;
for c=1:4
	aa=avgA(c,:);
	ii=avgI(c,:);
	aa(isnan(aa))=[];
	ii(isnan(ii))=[];
	
	plot(ii, aa)
	hold on
end

figure;
for c=1:4
	aa=avgA(c,:);
	vv=avgV(c,:);
	aa(isnan(aa))=[];
	vv(isnan(vv))=[];
	
	plot(vv, aa)
	hold on
end