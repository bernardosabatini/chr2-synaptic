% SPfiles=dir('SP*.mat');
% 
% allCells=cell(1,length(SPfiles));
% for counter=1:length(SPfiles)
% 	cc=SPfiles(counter);
% 	load(cc.name)
% 	allCells{counter}=newCell;
% end


pCount=17;

avgV=zeros(1,pCount);
avgI=zeros(1,pCount);
avgA=zeros(1,pCount);
avgN=zeros(1,pCount);

dvList=nan(1,pCount);
mlList=nan(1,pCount);

avgV0=zeros(1,pCount);
avgI0=zeros(1,pCount);
avgA0=zeros(1,pCount);
avgN0=zeros(1,pCount);

avgV1=zeros(1,pCount);
avgI1=zeros(1,pCount);
avgA1=zeros(1,pCount);
avgN1=zeros(1,pCount);

avgV2=zeros(1,pCount);
avgI2=zeros(1,pCount);
avgA2=zeros(1,pCount);
avgN2=zeros(1,pCount);

bLow=[-100 0 19 49 100 150 200 250 300];
bHigh=[-100 0 31 61 100 150 200 250 300];
nb=length(bLow);
bFound=bLow*0;

allCounter=0;
for counter=1:length(csAllCells)
	newCell=csAllCells(counter);
	if length(newCell.pulseI)>=pCount
		for bCounter=1:nb
			ff=find((newCell.pulseI>bLow(bCounter)) & (newCell.pulseI<bHigh(bCounter)));
			if ~isempty(ff)
				ff1=ff(1);
			end
		end
	end
end

		
% 		
% 		
% 		
% 		if length(find(isnan(newCell.pulseV(1:pCount))))<4
% 			if 1 %&& newCell.pulseI(1)==-100
% 				nn=double(~isnan(newCell.pulseV(1:pCount)));
% 				vv=newCell.pulseV(1:pCount);
% 				ii=newCell.pulseI(1:pCount);
% 				aa=newCell.nAP(1:pCount);
% 				vv(isnan(vv))=0;
% 				aa(isnan(aa))=0;
% 				avgV=avgV+vv;
% 				avgA=avgA+aa;
% 				avgI=avgI+ii;
% 				avgN=avgN+nn;
% 				newCell.Injection
% 				dvList(counter)=newCell.DV;
% 				mlList(counter)=newCell.ML;
% 				allCounter=allCounter+1;
% 				if newCell.Injection(1)=='L'
% 					avgV2=avgV1+vv;
% 					avgA2=avgA2+aa;
% 					avgI2=avgI1+ii;
% 					avgN2=avgN1+nn;
% 				elseif	newCell.Injection(1)=='C'		
% 					avgV1=avgV1+vv;
% 					avgA1=avgA1+aa;
% 					avgI1=avgI1+ii;
% 					avgN1=avgN1+nn;
% 				elseif newCell.Injection(1)=='M'
% 					avgV0=avgV0+vv;
% 					avgA0=avgA0+aa;
% 					avgI0=avgI0+ii;
% 					avgN0=avgN0+nn;
% 				end
% 				
% 			end
			%		figure; plot(newCell.pulseI(1:pCount), newCell.pulseV(1:pCount))
% 		end
% 	end
% end

avgI=avgI/allCounter;
avgV=avgV./avgN;
avgA=avgA./avgN;

avgI2=avgI2./avgN2;
avgV2=avgV2./avgN2;
avgA2=avgA2./avgN2;

avgI1=avgI1./avgN1;
avgV1=avgV1./avgN1;
avgA1=avgA1./avgN1;

avgI0=avgI0./avgN0;
avgV0=avgV0./avgN0;
avgA0=avgA0./avgN0;

figure;plot(avgI, avgV)
hold on
plot(avgI, avgV0)
plot(avgI, avgV1)
plot(avgI, avgV2)

figure;plot(avgI, avgA)
hold on
plot(avgI, avgA0)
plot(avgI, avgA1)
plot(avgI, avgA2)

figure;plot(avgV, avgA)
hold on
plot(avgV0, avgA0)
plot(avgV1, avgA1)
plot(avgV2, avgA2)