function [ output_args ] = csRunAnalysis( cellList )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	[~, savePath]=uiputfile('output.mat', 'Select output path');

	evalin('base', 'global newCell csAllCells csAllCellsLabels');
	global csTableRaw csTableSize newCell 

	newCell=[];

	if nargin<1 || isempty(cellList)
		cellList=1:csTableSize(2)-1;
	end
	
	if ~isnumeric(savePath) && ~isempty(savePath)
		save(fullfile(savePath, 'rawData.mat'), 'csTableRaw', 'csTableSize');
	else
		savePath=[];
	end
	
	keepOnlyFirst=1;

	ind=[];
	
	for counter=1:length(csTableRaw(1,:))
		if ~isempty(csTableRaw{1,counter}) & ~isnan(csTableRaw{1,counter})
			ind.(matlab.lang.makeValidName(csTableRaw{1,counter}))=counter;
		end
	end
	
%% set up the variables to process data

	pulseStart=500; % where is optogenetic pulse
	pulseEnd=pulseStart+3;
	anaStart=pulseStart; % where will we analyze post-synaptic currents
	anaEnd=pulseStart+10; 
	chargeAnaEnd=pulseStart+30;
	anaPeakPre=1;
	anaPeakPost=2;
	fakeShift=-50; % how will be shift the analysis window for fake data
	
	checkPulseSize=-5;
	checkPulseStart=2500;
	checkPulseEnd=2550;
%	prepath='/Users/Bernardo/Dropbox (HMS)/BernardoGilShare/(1)PFcircuitPaper/fig1analysisCellCharic/';
%	prepath='/Volumes/BS Office/Dropbox (HMS)/BernardoGilShare/(1)PFcircuitPaper/3.fig2analysisConnec.example/';
	prepath='/Volumes/BS Office/Dropbox (HMS)/BernardoGilShare/(1)PFcircuitPaper/5.fig6analysisConnec1.example/D1vsD2/';
%	prepath='/Volumes/BS Office/Dropbox (HMS)/BernardoGilShare/(1)PFcircuitPaper/fig2analysisCellCharic/';

	

%% nested function to return a subrange of the data 
    function rang=SR(startS, endS)
        rang=acqData(floor(startS*acqRate):floor(endS*acqRate));
	end

%% nested function to extract a number from string
    function ns=extractNum(s)
        if isnumeric(s) 
            ns=s;
        elseif ischar(s)
			si=strfind(s, '_');
            if isempty(si)
                ns=str2double(s);
            else
                ns=str2double(s(si+1:end));
            end
        end
	end

%% nested function to return a value from the headerstring
    function hv=headerValue(sString, conv)
        if nargin<2
            conv=0;
        end
        hv=csHeaderValue(a.(['AD0_' num2str(acqNum)]).UserData.headerString, ...
            sString, conv);
	end
	
%% nested function that returns if a value falls within limits
%
	function ww=within(x, lo, hi)
		ww=(x>=lo) & (x<=hi);
		return
	end

%% nested function to return only non-nan entries
	function ap=nonNan(a)
		ap=a(~isnan(a));
		return
	end
%%


%% run through the cells in the list

for cellCounter=cellList 
    rowCounter=cellCounter+1; % the first row has headers
    fullpath=[prepath num2str(csTableRaw{rowCounter, ind.Date}) '/WW_Gil'];
    sStart=extractNum(csTableRaw{rowCounter, ind.SweepStart});
    sEnd=extractNum(csTableRaw{rowCounter, ind.SweepEnd});
    
    nAcq=sEnd-sStart+1;
	if isnan(nAcq)
		nAcq=0;
	end
    
    %% initialize the data object
	for label=fieldnames(ind)'
		newCell.(label{1})=csTableRaw{rowCounter, ind.(label{1})};
	end
	
	newCell.QC=1; % assume passes QC
	newCell.acqRate=0;
	
	newCell.acq=cell(1,nAcq); % store the full object for that acq sweep
    
    newCell.acqNum=nan(1, nAcq);
    newCell.cycleName=cell(1, nAcq);
    newCell.cyclePosition=nan(1, nAcq);
    newCell.pulsePattern=nan(1, nAcq);
    newCell.extraGain=nan(1, nAcq);
    
   	newCell.traceQC=ones(1, nAcq);

	newCell.restMode=nan(1, nAcq);
    newCell.restMean=nan(1, nAcq);
    newCell.restMax=nan(1, nAcq);
    newCell.restMin=nan(1, nAcq);
	newCell.restSD=nan(1, nAcq);

	newCell.stepCm=nan(1, nAcq);
    newCell.stepRs=nan(1, nAcq);
    newCell.stepRm=nan(1, nAcq);
    newCell.stepRt=nan(1, nAcq);
		
		
%% run through the acquisitions and calculate passive parameters
% use to do a first pass QC 
% examine resting potential, RC
    		
 	for sCounter=1:nAcq
        acqNum=sCounter+sStart-1;
        sFile=fullfile(fullpath, ['AD0_' num2str(acqNum) '.mat']);
        a=load(sFile);
        
        newCell.acq{sCounter}=a.(['AD0_' num2str(acqNum)]);
        
        newCell.acqNum(sCounter)=acqNum;
        newCell.cycleName{sCounter}=headerValue('state.cycle.cycleName');
        
        newCell.cyclePosition(sCounter)=headerValue('state.cycle.currentCyclePosition', 1);
        newCell.pulsePattern(sCounter)=headerValue('state.cycle.pulseToUse0', 1);
        newCell.extraGain(sCounter)=headerValue('state.phys.settings.extraGain0', 1);
				
		acqData=a.(['AD0_' num2str(acqNum)]).data;
		if sCounter==1 % assume that the DAC sample rate doesn't change.
			acqRate=headerValue('state.phys.settings.inputRate', 1)/1000; % points per ms
			newCell.acqRate=acqRate;
			acqEndPt=length(acqData)-1;
			acqLen=length(acqData)/acqRate;
		end

		% define periods that are "baseline" and anylyze them 
		notPulse=[SR(1, pulseStart-10) SR(anaEnd+100, checkPulseStart-10)]; 
		newCell.restMode(sCounter)=mode(round(notPulse));
		newCell.restMean(sCounter)=mean(notPulse);
		newCell.restSD(sCounter)=std(notPulse);
		newCell.restMin(sCounter)=min(notPulse);
		newCell.restMax(sCounter)=max(notPulse);

		[newCell.stepRm(sCounter), ...
			newCell.stepRs(sCounter), ...
			newCell.stepCm(sCounter)] = ...
			csAnalyzeRC(SR(checkPulseStart,checkPulseEnd)-newCell.restMean(sCounter), ...
			checkPulseSize, acqRate);
	    newCell.stepRt(sCounter)=newCell.stepRm(sCounter)+newCell.stepRs(sCounter);
		
	end
			
	avgRest=median(nonNan(newCell.restMean));
	avgStepRm=median(nonNan(newCell.stepRm));
	avgStepRs=median(nonNan(newCell.stepRs));
	avgStepCm=median(nonNan(newCell.stepCm));


		
	% We'll set some QC based on holding current and how it changes
 	hiRest=300; 
 	loRest=-300;

	deltaR=0.15;
%  	loR=avgStepRm*(1-deltaR);
%  	hiR=avgStepRm*(1+deltaR);
 	loR=100;
 	hiR=500;

 
 	newCell.traceQC=...
 		within(newCell.restMean, loRest, hiRest) & ...
 		within(newCell.stepRs, 0, 30) & ...
 		within(newCell.stepRm, loR, hiR) & ...
 		within(newCell.restMax-newCell.restMin, 0, 100) & ...
 		within(newCell.restSD, 0, 10)...
 		;


 	for sCounter=1:nAcq
		acqData=newCell.acq{sCounter}.data;	
	
		newCell.pscPeakBL(sCounter)=mean(SR(anaStart+fakeShift,...
			anaEnd+fakeShift));	
		[pk, indp]=min(SR(anaStart, anaEnd));
		newCell.pscPeak(sCounter)=pk-newCell.pscPeakBL(sCounter);
		xStart=(indp-1)/acqRate+anaStart;
		newCell.pscPeriPeak(sCounter)=mean(SR(xStart+anaPeakPre, ...
			xStart+anaPeakPost)) - newCell.pscPeakBL(sCounter);

		newCell.pscFakePeakBL(sCounter)=mean(SR(anaStart+2*fakeShift,...
			anaEnd+2*fakeShift));
		[pk, indp]=min(...
			SR(anaStart+fakeShift, anaEnd+fakeShift));
		newCell.pscFakePeak(sCounter)=pk-newCell.pscFakePeakBL(sCounter);
		xStart=(indp-1)/acqRate+anaStart+fakeShift;
		newCell.pscPeriFakePeak(sCounter)=mean(SR(xStart+anaPeakPre, ...
			xStart+anaPeakPost)) - newCell.pscFakePeakBL(sCounter);

		newCell.pscCharge(sCounter)=(chargeAnaEnd-anaStart)*...
			(mean(SR(anaStart, chargeAnaEnd))-newCell.pscPeakBL(sCounter));
		newCell.pscFakeCharge(sCounter)=(chargeAnaEnd-anaStart)*...
			(mean(SR(anaStart+fakeShift, chargeAnaEnd+fakeShift))...
			-newCell.pscFakePeakBL(sCounter));	
	
	end


	
%% do some summary analysis
	goodTraces=find(newCell.traceQC);	
	badTraces=find(~newCell.traceQC);	
	
	nGood=length(goodTraces);
	if nGood>nAcq/2
		newCell.QC=1;
	else
		newCell.QC=0;
	end
	
	% calculate the average
	avgData=[];
	for sCounter=goodTraces
		acqData=newCell.acq{sCounter}.data;
		if isempty(avgData)
			avgData=acqData/nGood;
		else
			avgData=avgData+acqData/nGood;
		end
	end
	



%% extract PSC peaks based on where the peak of the average trace is
	acqData=avgData; % load up the data
	
	if ~isempty(avgData)
		newCell.avgData=avgData;
		newCell.avgRestMean=mean(newCell.restMean);

		[newCell.avgStepRs, ...
			newCell.avgStepRm, ...
			newCell.avgStepCm] = ...
			csAnalyzeRC(SR(checkPulseStart,checkPulseEnd)-newCell.avgRestMean, ...
			checkPulseSize, acqRate);

		[~, indp]=min(SR(anaStart, anaEnd));
		xStart=(indp-1)/acqRate+anaStart;

		[~, indp]=min(...
				SR(anaStart+fakeShift, anaEnd+fakeShift));
		xStartFake=(indp-1)/acqRate+anaStart+fakeShift;

		for sCounter=1:nAcq
			acqData=newCell.acq{sCounter}.data;
			newCell.pscPeriAvgPeak(sCounter)=mean(SR(xStart+anaPeakPre, ...
				xStart+anaPeakPost)) - newCell.pscPeakBL(sCounter);
			newCell.pscPeriAvgFakePeak(sCounter)=mean(SR(xStartFake+anaPeakPre, ...
				xStartFake+anaPeakPost)) - newCell.pscFakePeakBL(sCounter);		
		end

	else
		newCell.avgRestMean=nan;
		newCell.avgStepRs=nan;
		newCell.avgStepRm=nan;
		newCell.avgStepCm=nan;
		newCell.QC=0;
		newCell.pscPeriAvgPeak(1:nAcq)=NaN;
		newCell.pscPeriAvgFakePeak(1:nAcq)=NaN;
	end
	

	
%% plot data
	newName=[newCell.mouseID '_' num2str(newCell.cellID)];
	disp(['plotting ' newName]);
    fFig=figure('name', newName);
	hold on

%% plot the rejected traces
	a1n=subplot(5,   3, [1:3]);
	title(a1n, 'Bad acquisitions');
	xlabel('time (ms)') 
	ylabel('I (pA)')
	hold on	
	for sCounter=badTraces
		acqData=newCell.acq{sCounter}.data;
		plot([0:acqEndPt]/acqRate, acqData);
	end


%% plot the good traces	
	a1=subplot(5,   3, [4:6]);
	title(a1, ['Good acquisitions']);
	xlabel('time (ms)') 
	ylabel('I (pA)')
	hold on
	
	denom=nAcq;	
	nGood=length(goodTraces);

	avgData=[];
	if nGood>denom/2
		newCell.QC=1;
	else
		newCell.QC=0;
	end
	
	for sCounter=goodTraces
		acqData=newCell.acq{sCounter}.data;
		plot([0:acqEndPt]/acqRate, acqData);
	end
	
	% plot again with baselining
	a2=subplot(5,   3, [7:9]);
	title(a2, ['Good acquisitions BL']);
	xlabel('time (ms)') 
	ylabel('I (pA)')
	hold on

	for sCounter=goodTraces
		acqData=newCell.acq{sCounter}.data-newCell.restMean(sCounter);
		plot([0:acqEndPt]/acqRate, acqData);
	end
	set(a2, 'Ylim', [1.2*min(newCell.pscPeak) -1.2*min(newCell.pscPeak)]);
	set(a2, 'Xlim', [min(anaStart+fakeShift, anaEnd-fakeShift) max(anaStart+fakeShift, anaEnd-fakeShift)]);
	
	% plot cell params
	a4rs=subplot(5,   3, 10);
	title(a4rs, ['Rs']);
	xlabel(a4rs,'acq') 
	ylabel(a4rs,'MO')
	hold on
	a4rm=subplot(5,   3, 11);
	title(a4rm, ['Rm']);
	xlabel(a4rm, 'acq') 
	hold on
	a4im=subplot(5,   3, 12);
	title(a4im, ['Im']);
	xlabel(a4im,'acq') 
	ylabel(a4im,'pA')
	hold on
	if ~isempty(goodTraces)
		plot(a4rs, goodTraces,newCell.stepRs(goodTraces), 'go')
		plot(a4rm, goodTraces,newCell.stepRm(goodTraces), 'go')
		plot(a4im, goodTraces,newCell.restMean(goodTraces), 'go')
	end
	if ~isempty(badTraces)
		plot(a4rs, badTraces,newCell.stepRs(badTraces), 'ro')
		plot(a4rm, badTraces,newCell.stepRm(badTraces), 'ro')
		plot(a4im, badTraces,newCell.restMean(badTraces), 'ro')
	end
	
	
	% plot psc params
	aPscPeak=subplot(5,   3, 13);
	title(aPscPeak, ['peak']);
	xlabel(aPscPeak,'acq') 
	ylabel(aPscPeak,'pA')
	hold on
	aPeriPeak=subplot(5,  3, 14);
	title(aPeriPeak, ['peri peak']);
	xlabel(aPeriPeak, 'acq') 
	hold on
	aCharge=subplot(5,  3, 15);
	title(aCharge, ['charge']);
	xlabel(aCharge,'acq') 
	hold on
	if ~isempty(goodTraces)
		plot(aPscPeak, goodTraces, newCell.pscPeak(goodTraces), 'go')
		plot(aPeriPeak, goodTraces,newCell.pscPeriAvgPeak(goodTraces), 'go')
		plot(aCharge, goodTraces, newCell.pscCharge(goodTraces), 'go')
		plot(aPscPeak, goodTraces, newCell.pscFakePeak(goodTraces), 'bo')
		plot(aPeriPeak, goodTraces, newCell.pscPeriAvgFakePeak(goodTraces), 'bo')
		plot(aCharge, goodTraces, newCell.pscFakeCharge(goodTraces), 'bo')
	end		

%% Run through the good ones and extract data
	evalin('base', [newName '_ana=newCell;'])
% 

	fieldsToKeep={'newName', ...
		'newCell.ML', ...
		'newCell.DV', ...
		'newCell.avgRestMean', ...
		'newCell.avgStepRm', ...
		'newCell.avgStepRs', ...
		'newCell.avgStepCm', ...
		
 	csAllCells{cellCounter,1}=newName;
 	csAllCells{cellCounter,2}=newCell.ML;
 	csAllCells{cellCounter,3}=newCell.DV;
 	csAllCells{cellCounter,5}=newCell.avgRestMean;
 	csAllCells{cellCounter,6}=newCell.avgStepRm;
 	csAllCells{cellCounter,7}=newCell.avgStepRs;
 	csAllCells{cellCounter,8}=newCell.avgStepCm;
	
	allCp=sort(unique(newCell.cyclePosition));
	for cpCounter=1:length(allCp)
		cp=
		
	end
	
	if ~isempty(savePath)
		print(fullfile(savePath, [newName]),'-dpdf','-fillpage')
		save(fullfile(savePath, [newName '.mat']), 'newCell');	
	end
end

% if ~isempty(savePath)
% 	save(fullfile(savePath, 'csAllCells.mat'), 'csAllCells', 'csAllCellsLabels');	
% end


end




