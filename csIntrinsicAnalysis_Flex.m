function [ output_args ] = csIntrinsicAnalysis_Flex( cellList, varargin )
%csRunAnalysis_Flex Runs the EPSC analysis 
%	with flexible way to override variables

%% basic initialization and gui stuff
	% where is the data stored
	[~,prepath] = uiputfile('*.*','Select the path with the data folders', 'datapath.mat');
	if isnumeric(prepath)
		disp('Use must select an input path');
		return
	end
	
	% where do we write the output
	[~, savePath]=uiputfile('output.mat', 'Select output path');

	% some globals for posthoc analysis, if you want
	evalin('base', 'global newCell csAllCells csTableOut');
	global csTableRaw csTableSize newCell csTableOut csAllCells

	nameList={};
	nameListCounter=[];
	
    csTableOut={};
    csAllCells=[];
	newCell=[];
    	
	if nargin<1 || isempty(cellList)
		cellList=1:csTableSize(2)-1;
	end

	% save the input data
	if ~isnumeric(savePath) && ~isempty(savePath)
		save(fullfile(savePath, 'rawData.mat'), 'csTableRaw', 'csTableSize');
	else
		savePath=[];
	end
	
	% find the column names in the input and set up fields to hold them
	ind=[];
	for counter=1:length(csTableRaw(1,:))
		if ~isempty(csTableRaw{1,counter}) & ~isnan(csTableRaw{1,counter})
			ind.(matlab.lang.makeValidName(csTableRaw{1,counter}))=counter;
		end
	end
    
	% the fields that will go to the csv and xlsx tables at the end
    newCellFieldsToKeep={'acqNum', 'restMean', 'pulseI', 'pulseV', 'pulseRm', ...
		'nAP', 'sagV', 'reboundV', 'reboundAP'};
    
	
%% set up the variables to process data
	pulseList.p1=[-100 -75 -50 -25 10 20 30 40 50 60 70 80 90 100 150 200 250];
	pulseList.p2=[-100 0 20 40 60 80 100 150 200 250 300];

	pulseStart=1500;
	pulseEnd=2500;

	% where the RC check occurs
	checkPulseSize=-50; % pA we are in current clamp
	checkPulseStart=200;
	checkPulseEnd=500;
	usePulseListInExcel=0;

	% values for QC inclusion of individual sweeps
	maxRestSD=5;
 	maxRest=-30; 
 	minRest=-100;
 	minRm=50;
 	maxRm=1000;
	
	firstOnly=1;
	
	medianFilterSize=1;
	
    goodTracesToKeep=10;
    blockLength=goodTracesToKeep+4;
    colOffset=30;
    rowCounter=1;
    outputRowCounter=1;

	for c=1:2:length(varargin)
		if isempty(varargin{c+1})
			vStr='[]';
		else
			vStr=num2str(varargin{c+1});
		end
		
		disp(['Override: ' varargin{c} '=' vStr]);
		eval([varargin{c} '=' vStr ';']);
	end
	
    nCol=min(size(csTableRaw,2), colOffset);
    csTableOut(outputRowCounter, 1:nCol)=csTableRaw(rowCounter, 1:nCol);
	
    insertToCell('cycleName', size(csTableRaw,2)+1, 1)
    insertToCell('cyclePosition', size(csTableRaw,2)+2, 1)
    insertToCell('good_count', size(csTableRaw,2)+3, 1)
    insertToCell('bad_count', size(csTableRaw,2)+4, 1)
    insertToCell('keep_count', size(csTableRaw,2)+5, 1)

    for blockCounter=1:length(newCellFieldsToKeep)
        offset=colOffset+blockLength*(blockCounter-1);
        insertToCell(newCellFieldsToKeep{blockCounter},...
            offset, 1, goodTracesToKeep);
        insertToCell('Avg', offset+goodTracesToKeep, 1);
        insertToCell('Std', offset+goodTracesToKeep+1, 1);
    end
   
%% nested function to insert to cell
    function insertToCell(val, col, row, repeats)
        if nargin<3
            row=rowCounter;
        end
        if nargin<4
            repeats=1;
        end
        
        if length(val)==1 || ischar(val)
            if ~iscell(val)
                val={val};
            end
            if repeats==1
                csTableOut(row, col)=val;
            else
                csTableOut(row, col:col+repeats-1)=repmat(val, 1, repeats);
            end
        else
            if isnumeric(val)
                val=num2cell(val);
            end
            csTableOut(row, col:col+length(val)-1)=val;
        end
        
    end

%% nested function to insert into cell with avg and std
    function insertToCellWithStats(val, block, row)
        if nargin<3
            row=rowCounter;
        end
        
        offset=colOffset+blockLength*(block-1);
        if length(val)>goodTracesToKeep
            val=val(1:goodTracesToKeep);
        end
        
        insertToCell(val, offset, row);
        insertToCell(mean(val), offset+goodTracesToKeep, row);
        insertToCell(std(val), offset+goodTracesToKeep+1, row);
    end

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
    fullpathshort=[prepath num2str(csTableRaw{rowCounter, ind.Date})];
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
	newCell.firstOnly=firstOnly;
	
	newCell.acq=cell(1,nAcq); % store the full object for that acq sweep
    
    newCell.acqNum=nan(1, nAcq);
    newCell.cycleName=cell(1, nAcq);
    newCell.cyclePosition=nan(1, nAcq);
    newCell.pulsePattern=nan(1, nAcq);
    newCell.extraGain=nan(1, nAcq);
	if usePulseListInExcel
		newCell.pulseList=str2num(newCell.CurrentPulse);		
	else
		newCell.pulseList=pulseList.(['p' num2str(newCell.CurrentPulseID)]);
	end
    newCell.pulseListFirst=nan(1, length(newCell.pulseList));
	
   	newCell.traceQC=ones(1, nAcq);

	newCell.restMode=nan(1, nAcq);
    newCell.restMean=nan(1, nAcq);
    newCell.restMax=nan(1, nAcq);
    newCell.restMin=nan(1, nAcq);
	newCell.restSD=nan(1, nAcq);

	newCell.pulseRm=nan(1, nAcq);
    newCell.pulseV=nan(1, nAcq);
    newCell.nAP=nan(1, nAcq);
    newCell.traceQC=ones(1, nAcq);
	newCell.sagV=nan(1, nAcq);
	newCell.reboundV=nan(1, nAcq);
	newCell.reboundAP=nan(1,nAcq);

		
%% run through the acquisitions and calculate passive parameters
% use to do a first pass QC 
% examine resting potential, RC
    		
 	for sCounter=1:nAcq
        acqNum=sCounter+sStart-1;
        sFile=fullfile(fullpath, ['AD0_' num2str(acqNum) '.mat']);
		try
			a=load(sFile);
		catch
			disp('Trying without ww_gil');
			sFile=fullfile(fullpathshort, ['AD0_' num2str(acqNum) '.mat']);
			a=load(sFile);
		end
        
		if medianFilterSize>1
			a.(['AD0_' num2str(acqNum)]).data=...
				medfilt1(a.(['AD0_' num2str(acqNum)]).data, medianFilterSize);
		else
			a.(['AD0_' num2str(acqNum)]).data=a.(['AD0_' num2str(acqNum)]).data;
		end
		
        newCell.acq{sCounter}=a.(['AD0_' num2str(acqNum)]);
        
        newCell.acqNum(sCounter)=acqNum;
        newCell.cycleName{sCounter}=headerValue('state.cycle.cycleName');
        
        newCell.cyclePosition(sCounter)=headerValue('state.cycle.currentCyclePosition', 1);
        newCell.pulsePattern(sCounter)=headerValue('state.cycle.pulseToUse0', 1);
        newCell.extraGain(sCounter)=headerValue('state.phys.settings.extraGain0', 1);

		if newCell.cyclePosition(sCounter)<=length(newCell.pulseList)
			deltaI=newCell.pulseList(newCell.cyclePosition(sCounter))...
				*newCell.extraGain(sCounter);
		else
			disp([sFile ' ' num2str(acqNum) ' error in cycle']);
			deltaI=0;
		end
		newCell.pulseI(sCounter)=deltaI;
		
		acqData=a.(['AD0_' num2str(acqNum)]).data;
		if sCounter==1 % assume that the DAC sample rate doesn't change.
			acqRate=headerValue('state.phys.settings.inputRate', 1)/1000; % points per ms
			newCell.acqRate=acqRate;
			acqEndPt=length(acqData)-1;
			acqLen=length(acqData)/acqRate;
		end

		% define periods that are "baseline" and anylyze them 
		if isempty(checkPulseStart)
			notPulse=[SR(1, pulseStart-10) SR(pulseEnd+150, acqLen-1)]; 
		else
			if checkPulseStart>pulseStart % the RC check comes late
				notPulse=[SR(1, pulseStart-10) SR(pulseEnd+150, checkPulseStart-10)]; 
			else
				notPulse=[SR(checkPulseEnd+50, pulseStart-10) SR(pulseEnd+150, acqLen-1)]; 
			end
		end
		
		newCell.restMode(sCounter)=mode(round(notPulse));
		newCell.restMean(sCounter)=mean(notPulse);
		newCell.restMedian(sCounter)=median(notPulse);
		newCell.restSD(sCounter)=std(notPulse);
		newCell.restMin(sCounter)=min(notPulse);
		newCell.restMax(sCounter)=max(notPulse);

		% decide here reasons to reject a trace
		if newCell.restMode(sCounter)>maxRest
			newCell.traceQC(sCounter)=0;
		else
			ff=find(newCell.pulseList==deltaI);
			if isnan(newCell.pulseListFirst(ff))
				newCell.pulseListFirst(ff)=sCounter;
			end
		end
			

	end
			
	% We'll set some QC based on holding current and how it changes

%	deltaR=0.25;
% 	loR=max(avgStepRm*(1-deltaR), 50);
%  	hiR=min(avgStepRm*(1+deltaR), 500);

 	
% 	newCell.traceQC=...
% 		within(newCell.VrestMean, loV, hiV) && ...
% 		within(newCell.stepRmE, loR, hiR)...
% 		;
	
%  	newCell.traceQC=...
%  		within(newCell.restMean, minRest, maxRest) & ...
%  		within(newCell.stepRs, 0, maxRs) & ...
%  		within(newCell.stepRm, minRm, maxRm) & ...
%  		within(newCell.restMax-newCell.restMin, 0, 150) & ...
%  		within(newCell.restSD, 0, maxRestSD)...
%  		;


%% Run through the traces and analyze the data	

	if firstOnly
		acqToAna=sort(newCell.pulseListFirst(~isnan(newCell.pulseListFirst)));
	else
		acqToAna=1:nAcq;
	end
	
    for sCounter=acqToAna
        acqData=newCell.acq{sCounter}.data;
        acqRate=headerValue('state.phys.settings.inputRate', 1)/1000; % points per ms
		
		newCell.pulseV(sCounter)=mode(round(SR(pulseStart, pulseEnd)));

		deltaI=newCell.pulseI(sCounter);
  		if deltaI~=0
	        newCell.pulseRm(sCounter)=1000*...       % Rm in MOhm
	            (newCell.pulseV(sCounter)-newCell.restMean(sCounter))/deltaI;
			if deltaI<0
				minHyp=min(SR(pulseStart, pulseEnd));
				endHyp=mean(SR(pulseEnd-20,pulseEnd-1));
				newCell.sagV(sCounter)=endHyp-minHyp;
				newCell.reboundV(sCounter)=mean(SR(pulseEnd+20,pulseEnd+70)) ...
					-newCell.restMean(sCounter);
			end
		end
				
        newCell.pulseAP{sCounter}=ipAnalyzeAP(SR(pulseStart, pulseEnd));
		if isempty(newCell.pulseAP{sCounter})
			newCell.nAP(sCounter)=0;
		else
			newCell.nAP(sCounter)=newCell.pulseAP{sCounter}.nAP;
		end
		
        newCell.postAP{sCounter}=ipAnalyzeAP(SR(pulseEnd+1, floor(length(acqData)/acqRate)));
 		if isempty(newCell.postAP{sCounter})
			newCell.reboundAP(sCounter)=0;
		else
			newCell.reboundAP(sCounter)=newCell.postAP{sCounter}.nAP;
		end
		newCell.pulseAHP(sCounter)=min(SR(pulseEnd+1, floor(length(acqData)/acqRate)))-...
			newCell.restMean(sCounter);
		
	end


	
%% do some summary analysis
	goodTraces=find(newCell.traceQC);	
	badTraces=find(~newCell.traceQC);	
	
	nGood=length(goodTraces);

	
%% plot data
	newName=[newCell.mouseID '_' num2str(newCell.cellID)];
	disp(['plotting ' newName]);
    fFig=figure('name', newName);
	hold on

%% plot the rejected traces
	a1n=subplot(5,   3, [1:6]);
	title(a1n, 'Bad acquisitions');
	xlabel('time (ms)') 
	ylabel('I (pA)')
	hold on	
	for sCounter=badTraces
		acqData=newCell.acq{sCounter}.data;
		acqEndPt=length(acqData)-1;
		plot([0:acqEndPt]/acqRate, acqData);
	end


%% plot the good traces	
	a1=subplot(5,   3, [7:12]);
	title(a1, ['Good acquisitions']);
	xlabel('time (ms)') 
	ylabel('I (pA)')
	hold on
	
	denom=nAcq;	
	nGood=length(goodTraces);

	if nGood>denom/2
		newCell.QC=1;
	else
		newCell.QC=0;
	end
	
	for sCounter=goodTraces
		acqData=newCell.acq{sCounter}.data;
		acqEndPt=length(acqData)-1;		
		plot([0:acqEndPt]/acqRate, acqData);
	end
	

%%		
	a2=subplot(5, 4, [17:18]);
	yyaxis left
	scatter(newCell.pulseI(goodTraces), newCell.pulseV(goodTraces))
	set(gca, 'Ylim', [-100 0])
	
	ylabel('V (mV)')
	yyaxis right
	scatter(newCell.pulseI(goodTraces), newCell.nAP(goodTraces))
	title(a2, 'vs CURRENT')
	ylabel('# AP')
	xlabel('step (pA)') 

	a3=subplot(5, 4, [19:20]);
	scatter(newCell.pulseV(goodTraces), newCell.nAP(goodTraces));
	set(gca, 'Xlim', [-100 0])
	title(a3, 'vs VOLTAGE')
	xlabel('V (mV)') 
	ylabel('# AP')
	
	if ~isempty(savePath)
		saveas(fFig, fullfile(savePath, [newName 'Fig.fig']));
		print(fullfile(savePath, [newName 'FigPDF']),'-dpdf','-fillpage')
		save(fullfile(savePath, [newName '.mat']), 'newCell');
	end	
	
	
%% Run through the good ones and extract data
%	evalin('base', [newName '_ana=newCell;'])

    uniquePositions=sort(unique(newCell.cyclePosition));
    for cPos=uniquePositions
        outputRowCounter=outputRowCounter+1;
        csTableOut(outputRowCounter, 1:nCol)=csTableRaw(rowCounter, 1:nCol);

        bCount=0;
        gCount=0;
        keepCount=0;
        if ~isempty(badTraces)
            bCount=length(find(newCell.cyclePosition(badTraces)==cPos));
        end

        fPos=find(newCell.cyclePosition==cPos);
        insertToCell(newCell.cycleName{fPos(1)}, size(csTableRaw,2)+1, outputRowCounter)
        insertToCell(cPos, size(csTableRaw,2)+2, outputRowCounter)
       
        if ~isempty(goodTraces)
            gElements=goodTraces(newCell.cyclePosition(goodTraces)==cPos);
            gCount=length(gElements);
            if ~isempty(gElements)
  
                for blockCounter=1:length(newCellFieldsToKeep)
                    field=newCellFieldsToKeep{blockCounter};
                    insertToCellWithStats(newCell.(field)(gElements), blockCounter, outputRowCounter);
                end
            end
        end
        insertToCell(gCount,...
                    size(csTableRaw,2)+3, outputRowCounter)
        insertToCell(bCount,...
                    size(csTableRaw,2)+4, outputRowCounter)
        insertToCell(min(gCount, goodTracesToKeep),...
                    size(csTableRaw,2)+5, outputRowCounter)
	end
    

	nameIndex=find(ismember(nameList, newName));
	if isempty(nameIndex)
		nameList{end+1}=newName;
		nameListCounter(end+1)=1;
	else
		newName=[newName '_' num2str(nameListCounter(nameIndex))];
		nameListCounter(nameIndex)=nameListCounter(nameIndex)+1;
	end
	
	if ~isempty(savePath)
		print(fullfile(savePath, [newName]),'-dpdf','-fillpage')
		save(fullfile(savePath, [newName '.mat']), 'newCell');	
	end
	
    if isempty(csAllCells)
        csAllCells=newCell;
    else
        csAllCells(end+1)=newCell;
    end
end

T=cell2table(csTableOut);
if ~isempty(savePath)
 	save(fullfile(savePath, 'csTableOut.mat'), 'csTableOut');	
  	save(fullfile(savePath, 'csAllCells.mat'), 'csAllCells');	
    writetable(T, fullfile(savePath, 'analysisTable.csv'));
    writetable(T, fullfile(savePath, 'analysisTable.xlsx'));
end


end




