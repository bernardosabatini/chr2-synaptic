function [ output_args ] = phAnalyzeIntrinsicProperties( cellList, varargin )
%phAnalyzeIntrinsicProperties 

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
    newCellFieldsToKeep={'acqNum', 'cyclePosition', 'extraGain', 'restMean', 'pulseI', 'pulseV', 'pulseRm', ...
		'nAP', 'sagV', 'reboundV', 'reboundAP', 'checkPulseRend', 'checkPulseRpeak', 'checkPulseTau'};
    
	
%% set up the variables to process data

	% this is a very important selection that determines how the current
	% pulse amplitude will be determined for each trace
	pulseAmplitudeMode=3; 
	% 1 - use the values stored in pulseListbyCyclePosition and index them by the cycle
	% position. e.g. cyclePostition=2 means take the second value in the
	% list
	% 2 - use the values in the excel sheet. THis option will take the list 
	% from the "CurrentPulse" column of the excel sheet and index by the
	% cycle position.  i.e. the column can have "-100, 0, 50, 100" in it 
	% cycle posiiton=3 will retreive 50
	% 3 - use the values stored in pulseListByPatternNumber using the pulse
	% pattern number instead of the cycle position.  The pulsePatternNumber
	% will be used to look up an index in the first column of pulseListByPatternNumber
	% the value in the second column is the current pulse amplitude;
	% 4 - extract the value from the headerString.  This is obviously the
	% preferred method but it was not implemented correctly in some old
	% data
	
	pulseListByCyclePosition=[-100 -75 -50 -25 10 20 30 40 50 60 70 80 90 100 150 200 250];
	pulseListByPatternNumber=[...
		[50 -100]; ...
		[51 -75];
		[52 -50];
		[53 -25]; 
		[54 -10];
		[55 10];
		[56 -20];
		[57 30];
		[58 40];
		];
	% pulseListbyPosition=[-100 0 20 40 60 80 100 150 200 250 300];

	% Flag to determine if only the first time occurence of a pulse amplitude
	% should be used or if all repeated occurences should be processed 
	firstOnly=0; 

	pulseStart=1500;
	pulseEnd=2500;

	% where the RC check occurs
	% set checkPulseStart=[] if there is no standard RC check pulse in
	% every trace
	checkPulseSize=-50; % pA we are in current clamp
	checkPulseStart=200;
	checkPulseEnd=500;

	% values for QC inclusion of individual sweeps
	maxRestSD=5; % max SD of the resting voltage to pass QC
 	minRm=50; % min Rm to pass QC
 	maxRm=1000; % max Rm to pass QC
	
	maxVm=-50; % max Vm for inclusion and to pass QC
	minVm=-90; % min Vm for inclusion and to pass QC
	
	initialDeltaVmSweeps=5; % How many of the first sweeps should we use to 
							% determine the starting Vm
	maxDeltaVm=5; % How far can the Vm move from the above value and still be included
	maxFractionalDeltaRm=0.2; % What fractional change in the Rm will be accept
	
	minFractionGoodToKeepCell=0.5;	% What fraction of the traces need to pass QC for
									% the cell to pass QC
		
	medianFilterSize=1; % filter the data?
	
%% Set up some parameters to format the output .CSV file and Overwrite the 
% default parameters set above using any parameters values specified 
% by the user in the varargin

	% 20 is the assumed maximum number of traces - used to format the
	% output spreadsheet only
    goodTracesToKeep=20;
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
	if ~firstOnly
	    insertToCell('cyclePosition', size(csTableRaw,2)+2, 1)
	end
    insertToCell('good_count', size(csTableRaw,2)+3, 1)
    insertToCell('bad_count', size(csTableRaw,2)+4, 1)
    insertToCell('keep_count', size(csTableRaw,2)+5, 1)

    for blockCounter=1:length(newCellFieldsToKeep)
        offset=colOffset+blockLength*(blockCounter-1);
        insertToCell(newCellFieldsToKeep{blockCounter},...
            offset, 1, goodTracesToKeep);
		if ~firstOnly
			insertToCell('Avg', offset+goodTracesToKeep, 1);
			insertToCell('Std', offset+goodTracesToKeep+1, 1);
		end
    end
   
%% nested function to insert data to a cell in the output table
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

%% nested function to insert data into a cell with its avg and std too
    function insertToCellWithStats(val, block, row, stopStats)
        if nargin<3
            row=rowCounter;
		end
		if nargin<4
			stopStats=0;
		end
        
        offset=colOffset+blockLength*(block-1);
        if length(val)>goodTracesToKeep
            val=val(1:goodTracesToKeep);
        end
        
        insertToCell(val, offset, row);
		if ~stopStats
			insertToCell(mean(val), offset+goodTracesToKeep, row);
			insertToCell(std(val), offset+goodTracesToKeep+1, row);
		end
    end

%% nested function to return a subrange of the data 
    function rang=SR(startS, endS)
        rang=acqData(floor(startS*acqRate):floor(endS*acqRate));
	end

%% nested function to extract a number from string dealing with possible multiple formats
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
		newCell=[];
		for label=fieldnames(ind)'
			newCell.(label{1})=csTableRaw{rowCounter, ind.(label{1})};
		end

		newCell.QC=1; % assume cell passes QC
		newCell.acqRate=0;
		newCell.firstOnly=firstOnly;
		newName=[newCell.mouseID '_' num2str(newCell.cellID)];

		newCell.acq=cell(1,nAcq); % store the full object for that acq sweep
		newCell.acqNum=nan(1, nAcq);
		newCell.cycleName=cell(1, nAcq);
		newCell.cyclePosition=nan(1, nAcq);
		newCell.pulsePattern=nan(1, nAcq);
		newCell.extraGain=nan(1, nAcq);
		switch pulseAmplitudeMode
			case 1
				newCell.pulseList=pulseListByCyclePosition;
			case 2
				if isempty(newCell.CurrentPulse)
					disp('**** No pulse amplitude list found in spreadsheet')
					disp('   make sure there is a column titled CurrentPulse')
					error('can not process amplitude current pulses');
				else
					newCell.pulseList=str2num(newCell.CurrentPulse);
					if isempty(newCell.pulseList)
						newCell.pulseList=str2num(newCell.CurrentPulse(1:end-1));
					end
					if isempty(newCell.pulseList)
						disp(newCell.CurrentPulse);
						error('can not process amplitude current pulses');
					end

				end
			case 3
				newCell.pulseList=[]; % look it up dynamically below		
			case 4
				newCell.pulseList=[]; % look it up dynamically below		
			otherwise
				error(['pulseAmplitudeMode = ' num2str(pulseAmplitudeMode) ...
					' which is an unknown mode']);
		end
			
		newCell.pulseListFirst=[];

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
		newCell.pulseAHP=nan(1,nAcq);
		newCell.postAP=cell(1,nAcq);
		newCell.pulseAP=cell(1,nAcq);	

		newCell.checkPulseRend=nan(1,nAcq);
		newCell.checkPulseRpeak=nan(1,nAcq);
		newCell.checkPulseTau=nan(1,nAcq);


	%% run through the acquisitions and calculate passive parameters
	% use to do a first pass QC 
	% examine resting potential, RC

		foundSomething=0;
		for sCounter=1:nAcq
			acqNum=sCounter+sStart-1;
			sFile=fullfile(fullpath, ['AD0_' num2str(acqNum) '.mat']);
			failedLoad=1;
			try
				a=load(sFile);
				failedLoad=0;
			catch
				disp('First load failed. Trying without ww_gil');
			end
			if failedLoad
				try
					sFile=fullfile(fullpathshort, ['AD0_' num2str(acqNum) '.mat']);
					a=load(sFile);
					failedLoad=0;
				catch
					disp(['Second load failed.  Will skip ' sFile]);
				end
			end

			if ~failedLoad
				foundSomething=1;

				if medianFilterSize>1
					a.(['AD0_' num2str(acqNum)]).data=...
						medfilt1(a.(['AD0_' num2str(acqNum)]).data, medianFilterSize);
				else
					a.(['AD0_' num2str(acqNum)]).data=a.(['AD0_' num2str(acqNum)]).data;
				end
				
				a.(['AD0_' num2str(acqNum)]).headerString= ...
					a.(['AD0_' num2str(acqNum)]).UserData.headerString;
				
				newCell.acq{sCounter}=a.(['AD0_' num2str(acqNum)]);

				newCell.acqNum(sCounter)=acqNum;
				newCell.cycleName{sCounter}=headerValue('state.cycle.cycleName');

				newCell.cyclePosition(sCounter)=headerValue('state.cycle.currentCyclePosition', 1);
				newCell.pulsePattern(sCounter)=headerValue('state.cycle.pulseToUse0', 1);
				newCell.extraGain(sCounter)=headerValue('state.phys.settings.extraGain0', 1);
				
				if (pulseAmplitudeMode==1) || (pulseAmplitudeMode==2) % value is stored in pulse list by now
					if newCell.cyclePosition(sCounter)<=length(newCell.pulseList)
						deltaI=newCell.pulseList(newCell.cyclePosition(sCounter))...
							*newCell.extraGain(sCounter);
					else
						disp(['****' sFile ' ' num2str(acqNum) ' error in cycle']);
						disp(['   cycle position is beyond the entried in the current pulse list']);
						disp(['   Setting pulse amplitude to 0']);
						deltaI=0;
					end
				elseif (pulseAmplitudeMode==3) % look it up by pulse pattern in the list set by user above
					ff=find(pulseListByPatternNumber(:,1)==newCell.pulsePattern(sCounter));
					if isempty(ff)
						disp(['****' sFile ' ' num2str(acqNum) ' error in cycle']);
						disp(['   pulse pattern #' num2str(newCell.pulsePattern(sCounter)) ...
							' is not in pulseListByPatternNumber'])
						disp(['   Setting pulse amplitude to 0']);
						deltaI=0;						
					else
						if length(ff)>1
							disp(['****' sFile ' ' num2str(acqNum) ' error in cycle']);
							disp(['   pulse pattern #' num2str(newCell.pulsePattern(sCounter)) ...
								' has multiple entreis in pulseListByPatternNumber'])
							disp(['   Setting pulse amplitude to the first entry']);
						end
						deltaI=pulseListByPatternNumber(ff(1), 2)...
							*newCell.extraGain(sCounter);
					end
				elseif (pulseAmplitudeMode==4)
					disp(['**** pulseAmplitudeMode 4 is not implemented']);
					disp(['   Setting pulse amplitude to 0']);
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
					rPeak=NaN;
					rEnd=NaN;
					tau=NaN;
				else
					if checkPulseStart>pulseStart % the RC check comes late
						notPulse=[SR(1, pulseStart-10) SR(pulseEnd+150, checkPulseStart-10)]; 
					else
						notPulse=[SR(1, checkPulseStart-1) SR(checkPulseEnd+50, pulseStart-10) SR(pulseEnd+150, acqLen-1)]; 
					end
					[rPeak, rEnd, tau]=csCurrentClampPulseAnalysis(SR(checkPulseStart, checkPulseEnd)- mean(notPulse), acqRate, checkPulseSize);
				end

				newCell.restMode(sCounter)=mode(round(notPulse));
				newCell.restMean(sCounter)=mean(notPulse);
				newCell.restMedian(sCounter)=median(notPulse);
				newCell.restSD(sCounter)=std(notPulse);
				newCell.restMin(sCounter)=min(notPulse);
				newCell.restMax(sCounter)=max(notPulse);
				newCell.checkPulseRpeak(sCounter)=rPeak;
				newCell.checkPulseRend(sCounter)=rEnd;
				newCell.checkPulseTau(sCounter)=tau;

				% We'll set some QC 
				newCell.traceQC(sCounter)=...
					within(newCell.restMode(sCounter), minVm, maxVm) & ...
					within(newCell.restSD(sCounter), 0, maxRestSD);
				
				if ~isempty(checkPulseStart) 
					newCell.traceQC(sCounter)=...
						newCell.traceQC(sCounter) & ...
						within(newCell.checkPulseRend(sCounter), minRm, maxRm);
				end	

				% ths trace passes QC.  Is it the first one of that pulse amplitude?
				if newCell.traceQC(sCounter)
					if isempty(newCell.pulseListFirst)
						newCell.pulseListFirst=[newCell.pulseI(sCounter) sCounter];
					else						
						if isempty(find(newCell.pulseListFirst(:,1)==newCell.pulseI(sCounter),1))
							newCell.pulseListFirst(end+1, :)=[newCell.pulseI(sCounter) sCounter];
						end
					end
				end
			end
		end

		if ~foundSomething
			disp(['  NO DATA FOR ' newName]);
		else

			% We'll do some more QC 
			% if volage more than maxDeltaVm mV from that of the first sweep, reject
			ppp=find(newCell.traceQC);
			ppp=ppp(1:min(initialDeltaVmSweeps, length(ppp)));
			v0=median(newCell.restMode(ppp));
			
			newCell.traceQC=...
				newCell.traceQC & ...
				within(newCell.restMode, v0-maxDeltaVm, v0+maxDeltaVm);
			if ~isempty(checkPulseStart) 
				r0=median(newCell.checkPulseRend(ppp));
				newCell.traceQC=...
					newCell.traceQC & ...
					within(newCell.checkPulseRend, ...
						(1-maxFractionalDeltaRm)*r0, ...
						(1+maxFractionalDeltaRm)*r0);		
			end

		%% do some summary analysis
			goodTraces=find(newCell.traceQC);	
			badTraces=find(~newCell.traceQC);	

			nGood=length(goodTraces);

		%% Run through the traces and analyze the data	

			allAcq=find(~isnan(newCell.restMedian));

			if firstOnly
				if isempty(newCell.pulseListFirst)
					acqToAna=[];
				else
					acqToAna=sort(newCell.pulseListFirst(:,2))';
				end
				goodTraces=intersect(acqToAna, goodTraces);
				badTraces=setdiff(allAcq, goodTraces);
			else
				acqToAna=allAcq;
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

			newCell.checkPulseRpeakMean=mean(newCell.checkPulseRpeak(goodTraces));
			newCell.checkPulseRendMean=mean(newCell.checkPulseRend(goodTraces));
			newCell.checkPulseTauMean=mean(newCell.checkPulseTau(goodTraces));



		%% plot data
			disp(['plotting ' newName]);
			fFig=figure('name', newName);
			hold on

		%% plot the rejected traces
			a1n=subplot(2,   2, 1);
			title(a1n, 'Not used');
			xlabel('time (ms)') 
			ylabel('I (pA)')
			hold on	
			for sCounter=badTraces
				acqData=newCell.acq{sCounter}.data;
				acqEndPt=length(acqData)-1;
				plot([0:acqEndPt]/acqRate, acqData);
			end


		%% plot the good traces	
			a1=subplot(2,   2, 2);
			title(a1, ['Used']);
			xlabel('time (ms)') 
			ylabel('I (pA)')
			hold on

			denom=length(acqToAna);	
			nGood=length(goodTraces);

			if nGood>denom*minFractionGoodToKeepCell
				newCell.QC=1;
			else
				newCell.QC=0;
			end

			if ~isempty(goodTraces)
				for sCounter=goodTraces
					acqData=newCell.acq{sCounter}.data;
					acqEndPt=length(acqData)-1;		
					plot([0:acqEndPt]/acqRate, acqData);
				end
			end

		%%	resistances
			% plot cell params
			a4rPeak=subplot(4,  3, 7);
			title(a4rPeak, ['RC R']);
			xlabel(a4rPeak,'acq') 
			ylabel(a4rPeak,'MO')
			hold on
			a4rPulse=subplot(4, 3, 8);
			title(a4rPulse, ['Pulse R']);
			xlabel(a4rPulse, 'acq') 
			hold on
			a4tau=subplot(4, 3, 9);
			title(a4tau, ['RC Tau']);
			xlabel(a4tau,'acq') 
			ylabel(a4tau,'ms')
			hold on

			if ~isempty(goodTraces)
				plot(a4rPeak, goodTraces, newCell.checkPulseRpeak(goodTraces), 'go')
				plot(a4rPeak, goodTraces, newCell.checkPulseRend(goodTraces), 'gx')
				plot(a4rPulse, goodTraces, newCell.pulseRm(goodTraces), 'go')		
				plot(a4tau, goodTraces, newCell.checkPulseTau(goodTraces), 'go')
			end
			if ~isempty(badTraces)
				plot(a4rPeak, badTraces,newCell.checkPulseRpeak(badTraces), 'ro')
				plot(a4rPeak, badTraces,newCell.checkPulseRend(badTraces), 'rx')
				plot(a4rPulse, badTraces, newCell.pulseRm(badTraces), 'go')		
				plot(a4tau, badTraces,newCell.checkPulseTau(badTraces), 'ro')
			end


			a3=subplot(5, 3, [14]);
			title(a3, 'vs CURRENT')
			scatter(newCell.pulseI(goodTraces), newCell.nAP(goodTraces))
			hold on
			scatter(newCell.pulseI(goodTraces), newCell.reboundAP(goodTraces), 'rx');	
			ylabel('# AP')
			xlabel('step (pA)') 

			a4=subplot(5, 3, [15]);
			title(a4, 'vs VOLTAGE')
			scatter(newCell.pulseV(goodTraces), newCell.nAP(goodTraces));
			hold on
			scatter(newCell.pulseV(goodTraces), newCell.reboundAP(goodTraces), 'rx');
			set(gca, 'Xlim', [-100 0])
			xlabel('V (mV)') 
			ylabel('# AP')

			a2=subplot(5, 3, [13]);
			title(a2, 'vs CURRENT')
			yyaxis left
			scatter(newCell.pulseI(goodTraces), newCell.pulseV(goodTraces))
			set(gca, 'Ylim', [-100 0])
			ylabel('V (mV)')
			yyaxis right
			scatter(newCell.pulseI(goodTraces), newCell.pulseRm(goodTraces))
			ylabel('R (MO)')

			%% Run through the good ones and extract data

			evalin('base', [newName '_ana=newCell;'])
			drawnow
			
			uniquePositions=sort(unique(newCell.cyclePosition));

			if firstOnly % if we only keep the first time through, then put it all on a line
				outputRowCounter=outputRowCounter+1; % one row per cell
				csTableOut(outputRowCounter, 1:nCol)=csTableRaw(rowCounter, 1:nCol);
				insertToCell(newCell.cycleName{1}, nCol+1, outputRowCounter)
				bCount=length(find(newCell.traceQC==0));
				gCount=length(find(newCell.traceQC==1));
				keepCount=length(goodTraces);
				insertToCell(gCount, nCol+3, outputRowCounter)
				insertToCell(bCount, nCol+4, outputRowCounter)
				insertToCell(keepCount, nCol+5, outputRowCounter)
				
				if ~isempty(goodTraces)
					for blockCounter=1:length(newCellFieldsToKeep)
						field=newCellFieldsToKeep{blockCounter};
						insertToCellWithStats(newCell.(field)(goodTraces), blockCounter, outputRowCounter);
					end
				end					
			else
				for cPos=uniquePositions
					outputRowCounter=outputRowCounter+1; % one row person cycle position
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
			end

			%% Do some saving
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

% stupid that they make me add this follwoing end
end




