function [ results ] = ipAnalyzeAP(dData, acqRate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    [gUp, gDown]=ipFindXings(dData, 0, 1); % find the zero crossings
    gUp=floor(gUp);
    gDown=ceil(gDown);
    
	if nargin<2
		acqRate=10;
	end
	
    if isempty(gUp)
        results=[];
%        disp([   'No action potentials found']);
        return
    end
    
    xNum=min(length(gUp), length(gDown));
%    disp([   num2str(xNum) ' action potentials found']);

    results.nAP=xNum;
    
    results.AP_peak_V=zeros(1, xNum);
    results.AP_peak_time=zeros(1, xNum);
    results.AP_AHP_V=zeros(1,xNum);
    results.AP_thresh_V=zeros(1, xNum);
    results.AP_thresh_time=zeros(1, xNum);   
    results.AP_HW_V=zeros(1, xNum);
    results.AP_HW=zeros(1, xNum);
    results.AP_max_dVdT=zeros(1, xNum);
	
    g2=gradient(dData);
    g3=gradient(g2);
    
	if length(gDown)>length(gUp)
		if gDown(1)<gUp(1)
			gDown=gDown(2:(length(gUp)+1));
		else
			disp('problem with gDown');
		end
	end
    gDown(end+1)=length(dData);
    gUp(end+1)=length(dData);	
    

    lastMin=1;
    for counter=1:xNum
        [results.AP_peak_V(counter), Imax]=max(dData(gUp(counter):gDown(counter)));
        Imax=Imax+gUp(counter)-1;
        results.AP_peak_time(counter)=Imax;
        [results.AP_AHP_V(counter), Imin]=...
			min(...
			dData(gDown(counter):...
			min(gDown(counter)+10*acqRate,gUp(counter+1)))); % find the min between this AP and the next or 10 ms later, whichever comes first
        Imin=Imin+gDown(counter)-1;
        HW_V=(results.AP_peak_V(counter)-results.AP_AHP_V(counter))/2+results.AP_AHP_V(counter);
        [ggUp, ggDown]=ipFindXings(dData(lastMin:Imin), 0, 1);
		if length(ggDown)>length(ggUp)
			if ggDown(1)<ggUp(1)
				ggDown=ggDown(2:(length(ggUp)+1));
			else
				disp('problem with ggDown');
			end
		end
        results.AP_HW(counter)=ggDown-ggUp;
        results.AP_HW_V(counter)=HW_V;
        results.AP_max_dVdT(counter)=max(g2(lastMin:Imax));
        [~, I]=max(g3(lastMin:Imax));
		I=lastMin+I-1-1; % the extra -1 is because of a shift in points taking the derivative
        results.AP_thresh_V(counter)=dData(I); 
        results.AP_thresh_time(counter)=I;
        lastMin=Imin;
	end
    results.AP_thresh_time=results.AP_thresh_time/acqRate;
    results.AP_HW=results.AP_HW/acqRate;
	results.AP_peak_time=results.AP_peak_time/acqRate;
    results.AP_max_dVdT=results.AP_max_dVdT*acqRate;
end

