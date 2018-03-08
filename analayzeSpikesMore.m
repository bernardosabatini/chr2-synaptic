
DV=nan(1,85);
ML=nan(1,85);

ap200Thresh=nan(1,85);
ap200Dt=nan(1,85);
ap200HW=nan(1,85);
ap200N=nan(1,85);
ap200PulseISI=nan(1,85);

apLowThresh=nan(1,85);
apLowDT=nan(1,85);
apLowHW=nan(1,85);
apLowPulse=nan(1,85);
apLowN=nan(1,85);
apLowPulseISI=nan(1,85);

apHiThresh=nan(1,85);
apHiDT=nan(1,85);
apHiHW=nan(1,85);
apHiPulse=nan(1,85);
apHiN=nan(1,85);
apHiPulseISI=nan(1,85);

for c=[1:30 32:37 39:85]
	DV(c)=csAllCells(c).DV;
	ML(c)=csAllCells(c).ML;
	
	fff=find(csAllCells(c).pulseI==200);
	if ~isempty(fff)
		fff=fff(1);
		ap200N(c)=csAllCells(c).nAP(fff);
		if ap200N(c)>0
			ap200Thresh(c)=mean(csAllCells(c).pulseAP{fff}.AP_thresh_V);
			ap200HW(c)=mean(csAllCells(c).pulseAP{fff}.AP_HW);
			ap200Dt(c)=...
				 csAllCells(c).pulseAP{fff}.AP_peak_time(end)-csAllCells(c).pulseAP{fff}.AP_peak_time(1);
			if ap200N(c)>2
				ap200PulseISI(c)=ap200Dt(c)/(ap200N(c)-2);
			end
		end
	end
	
	fff=find(csAllCells(c).nAP>0);
	if ~isempty(fff)
		fff=fff(1);
		apLowN(c)=csAllCells(c).nAP(fff);
		if apLowN(c)>0
			apLowThresh(c)=mean(csAllCells(c).pulseAP{fff}.AP_thresh_V);
			apLowHW(c)=mean(csAllCells(c).pulseAP{fff}.AP_HW);
			apLowDT(c)=...
				 csAllCells(c).pulseAP{fff}.AP_peak_time(end)-csAllCells(c).pulseAP{fff}.AP_peak_time(1);
			apLowPulse(c)=csAllCells(c).pulseI(fff);
			if apLowN(c)>2
				apLowPulseISI(c)=apLowDT(c)/(apLowN(c)-2);
			end
		end
	end
	
	fff=find(csAllCells(c).nAP>0);
	if ~isempty(fff)
		fff=fff(end);
		apHiN(c)=csAllCells(c).nAP(fff);
		if apHiN(c)>0
			apHiThresh(c)=mean(csAllCells(c).pulseAP{fff}.AP_thresh_V);
			apHiHW(c)=mean(csAllCells(c).pulseAP{fff}.AP_HW);
			apHiDT(c)=...
				 csAllCells(c).pulseAP{fff}.AP_peak_time(end)-csAllCells(c).pulseAP{fff}.AP_peak_time(1);
			apHiPulse(c)=csAllCells(c).pulseI(fff);
			if apHiN(c)>2
				apHiPulseISI(c)=apHiDT(c)/(apHiN(c)-2);
			end
		end
	end

	
	
end