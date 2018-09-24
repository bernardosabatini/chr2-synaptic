function [ meanData, stdData ] = processWindowData
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	ll=get(gca, 'Children');
	dCount=length(ll);
	for counter=1:dCount
		newData=ll(counter).YData;
		if counter==1
			meanData=newData/dCount;
		else
			meanData=meanData+newData/dCount;
		end
	end
	stdData=meanData.*0;
	for counter=1:dCount
		newData=ll(counter).YData;
		stdData=stdData+(newData-meanData).^2/dCount;
	end
	stdData=sqrt(stdData);
	
	
end

