function [ goingUp, goingDown ] = ipFindXings( dData, level, interp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if nargin<3
        interp=0;
    end
    
    if nargin<2
        level=0;
    end
    
    isAbove=dData>level;
    goingUp=find(isAbove(2:end).*~isAbove(1:end-1));
    goingDown=find(isAbove(1:end-1).*~isAbove(2:end));
    
    if interp
        delta=dData(goingUp+1)-dData(goingUp);
        dT=(level-dData(goingUp))./delta;
        goingUp=goingUp+dT;

        delta=dData(goingDown)-dData(goingDown+1);
        dT=(dData(goingDown)-level)./delta;
        goingDown=goingDown+dT;
    end
end

