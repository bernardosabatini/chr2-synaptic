function [ output_args ] = csLoadExcel( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    evalin('base', 'global csTableNum csTableTxt csTableRaw csTableSize')
    global csTableNum csTableTxt csTableRaw csTableSize

 	clear('csTableNum', 'csTableTxt', 'csTableRaw', 'csTableSize');
	evalin('base', 'global csTableNum csTableTxt csTableRaw csTableSize')
    global csTableNum csTableTxt csTableRaw csTableSize
   
    [filename,pathname,~] = uigetfile('*.*','Select excel file with annotations');
    if isequal(filename,0)
        disp('User selected Cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
        [csTableNum, csTableTxt, csTableRaw] = xlsread(fullfile(pathname, filename));
    end
    
    aa=cellfun(@isnan, csTableRaw(1,:), 'un',0);
    lastCol=find(cell2mat(cellfun(@all, aa, 'un',0)));
    if isempty(lastCol)
        lastCol=size(csTableRaw,2);
    else
        lastCol=lastCol(1)-1;
    end
    
    aa=cellfun(@isnan, csTableRaw(:,1), 'un',0);
    lastRow=find(cell2mat(cellfun(@all, aa, 'un',0)));
    if isempty(lastRow)
        lastRow=size(csTableRaw,1);
    else
        lastRow=lastRow(1)-1;
    end   
    disp([num2str(lastCol) ' columns of data loaded for ' ...
        num2str(lastRow-1) ' cells'])
    csTableSize=[lastCol, lastRow];

    
    
    
    

