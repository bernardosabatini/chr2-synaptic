function [ output_args ] = csRunAnalysis_Fig5_158on( cellList )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	csRunAnalysis_Flex(cellList, 'pulseStart', 500,	'maxRestSD', 10, ...
		'checkPulseStart', 2800, ...
		'checkPulseEnd', 2850);

