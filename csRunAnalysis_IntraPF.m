function [ output_args ] = csRunAnalysis_IntraPF( cellList )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	csRunAnalysis_Flex(cellList, ...
		'pulseStart', 500,...
		'maxRs', 60);
