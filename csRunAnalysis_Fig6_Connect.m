function [ output_args ] = csRunAnalysis_Fig6_Connect( cellList )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

	csRunAnalysis_Flex(cellList, 'maxRs', 50, ...
		'maxRestSD', 6, ...
		'maxRm', 2000);
