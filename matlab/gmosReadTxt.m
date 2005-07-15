function [probes,X,annotation, geneName]=gmosReadTxt(file)

% GMOSREADTXT reads TXT file for the OC1 data files.

% NPPCA

% file ia a string containing the file name and the extension.
% X is a matrix with 12 columns and N( number of genes) rows   

[probes,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,annotation, geneName]=...
    textread(file,'%q %f %f %f %f %f %f %f %f %f %f %f %f %q %q',...
    'headerlines',1,'whitespace','','delimiter','\t');

X=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12];

