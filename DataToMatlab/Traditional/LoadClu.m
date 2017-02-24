% [Clu, nClusters] = LoadClu(FileName)
%
% A simple matlab function to load a .clu file

function [Clu, nClusters] = LoadClu(FileName)

Fp = fopen(FileName, 'r');

if Fp==-1
    error(['Could not open file ' FileName]);
end

nClusters = fscanf(Fp, '%d', 1);
Clu = fscanf(Fp, '%d');
nClusters = max(Clu);
fclose(Fp);