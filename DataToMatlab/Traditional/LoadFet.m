% [Fet, nFeatures] = LoadFet(FileName, BufSize, CluSubset)
%
% A simple matlab function to load a .fet file

function [Fet, nFeatures] = LoadFet(FileName,varargin);
[BufSize, CluSubset] = DefaultArgs(varargin,{inf, []});


Fp = fopen(FileName, 'r');

if Fp==-1
    error(['Could not open file ' lsFileName]);
end


nFeatures = fscanf(Fp, '%d', 1);

if isinf(BufSize)
    Fet = fscanf(Fp, '%f', [nFeatures, inf])';
else
    Fet = [];
    while ~feof(Fp)
        %TODO : add subsetting inside the buffer loop
        FetBuf = fscanf(Fp, '%f', [nFeatures, BufSize])';
        Fet = [Fet; FetBuf];
    end
end

fclose(Fp);


FileBase = FileName(1:strfind(FileName,'.fet.')-1);
El = str2num(FileName(strfind(FileName,'.fet.')+5:end));
Subset = 0;
if ~isempty(CluSubset) & FileExists([FileBase '.clu.' num2str(El)])
    Clu = LoadClu([FileBase '.clu.' num2str(El)]);
    in = ismember(Clu,CluSubset);
    Subset = 1;
    Fet = Fet(in,:);

end

