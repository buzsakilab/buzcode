function [Fet, nFeatures] = LoadFeatures(FileName,varargin);
% [Fet, nFeatures] = LoadFet(FileName, BufSize, CluSubset)
%
% A simple matlab function to load a .fet file

[BufSize, CluSubset] = DefaultArgs(varargin,{inf, []});

Fp = fopen(FileName, 'r');

if Fp==-1
    error(['Could not open file ' FileName]);
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



% function [fet,dim] = LoadFeatures(fbasename,fix,N)
% 
% fet = [];
% fname = [fbasename '.fet.' num2str(fix)];
% if ~N
%     fid = fopen(fname,'r');
%     N = str2num(fgetl(fid));
%     N = N-5;
%     fclose(fid);
% else
%     N = N*3;
% end
% fet_tmp = dlmread(fname,' ',1,0);    
% fet_tmp = fet_tmp(:,1:N);
% fet = [fet;fet_tmp];
% dim = N/3;