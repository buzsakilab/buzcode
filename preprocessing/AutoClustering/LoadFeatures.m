function [fet,dim] = LoadFeatures(fbasename,fix,N)

fet = [];
fname = [fbasename '.fet.' num2str(fix)];
if ~N
    fid = fopen(fname,'r');
    N = str2num(fgetl(fid));
    N = N-5;
    fclose(fid);
else
    N = N*3;
end
fet_tmp = dlmread(fname,' ',1,0);    
fet_tmp = fet_tmp(:,1:N);
fet = [fet;fet_tmp];
dim = N/3;