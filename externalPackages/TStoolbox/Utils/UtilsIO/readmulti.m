% reads multi-channel recording file to a matrix
% function [eeg] = function readmulti(fname,numchannel,chselect)
% last argument is optional (if omitted, it will read all the 
% channels

function [eeg] = readmulti(fname,numchannel,chselect)

if nargin == 2
  datafile = fopen(fname,'r');
  eeg = fread(datafile,[numchannel,inf],'int16');
  fclose(datafile);
  eeg = eeg';
  return
end

if nargin == 3

  % the real buffer will be buffersize * numch * 2 bytes
  % (short = 2bytes)
  
  buffersize = 4096;
  
  % get file size, and calculate the number of samples per channel
  fileinfo = dir(fname);
  N_EL = ceil(fileinfo(1).bytes / 2 / numchannel);
  
  datafile = fopen(fname,'r');
  
  mmm = sprintf('%d elements',N_EL);
%  disp(mmm);  
  
  eeg=zeros(length(chselect),N_EL);
  N_EL=0;
  numelm=0;
  while ~feof(datafile),
    [data,count] = fread(datafile,[numchannel,buffersize],'int16');
    numelm = count/numchannel;
    if numelm>0 % Kenji modified 061009.Otherwise if numelm == 0 an error occur.
        eeg(:,N_EL+1:N_EL+numelm) = data(chselect,:);
        N_EL = N_EL+numelm;
    end % Kenji modified 061009.Otherwise if numelm == 0 an error occur.
    
end
fclose(datafile);

end

eeg = eeg';
