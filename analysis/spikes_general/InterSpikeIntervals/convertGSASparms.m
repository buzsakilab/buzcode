function [parmstype2] = convertGSASparms(parmstype1,numcells,numAS)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% If structure, returns vector (for use in fitting algorithms).
% If vector, returns structure, (needs numcells, numAS)
%Note: vectr form uses logCV... for fitting bounds and no infs. CV=0 is
%turned into 0.005
%... could just do this in the fitting function instead... set min to 0.005
zeroCV = 0.001;

if isstruct(parmstype1)
    parmstype1.GSCVs(parmstype1.GSCVs<=0) = zeroCV;
    parmstype1.ASCVs(parmstype1.ASCVs<=0) = zeroCV;
    parmstype2 = [parmstype1.GSlogrates'; log10(parmstype1.GSCVs)'; parmstype1.GSweights';...
        parmstype1.ASlogrates'; log10(parmstype1.ASCVs)'; parmstype1.ASweights(:)];
else
    parmstype2.GSlogrates = parmstype1(1:numcells)';
    parmstype2.GSCVs = 10.^parmstype1(numcells+1 : 2*numcells)';
    parmstype2.GSweights = parmstype1(2*numcells+1 : 3*numcells)';
    
    parmstype2.ASlogrates = parmstype1(3*numcells+1 : 3*numcells+numAS)';
    parmstype2.ASCVs = 10.^parmstype1(3*numcells+numAS+1 : 3*numcells+2*numAS)';
    parmstype2.ASweights = parmstype1(3*numcells+2*numAS+1 : end)';
    parmstype2.ASweights = reshape(parmstype2.ASweights,numcells,numAS);
end

end

