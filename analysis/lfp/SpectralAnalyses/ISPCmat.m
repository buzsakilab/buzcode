function [ISPCmat] = ISPCmat(sigphasegroup1,sigphasegroup2)
%[ISPCmat] = ISPCmat(sigphasegroup1,sigphasegroup2) calculates the matrix  
%of Inter-Site Phase Clustering (ISPC) between multiple signals
%
%INPUT
%   sigphasegroup1/2    [t x Nsigs] matrices of signals
%
%
%%
Nsig1 = size(sigphasegroup1,2);
Nsig2 = size(sigphasegroup2,2);
ISPCmat = zeros(Nsig1,Nsig2);

for ss1 = 1:Nsig1
    for ss2 = 1:Nsig2
        ISPCmat(ss1,ss2) = ISPC(sigphasegroup1(:,ss1),sigphasegroup2(:,ss2));
    end
end


end

