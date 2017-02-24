function F = FischerInfo(h,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

db = median(diff(B));
F = 0;
for ii=1:length(h)
    f = (((diff(h{ii}).^2)./h{ii}(1:end-1)));
    f(isnan(f) | isinf(f))=0;
    F = F+f;
end
F = F/length(h);
end

