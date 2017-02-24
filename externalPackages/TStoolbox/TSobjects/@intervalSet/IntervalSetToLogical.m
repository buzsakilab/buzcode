function boolean = IntervalsToLogical(is)
% Converts intervalset data into a long logical array where periods inside
% any interval are 1 and outside interval are 0.  


b = FirstTime(is);
e = LastTime(is);

boolean = b:e;

for a = 1:length(length(is));
    boolean(Start(subset(is,a)):End(subset(is,a))) = 1;
end