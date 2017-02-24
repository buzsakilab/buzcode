function [Ampc, Tc] = findComplexSpikes(TT, T, varargin)


minInt = 2;
maxInt = 6;

Extract_varargin;

minInt = minInt * 10;
maxInt = maxInt * 10;




t = Data(T);



%TTr = Restrict(TT, Data(T));
TTr = TT;
[A, pn] = feature_peak(TTr, [1 1 1 1]);

Tc = S2M2(Data(T), Data(T), 15, minInt, maxInt);

Ampc = S2M2(A, Data(T), 15, minInt, maxInt);
