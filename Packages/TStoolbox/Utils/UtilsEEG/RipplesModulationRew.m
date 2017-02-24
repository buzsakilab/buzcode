function A = truc(A);

A = getResource(A,'MidRipRew');
t = Range(midRipRew{1});
A = getResource(A,'SpikeData');

nbC = length(S);


A = registerResource(A, 'RatioRipRew', 'numeric', {[], []}, ...
    'ratioRipRew', ...
    ['ratio of firing rate between ripples and control']);

A = registerResource(A, 'ProbaRipRew', 'numeric', {[], []}, ...
    'probaRipRew', ...
    ['probability of NULL hypothesis = {firing rate doens t change between ripples and control}']);



ripInt = intervalSet(t-1000,t+1000);

ripCtl1 = intervalSet(t-10000,t-1000);
ripCtl2 = intervalSet(t+1000,t+10000);

ratioRipRew = zeros(nbC,1);

probaRipRew = zeros(nbC,1);

percentMin = 0.05; %sfn abstract : 0.05

for i=1:nbC

	ripRate = Data(intervalRate2(S{i},ripInt));

	noRipRate = (Data(intervalRate2(S{i},ripCtl1))+Data(intervalRate2(S{i},ripCtl2)))/2;

	n = length(t);

	if (sum(ripRate>0)>percentMin*n) & (sum(noRipRate>0)>percentMin*n)

		mRip = mean(ripRate);
		mCtl = mean(noRipRate);
		if mCtl==0, mCtl=1;end
		ratioRipRew(i) = mRip/mCtl;
		
		[H,P] = ttest(ripRate,noRipRate);
		probaRipRew(i) = P;

	else
		ratioRipRew(i) = 0;
		probaRipRew(i) = 1;

	end

end


A = saveAllResources(A);


