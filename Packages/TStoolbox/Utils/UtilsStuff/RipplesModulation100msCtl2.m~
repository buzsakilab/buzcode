function A = truc(A);

A = getResource(A,'MidRipS1SWS');
A = getResource(A,'MidRipS2SWS');
A = getResource(A,'SpikeData');
A = getResource(A,'HcTraceRipples');
nbC = length(S);

t1 = Range(midRipS1SWS{hcTraceRipples-4});
t2 = Range(midRipS2SWS{hcTraceRipples-4});


A = registerResource(A, 'RatioRip100_2', 'numeric', {[], []}, ...
    'ratioRip100_2', ...
    ['ratio of firing rate between ripples and control']);

A = registerResource(A, 'ProbaRip100_2', 'numeric', {[], []}, ...
    'probaRip100_2', ...
    ['probability of NULL hypothesis = {firing rate doens t change between ripples and control}']);

A = registerResource(A, 'RatioRipS1100_2', 'numeric', {[], []}, ...
    'ratioRipS1100_2', ...
    ['ratio of firing rate between ripples and control']);

A = registerResource(A, 'RatioRipS2100_2', 'numeric', {[], []}, ...
    'ratioRipS2100_2', ...
    ['ratio of firing rate between ripples and control']);

A = registerResource(A, 'ProbaRipS1100_2', 'numeric', {[], []}, ...
    'probaRipS1100_2', ...
    ['probability of NULL hypothesis = {firing rate doens t change between ripples and control}']);

A = registerResource(A, 'ProbaRipS2100_2', 'numeric', {[], []}, ...
    'probaRipS2100_2', ...
    ['probability of NULL hypothesis = {firing rate doens t change between ripples and control}']);

t1 = t1+1000;
t2 = t2+1000;
t = [t1;t2];

ripInt1 = intervalSet(t1-250,t1+250);
ripInt2 = intervalSet(t2-250,t2+250);
ripInt = intervalSet(t-250,t+250);

ripCtl11 = intervalSet(t1-1250,t1-250);
ripCtl12 = intervalSet(t1+250,t1+1250);
ripCtl21 = intervalSet(t2-1250,t2-250);
ripCtl22 = intervalSet(t2+250,t2+1250);
ripCtl1 = intervalSet(t-1250,t-250);
ripCtl2 = intervalSet(t+250,t+1250);

ratioRipS1_2 = zeros(nbC,1);
ratioRipS2_2 = zeros(nbC,1);
ratioRip_2 = zeros(nbC,1);

probaRipS1_2 = zeros(nbC,1);
probaRipS2_2 = zeros(nbC,1);
probaRip_2 = zeros(nbC,1);

percentMin = 0.05; %sfn abstract : 0.05

for i=1:nbC

	ripRate1 = Data(intervalRate2(S{i},ripInt1));
	ripRate2 = Data(intervalRate2(S{i},ripInt2));
	ripRate = Data(intervalRate2(S{i},ripInt));

	noRipRate1 = (Data(intervalRate2(S{i},ripCtl11))+Data(intervalRate2(S{i},ripCtl12)))/2;
	noRipRate2 = (Data(intervalRate2(S{i},ripCtl21))+Data(intervalRate2(S{i},ripCtl22)))/2;
	noRipRate = (Data(intervalRate2(S{i},ripCtl1))+Data(intervalRate2(S{i},ripCtl2)))/2;

	n1 = length(t1);
	n2 = length(t2);
	n = length(t);

	if (sum(ripRate1>0)>percentMin*n1) & (sum(noRipRate1>0)>percentMin*n1)

		mRip1 = mean(ripRate1);
		mCtl1 = mean(noRipRate1);
		if mCtl1==0, mCtl1=1;end
		ratioRipS1_2(i) = mRip1/mCtl1;
		
		[H,P] = ttest(ripRate1,noRipRate1);
		probaRipS1_2(i) = P;

	else
		
		ratioRipS1_2(i) = 0;
		probaRipS1_2(i) = 1;

	end

	if (sum(ripRate2>0)>percentMin*n2) & (sum(noRipRate2>0)>percentMin*n2)

		mRip2 = mean(ripRate2);
		mCtl2 = mean(noRipRate2);
		if mCtl2==0, mCtl2=1;end
		ratioRipS2_2(i) = mRip2/mCtl2;
		
		[H,P] = ttest(ripRate2,noRipRate2);
		probaRipS2_2(i) = P;

	else
		ratioRipS2_2(i) = 0;
		probaRipS2_2(i) = 1;

	end

	if (sum(ripRate>0)>percentMin*n) & (sum(noRipRate>0)>percentMin*n)

		mRip = mean(ripRate);
		mCtl = mean(noRipRate);
		if mCtl==0, mCtl=1;end
		ratioRip_2(i) = mRip/mCtl;
		
		[H,P] = ttest(ripRate,noRipRate);
		probaRip_2(i) = P;

	else
		ratioRip_2(i) = 0;
		probaRip_2(i) = 1;

	end

end

ratioRip100_2 = ratioRip_2;
probaRip100_2 = probaRip_2;
ratioRipS1100_2 = ratioRipS1_2;
probaRipS1100_2 = probaRipS1_2;
ratioRipS2100_2 = ratioRipS2_2;
probaRipS2100_2 = probaRipS2_2;



A = saveAllResources(A);


