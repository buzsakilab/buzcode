function C = corrcoefsm(S,sd)


nv = length(S);
C = zeros(nv);
SCov = zeros(nv,1);

%  T_start_init = inf;
%  T_end_init = -inf;
%  
%  for iC = 1:length(S)
%     if ~isempty(Data(S{iC})) %class(ts)
%        T_start_init = min(T_start_init, StartTime(S{iC})); %class(ts)
%        T_end_init = max(T_end_init, EndTime(S{iC})); %class(ts)
%     end
%  end

%  T = T_end_init - T_start_init;

for x=1:nv

	ST1 = Range(S{x});
	N1 = length(ST1);
	CovG = 0;
	
	lbound=1;

	for i = 1:N1
		while (ST1(lbound) < ST1(i) - 6*sd) & (lbound<N1)
			lbound = lbound+1;
		end

		ix=0;

		while (ST1(lbound+ix) < ST1(i) + 6*sd) & (lbound+ix<N1)
			CovG = CovG + exp((ST1(i)-ST1(lbound+ix))^2 /(-4*sd^2));
			ix = ix+1;
		end
			
	end

	T = ST1(end)-ST1(1);
	CovG = CovG/(2*sd*sqrt(pi)*T) - N1*N1/T^2;
	SCov(x) = CovG;


end


for x=1:nv

	ST1 = Range(S{x});
	N1 = length(ST1);

	for y=x:nv

		if x==y
			C(x,y) = 1;
			C(y,x) = 1;
		else
			
			ST2 = Range(S{y});
			N2 = length(ST2);

			% Calculate convolution integrals:
	
			CovG = 0; % Gaussian-smoothed covariances

			lbound = 1;
			i=1;
			T_start = 0;	
			T_end = 0;

			while i<N1 & lbound<N2

				while (ST2(lbound) < ST1(i) - 6*sd) & (lbound<N2)
					lbound = lbound+1;
				end
				if T_start==0
					T_start=ST2(lbound)
				end

				j=0;
	
				while (ST2(lbound+j) < ST1(i) + 6*sd) & (lbound+j<N2)
					CovG = CovG + exp((ST1(i)-ST2(lbound+j))^2 /(-4*sd^2));
					j = j+1;
				end
				T_end = ST2(lbound+j);
			
				i=i+1;
				display(j)

			end

			T = T_end-T_start;	

			% Normalize integral and subtract firing rate term:
			CovG = CovG/(2*sd*sqrt(pi)*T) - N1*N2/T^2;
			CovG = CovG/sqrt(SCov(x)*SCov(y));
			C(x,y) = CovG;
			C(y,x) = CovG;
		end

	end
end