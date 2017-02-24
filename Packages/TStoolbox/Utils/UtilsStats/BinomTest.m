function pout=BinomTest(s,n,p,Sided)
%function pout=myBinomTest(s,n,p,Sided)
%
% Performs a binomial test of the number of successes given a total number
% of outcomes and a probability of success. Can be one or two-sided.
%
% Inputs:
% s- (Scalar) The observed number of successful outcomes
% n- (Scalar) The total number of outcomes (successful or not)
% p- (Scalar) The proposed probability of a successful outcome
% Sided- (String) can be 'Two','Greater' or 'Lesser'. A value of
% 'Two' will perform a two-sided test, that the actual number
% of success is different from the expected number of
% successes in any direction. 'Greater' or 'Lesser' will
% perform a 1-sided test to examine if the observed number of
% successes are either significantly greater than or less than
% (respectively) the expected number of successes.
%
% Outputs:
% pout- The probability of observing the resulting value of s or a
% value more extreme (the precise meaning of which depends on
% the value of Sided) given n total outcomes with a
% probability of success of p.
%
% For example, the signtest is a special case of this where the value of p
% is equal to 0.5 (and a 'success' is defined by whether or not a given
% sample is of a particular sign.), but the binomial test and this code is
% more general allowing the value of p to be any value between 0 and 1.
%
% References:
% http://en.wikipedia.org/wiki/Binomial_test
%
% by Matthew Nelson July 21st, 2009

if nargin<4 || isempty(Sided); Sided='Two'; end
if nargin<3 || isempty(p); p=0.5; end
    
switch lower(Sided)
    case 'two'
        %note that matlab's binocdf(s,n,p) gives the prob. of getting up to AND INCLUDING s # of successes...
        E=p*n;
        if s>=E
            pout=1-binocdf(s-1,n,p); %start with the prob of getting >= s # of successes
            
            %now figure the difference from the expected value, and figure the prob of getting lower than that difference from the expected value # of successes
            dE=s-E;
            pout=pout+ binocdf(floor(E-dE),n,p); %the binonmial is a discrete dist. ... so it's value over non-integer args has no menaing... this flooring of E-dE actually doesn't affect the outcome (the result is the same if the floor was removed) but it's included here as a reminder of the discrete nature of the binomial
        else
            pout=binocdf(s,n,p); %start with the prob of getting <= s # of successes
            
            %now figure the difference from the expected value, and figure the prob of getting greater than that difference from the expected value # of successes
            dE=E-s;
            pout=pout+ 1-binocdf(ceil(E+dE)-1,n,p); %Here the ceiling is needed b/c of the -1 following it, so that integer and non-integer vals of E+dE will bothe give the correct value with the same line of code
        end
    case 'greater' %one-sided
        pout=1-binocdf(s-1,n,p); %just report the prob of getting >= s # of successes
    case 'lesser' %one-sided
        pout=binocdf(s,n,p); %just report the prob of getting <= s # of successes
    otherwise
        error(['In myBinomTest, Sided variable is: ' Sided '. Unknown sided value.'])
end