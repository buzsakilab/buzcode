% mins = LocalMinima(x, NotCloserThan, LessThan)
%
% finds positions of all local minima in input array
% in the case that it goes down, stays level, then goes up,
% the function finds the earliest point of equality
%
% second optional argument gives minimum distance between
% two minima - if 2 are closer than this, it choses the lower between them.
%
% 3rd optional argument lets you only have minima less than a certain number
% use this option if you are computing minima of a long array - it'll take way
% less time and memory.
%
% endpoints will not be counted as minima.

% This program is the curse of my life.  Why can't things be simple?

function Mins = LocalMinima(x, NotCloserThan, LessThan)

if min(size(x))>1
    error('x should be a vector');
end
x = x(:);
nPoints = length(x);


% only look at those below LessThan
if nargin<3
    BelowThresh = (1:nPoints)';
else
    BelowThresh = find(x<LessThan);
end
xBelow = x(BelowThresh);
GapToLeft = find(diff([0; BelowThresh])>1);
GapToRight = find(diff([BelowThresh; nPoints+1])>1);

% compute left and right signs, bearing in mind that some points are missing
sDiff = sign(diff(xBelow));
LeftSign = [1; sDiff];
LeftSign(GapToLeft) = -1;
RightSign = [sDiff; -1];
RightSign(GapToRight) = 1;

% OK, now any zero right signs need to be replaced with next non-zero ...
Zeros = find(RightSign==0);
for i=fliplr(Zeros(:)')
    RightSign(i) = RightSign(i+1);
end
    
% now we can find local minima
Mins = BelowThresh(find(LeftSign<0 & RightSign>0));

% now remove any that are too close

if nargin>=2
    while 1
        TooClose = find(diff(Mins)<NotCloserThan);
        if isempty(TooClose)
            break;
        end
%        Vals = x(Mins(TooClose:TooClose+1));
        Vals = [x(Mins(TooClose)) , x(Mins(TooClose+1))];
        [dummy Offset] = max(Vals,[],2);
        Delete = TooClose + Offset -1;
        Mins(unique(Delete)) = [];
    end
end

return

if nargin<3
    LessThan = inf;
end



nPoints = length(x);
nMins = 0;
ArraySize = floor(nPoints/10);
Mins = zeros(ArraySize,1);
PrevSign = 1;
for i=1:length(x)-1
    NextSign = sign(x(i+1)-x(i));
    
    % do we have a minimum?
    if (PrevSign<0 & NextSign>0 & x(i)<LessThan)
        nMins = nMins+1;
        Mins(nMins) = i;
    end

    % reset PrevSign, if we are not in equality situation
    if NextSign
        PrevSign=NextSign;
    end
end

% use only those we have
if nMins<ArraySize
    Mins(nMins+1:ArraySize) = [];
end
    
% look for duplicates    

if nargin>=2
    while 1
        TooClose = find(diff(Mins)<NotCloserThan);
        if isempty(TooClose)
            break;
        end
        Vals = x(Mins(TooClose:TooClose+1));
        [dummy Offset] = max(Vals,[],2);
        Delete = TooClose + Offset -1;
        Mins(unique(Delete)) = [];
    end
end

return


nPoints = length(x);

% use a trick to save memory - with findstr
s = int8([1 sign(diff(x))]);


Zeros = find(s==0);
NonZeros = uint32(find(s~=0));

% wipe zeros ... 
s(Zeros) = [];

%mins = find(s(1:nPoints-1)==-1 & s(2:end)==1);
mins = double(NonZeros(findstr(s, [-1 1])));

if nargin>=3
    mins = mins(find(x(mins)<LessThan));
end

if nargin>=2
    while 1
        TooClose = find(diff(mins)<NotCloserThan);
        if isempty(TooClose)
            break;
        end
        Vals = x(mins(TooClose:TooClose+1));
        [dummy Offset] = max(Vals,[],2);
        Delete = TooClose + Offset -1;
        mins(unique(Delete)) = [];
    end
end