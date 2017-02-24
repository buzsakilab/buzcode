% Dependency of LoadPosition.
% [CleanWhl GoodRanges] = CleanWhlFR(Whl, StretchLen, JumpSize,Gap);
%
% "cleans up" a wheel file by interpolating missing stretches
% up to StretchLen long (default 20), for which the endpoints 
% don't differ by more than JumpSize, (default 30).
%
% also returns the ranges where the whl file is valid (in .whl units)
% GoodRanges which gives start and end samples of the good ranges
% (so Whl(GoodRanges) has no -1 values).
%
% if there are any very high derivatives left over, it warns you

% This is a modified version of CleanWhl. If the discrimination of Front and 
% Rear LED is imcomplete in MakeWhlFile_3spots program, uncleaned Whl file has 
% -1 -1 in only either whl(:,1:2) or whl(:,3:4). If those Whl file are used in Original 
% version of CleanWhl, the resultant cWhl file use -1 in whl(:,3:4) for interpolation.
% To prevent this problem, the value of whl(:,3:4) is also taken into account for the rang 
% of interpolation.

% Before useing interp1, remove the data which has (-1,-1) or big Gap between two continous rows.
% Assing the (-Gap, -Gap) to the row which has big Gap or (-1,-1).
% interporate.
% if the Gap is longer than StretchLen in terms of number of .whl rows or longer than JumpSize in terms of distance, remove the interporated values.
%
%
% change by A Peyrache, @2012

function [cWhl, GoodRanges_F] = CleanWhlForR(Whl, StretchLen, JumpSize, Gap)

% If the gap between the good strech is more than StrethcLen in terms of Whl row number,remove interporated values.
if nargin<2
	StretchLen = 30;
end

% If the Gap between the good strech is more than JumpSize, remove interporated values.
if nargin<3
    JumpSize = 30;
end

% if the distance between the two contimous rows are more than Gap centimeter, It's a big jump and do not use as an input for inpterp1.
if nargin<4,
	Gap = 30;
end

nWhl = size(Whl,1);

% interpolate missing values or large jumps.
% the value of whl(:,3:4) is also taken into account for the rang of interpolation.
% A transision to and form (-1,-1) should be taken as a BigJump.

% I hsould use distance, not the one dimentinal projection of trajectory, by the way.

whltemp = Whl;
whltemp(find(whltemp)==-1) = -Gap;
dist_F = sqrt(diff(whltemp(:,1)).^2+diff(whltemp(:,2)).^2);
dist_R = sqrt(diff(whltemp(:,3)).^2+diff(whltemp(:,4)).^2);
BigJump_F = dist_F>Gap;
BigJump_R = dist_R>Gap;

Good_F = find(Whl(:,1)>-1 & ~([BigJump_F;0] | [0;BigJump_F]));
Bad_F = find(~(Whl(:,1)>-1 & ~([BigJump_F;0] | [0;BigJump_F])));
Good_R = find(Whl(:,3)>-1 & ~([BigJump_R;0] | [0;BigJump_R]));
Bad_R = find(~(Whl(:,3)>-1 & ~([BigJump_R;0] | [0;BigJump_R])));

whltemp(Bad_F,1:2) = -Gap;
whltemp(Bad_R,3:4) = -Gap;

WhlNaN = Whl;
WhlNaN(find(Whl==-1)) = NaN;

% Give -1 outside of the interpolation.

if length(Good_F)<2 || length(Good_R)<2;
    cWhl(:,1:2) = -ones(size(Whl,1),2);
else
	cWhl(:,1:2) = interp1(Good_F, Whl(Good_F,1:2), 1:nWhl, 'linear', -1);
	cWhl(:,3:4) = interp1(Good_R, Whl(Good_R,3:4), 1:nWhl, 'linear', -1);
end


% find missing stretches for Front LED
dGoodF = [-(whltemp(1,1)==-Gap) ; diff(whltemp(:,1)>-Gap)];
BadStartF = find(dGoodF<0);
BadEndF = find(dGoodF>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,1)==-1
	BadEndF = [BadEndF; nWhl];
end

if length(BadStartF)>length(BadEndF)
	BadEndF = [BadEndF; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStartF>0)

    StartIndF = clip(BadStartF-1, 1, nWhl); % StartInd and EndInd give the 
    EndIndF = clip(BadEndF+1, 1, nWhl);     % points you are interpolating between
    
	dist_F = sqrt((Whl(StartIndF,1)-Whl(EndIndF,1)).^2+(Whl(StartIndF,2)-Whl(EndIndF,2)).^2);
    ToChuckF = find(BadEndF-BadStartF>=StretchLen ...
		| dist_F > JumpSize);
	% chuck em

	for i=ToChuckF(:)'
       	cWhl(BadStartF(i):BadEndF(i),1:2) = NaN;
	end
end

% find missing stretches for Rear LED
dGoodR = [-(whltemp(1,3)==-Gap) ; diff(whltemp(:,3)>-Gap)];
BadStartR = find(dGoodR<0);
BadEndR = find(dGoodR>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,3)==-1
	BadEndR = [BadEndR; nWhl];
end

if length(BadStartR)>length(BadEndR)
	BadEndR = [BadEndR; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStartR>0)
    StartIndR = clip(BadStartR-1, 1, nWhl); % StartInd and EndInd give the 
    EndIndR = clip(BadEndR+1, 1, nWhl);     % points you are interpolating between

  	dist_R = sqrt((Whl(StartIndR,3)-Whl(EndIndR,3)).^2+(Whl(StartIndR,4)-Whl(EndIndR,4)).^2);
	ToChuckR = find(BadEndR-BadStartR>=StretchLen ...
		| dist_R > JumpSize);
	
	% chuck em
	for i=ToChuckR(:)'
       	cWhl(BadStartR(i):BadEndR(i),3:4) = NaN;
	end
end


if 0 % OLD VERSION (BUG?)
% % now find good ranges
% dcGood = [-(Whl(1,1)==1) ; diff(cWhl(:,1)>-1)];
% GoodStart = find(dcGood>0);
% GoodEnd = find(dcGood<0)-1;
% % if last point is good, need to finish GoodEnd
% if cWhl(end,1)>-1
%     GoodEnd = [GoodEnd; nWhl];
% end
% GoodRanges = [GoodStart, GoodEnd];
else
    dcGood_F = diff([0; cWhl(:,1)>-1; 0]);
    GoodStart_F = find(dcGood_F>0);
    GoodEnd_F = find(dcGood_F<0)-1;
    GoodRanges_F = [GoodStart_F, GoodEnd_F];
end


% delete singletons
%%% I don't think that I need to remove singletons here.(Kenji 060405)
%%if length(GoodStart>0)
%   Singletons = find(GoodStart==GoodEnd);
%    cWhl(GoodStart(Singletons),:) = -1;
%    GoodRanges(Singletons,:) = [];
% end

%keyboard


return
