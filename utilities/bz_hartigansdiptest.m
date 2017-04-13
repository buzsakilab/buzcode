function	[dip,xl,xu, ifault, gcm, lcm, mn, mj] = bz_hartigansdiptest(xpdf)
% USAGE
% [dip,xl,xu, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf)
%
% INPUTS
%
%
%
%
% OUTPUTS
%
%
% DESCRIPTION
%
% This is a direct translation by F. Mechler (August 27 2002)
% into MATLAB from the original FORTRAN code of Hartigan's Subroutine DIPTST algorithm 
% Ref: Algorithm AS 217 APPL. STATIST. (1985) Vol. 34. No.3 pg 322-32
% Appended by F. Mechler (September 2 2002) to deal with a perfectly unimodal input
% This check the original Hartigan algorithm omitted, which leads to an infinite cycle
%
% HartigansDipTest, like DIPTST, does the dip calculation for an ordered vector XPDF using
% the greatest convex minorant (gcm) and the least concave majorant (lcm),
% skipping through the data using the change points of these distributions.
% It returns the 'DIP' statistic, and 7 more optional results, which include
% the modal interval (XL,XU), ann error flag IFAULT (>0 flags an error)
% as well as the minorant and majorant fits GCM, LCM, and the corresponding support indices MN, and MJ

% sort X in increasing order in column vector
x=sort(xpdf(:));
N=length(x);
mn=zeros(size(x));
mj=zeros(size(x));
lcm=zeros(size(x));
gcm=zeros(size(x));
ifault=0;

% Check that N is positive
if (N<=0) 
   ifault=1;
   error(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
end;

% Check if N is one
if (N==1)
   xl=x(1);
   xu=x(N);
   dip=0.0;
   ifault=2;
   error(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);
end;

if (N>1)
   % Check that X is sorted
   if (x ~= sort(x))
      ifault=3;
      error(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);

   end;
   % Check for all values of X identical OR for case 1<N<4
   if ~((x(N)>x(1)) & (4<=N))
      xl=x(1);
      xu=x(N);
      dip=0.0;
      ifault=4;
      error(1,'\nHartigansDipTest.    InputError :  ifault=%d\n',ifault);

   end;
end;

% Check if X is perfectly unimodal
% Hartigan's original DIPTST algorithm did not check for this condition
% and DIPTST runs into infinite cycle for a unimodal input
% The condition that the input is unimodal is equivalent to having 
% at most 1 sign change in the second derivative of the input p.d.f.
xsign=-sign(diff(diff(x)));
% This condition check below works even 
% if the unimodal p.d.f. has its mode in the very first or last point of the input 
% because then the boolean argument is Empty Matrix, and ANY returns 1 for an Empty Matrix
posi=find(xsign>0);
negi=find(xsign<0);
if isempty(posi) | isempty(negi) | all(posi<min(negi))
   % A unimodal function is its own best unimodal approximation, with a zero corresponding dip
   xl=x(1);
   xu=x(N);
   dip=0.0;
   ifault=5;
	%fprintf(1,'\n  The input is a perfectly UNIMODAL input function\n');
   return;
end;

% LOW  contains the index of the current estimate of the lower end of the modal interval
% HIGH contains the index of the current estimate of the upper end of the modal interval
fn=N;
low=1;
high=N;
dip=1/fn;
xl=x(low);
xu=x(high);

% establish the indices over which combination is necessary for the convex minorant fit
mn(1)=1;
for j=2:N
   mn(j)=j-1;
   % here is the beginning of a while loop
   mnj=mn(j);
   mnmnj=mn(mnj);
   a=mnj-mnmnj;
   b=j-mnj;
   while ~( (mnj==1) | ((x(j)-x(mnj))*a < (x(mnj)-x(mnmnj))*b))
      mn(j)=mnmnj;
      mnj=mn(j);
      mnmnj=mn(mnj);
      a=mnj-mnmnj;
      b=j-mnj;
   end;   % here is the end of the while loop
end; % end  for j=2:N

% establish the indices over which combination is necessary for the concave majorant fit
mj(N)=N;
na=N-1;
for jk=1:na
   k=N-jk;
   mj(k)=k+1;
   % here is the beginning of a while loop
   mjk=mj(k);
   mjmjk=mj(mjk);
   a=mjk-mjmjk;
   b=k-mjk;
   while ~( (mjk==N) | ((x(k)-x(mjk))*a < (x(mjk)-x(mjmjk))*b))
      mj(k)=mjmjk;
      mjk=mj(k);
      mjmjk=mj(mjk);
      a=mjk-mjmjk;
      b=k-mjk;
   end;   % here is the end of the while loop
end; % end  for jk=1:na

itarate_flag = 1;

% start the cycling of great RECYCLE
while itarate_flag 

% collect the change points for the GCM from HIGH to LOW
% CODE BREAK POINT 40
ic=1;
gcm(1)=high;
igcm1=gcm(ic);
ic=ic+1;
gcm(ic)=mn(igcm1);
while(gcm(ic) > low)
   igcm1=gcm(ic);
   ic=ic+1;
   gcm(ic)=mn(igcm1);
end;
icx=ic;

% collect the change points for the LCM from LOW to HIGH
ic=1;
lcm(1)=low;
lcm1=lcm(ic);
ic=ic+1;
lcm(ic)=mj(lcm1);
while(lcm(ic) < high)
   lcm1=lcm(ic);
   ic=ic+1;
   lcm(ic)=mj(lcm1);
end;
icv=ic;

% ICX, IX, IG are counters for the convex minorant
% ICV, IV, IH are counters for the concave majorant
ig=icx;
ih=icv;

% find the largest distance greater than 'DIP' between the GCM and the LCM from low to high
ix=icx-1;
iv=2;
d=0.0;

% Either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50;
if ~(icx~=2 | icv~=2)
   d=1.0/fn;
else
   iterate_BP50=1;
   while iterate_BP50
		% CODE BREAK POINT 50
		igcmx=gcm(ix);
      lcmiv=lcm(iv);
      if ~(igcmx > lcmiv)
         % if the next point of either the GCM or LCM is from the LCM then calculate distance here
         % OTHERWISE, GOTO BREAK POINT 55
         lcmiv1=lcm(iv-1);
         a=lcmiv-lcmiv1;
         b=igcmx-lcmiv1-1;
         dx=(x(igcmx)-x(lcmiv1))*a/(fn*(x(lcmiv)-x(lcmiv1)))-b/fn;
         ix=ix-1;
         if(dx < d) 
            goto60 = 1; 
         else
            d=dx;
            ig=ix+1;
            ih=iv;
            goto60 = 1;
         end;
      else
         % if the next point of either the GCM or LCM is from the GCM then calculate distance here
         % CODE BREAK POINT 55
         lcmiv=lcm(iv);
         igcm=gcm(ix);
         igcm1=gcm(ix+1);
         a=lcmiv-igcm1+1;
         b=igcm-igcm1;
         dx=a/fn-((x(lcmiv)-x(igcm1))*b)/(fn*(x(igcm)-x(igcm1)));
         iv=iv+1;
         if ~(dx < d) 
            d=dx;
            ig=ix+1;
            ih=iv-1;
         end;
         goto60 = 1;
      end;
      
      if goto60
         % CODE BREAK POINT 60
         if (ix < 1) ix=1; end;
         if (iv > icv) iv=icv; end;
         iterate_BP50 = (gcm(ix) ~= lcm(iv)); 
      end;
   end; % End of WHILE iterate_BP50
end; % End of ELSE (IF ~(icx~=2 | icv~=2)) i.e., either GOTO CODE BREAK POINT 65 OR ELSE GOTO CODE BREAK POINT 50

% CODE BREAK POINT 65
itarate_flag = ~(d < dip);
if itarate_flag
% if itarate_flag is true, then continue calculations and the great iteration cycle
% if itarate_flag is NOT true, then stop calculations here, and break out of great iteration cycle to BREAK POINT 100
   
% calculate the DIPs for the corrent LOW and HIGH

% the DIP for the convex minorant
dl=0.0;
% if not true, go to CODE BREAK POINT 80
if (ig ~= icx)
   icxa=icx-1;
   for j=ig:icxa
      temp=1.0/fn;
   	jb=gcm(j+1);
      je=gcm(j);
      % if not true either, go to CODE BREAK POINT 74
      if ~(je-jb <= 1)
         if~(x(je)==x(jb))
            a=(je-jb);
            const=a/(fn*(x(je)-x(jb)));
            for jr=jb:je
               b=jr-jb+1;
               t=b/fn-(x(jr)-x(jb))*const;
               if (t>temp) temp=t; end;
            end;
         end;
      end;
      % CODE BREAK POINT 74
      if (dl < temp) dl=temp; end;
   end;
end;

% the DIP for the concave majorant
% CODE BREAK POINT 80
du=0.0;
% if not true, go to CODE BREAK POINT 90
if ~(ih==icv)
   icva=icv-1;
   for k=ih:icva
      temp=1.0/fn;
      kb=lcm(k);
      ke=lcm(k+1);
      % if not true either, go to CODE BREAK POINT 86
      if ~(ke-kb <= 1)
         if ~(x(ke)==x(kb))
            a=ke-kb;
            const=a/(fn*(x(ke)-x(kb)));
            for kr=kb:ke
               b=kr-kb-1;
               t=(x(kr)-x(kb))*const-b/fn;
               if (t>temp) temp=t; end;
            end;
         end;
      end;
      % CODE BREAK POINT 86
      if (du < temp) du=temp; end;
   end;
end;

% determine the current maximum
% CODE BREAK POINT 90
dipnew=dl;
if (du > dl) dipnew=du; end;
if (dip < dipnew) dip=dipnew; end;
low=gcm(ig);
high=lcm(ih);      

end; % end of IF(itarate_flag) CODE from BREAK POINT 65

% return to CODE BREAK POINT 40 or break out of great RECYCLE;
end; % end of WHILE of great RECYCLE

% CODE BREAK POINT 100
dip=0.5*dip;
xl=x(low);
xu=x(high);


