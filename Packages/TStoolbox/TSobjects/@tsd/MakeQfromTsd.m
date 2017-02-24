function Q = MakeQfromTsd(ts,bin,varargin)

% Spike train binning
%  	
%  	USAGE:
%  	Q = MakeQfromTsd(tsd, bin,options)
% 
%  	INPUTS:
%  	tsd - a tsd, single timeseries
%  	bin 	- a bin size (in the same time_units as S).	
%  		  Bins are generated automatically and begin with first spike occurence
%  		- an intervalSet object fixing the bounds of each bin
%  		- a TS whose timestamps are the bounds of the bins 
%  		   example:	
%  			t = [0;10;20], bin = ts(t);
%  			S is binned on [0 10], [10 20].	
%  		  This last calling method is useful, for example, 
%  		  when bins have to be phase-locked to an oscillation:
%  		  bin = ts(osc_peaks)
%  		  it then becomes a 'phase' binning 
%  		  (see 'Subbins' in OPTIONS to see how to bin each cycle)
%  
%  		  Warning when bins are defined with TS:
%  		  the 'rate' option (see below) is most certainly required as the bin size is not constant.
%
%	OPTIONS:
%	'Rate'    - counts are expressed in instantaneous rate (Hz) if 1. 
%  		    Default is 0 (that is simple spike counts).
%	'Sparse'  - 1 to output a sparse matrix (useful for huge datasets), 0 otherwise (default).
%	'Subbins' - number of subbins, used only when 'bin' is a TS object.
%  		    Allows to sub-bin each interval bounded by the timestamps of the TS. Default is 1.
%  	'T_start' - StartTime for the Q matrix
%  		    defaults to min(StartTime(S))
%  		    used only when 'bin' is a numeric
%    	'T_end'   - EndTime for the Q matrix
%  		    defaults to max(EndTime(S))
%    		    used only when 'bin' is a numeric
%  	'Sparse'  - 
%  	OUTPUTS:
%  	Q - a ctsd in which the main structure is a |t| x nCells histogram of firing rates
%
%  Inspired from the previous version by David Redish (1998) and F P. Battaglia (2004).
%  These previous versions of the code used the SPARSE command 
%  to automatically generate the binning. However, this requires 
%  unnecessary computation time as it can be more efficiently 
%  coded (and runs much faster). T_start and T_end options no
%  longer exist as the tsdArray can be restricted before 
%  or the Q matrix can be itself restricted as a any other TSD.

% copyright (c) 2009 Adrien Peyrache adrien.peyrache@gmail.com
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

opt_varargin = varargin;

defined_options = dictArray({	{ 'Rate', {0, {'numeric'} } }
    				{ 'Sparse', {0, {'numeric'} } }
    				{ 'Subbins', {1, {'numeric'} } }    
    				{ 'Sparse', {0, {'numeric'} } }    

    			    });

getOpt;

nbCells = 1;
% nbCells = length(tsd.C); % Number of Cells

% in function of the type of 'bin', generates a vector of bin start time 'st' and end times 'en'

if isa(bin,'numeric')

	is = timeSpan(ts);
	T_start = Start(is);
	T_end = End(is);
	l = (T_end-T_start)/bin;
	T_end = T_start+bin*ceil(l);
	st = [T_start:bin:T_end-bin]';
	en = st+bin;

elseif	isa(bin,'ts')

	rg = Range(bin);
	st = [];

	for i=1:length(rg)-1 %now define each bin
	
		ts = rg(i);
		te = rg(i+1);
		dt = (te-ts)/Subbins;
		dt = [ts:dt:te-dt];
		st = [st;dt'];
	
	end

	en = [st(2:end);rg(end)];

elseif isa(bin,'intervalSet')
	st = Start(bin);
	en = End(bin);
end

% Creates the bin matrix
dQ = zeros(length(st),nbCells);

% for i=1:length(tsd)
    i = 1;
	spk = Range(ts);
	l = length(spk);
      
	if l
	  ix =1; %index of spike times. Begins at 1!
	  % runs through all bins
	  for j=1:length(st)
		  
		  % the idea in this piece of code is to scan all spikes 
		  % and check if it is in the given bin. as in the
		  % intervalCount method of TSD objects.
		  % but calling this method each time slows down this function.

		  while spk(ix)<st(j) & ix<l 
			  ix = ix+1;
		  end
		  while spk(ix)<=en(j) & ix<l 
			  dQ(j,i) = dQ(j,i)+1;
			  ix = ix+1;
		  end
		  if ix==l & spk(ix)<=en(j) & spk(ix)>st(j) 
			  dQ(j,i) = dQ(j,i)+1;
		  end

	  end

	end
% end

if Sparse
	dQ = sparse(dQ);
end

if Rate
	duration = (en-st)/10000;
	dQ = dQ./(duration*ones(1,size(dQ,2)));
end

Q = tsd((en+st)/2,dQ);
