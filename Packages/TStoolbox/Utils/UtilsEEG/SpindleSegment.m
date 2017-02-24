function [t_start, t_end, gaps] = SpindleSegment(SQ, HighThresh, LowThresh)
  
% [t_start, t_end] = SpindleSegment(SQ, HighThresh, LowThresh)
%
% Detects spindles and returns the times of the spindles through
% INPUTS:
% SQ: the "sum of Q" tsd
% HighThresh, LowThresh: thresholds for spindle detection 
% OUTPUTS: 
% t_start and t_end:   the throughs times, as beginnings and ends of
% spindle cycles 
% gaps: the indices corresponding to the gaps between spindles  

% fpbatta 2003
  
  
  
  
  do_plot = 0;
  [S, E, M] = findSpindles(SQ, HighThresh, LowThresh);
  
% each spindle oscillation (trough to trough) is a time bin

  s = Data(SQ);
  ds = diff(s);
  ds1 = [0 ; ds(:)];
  ds2 = [ds(:) ; 0];
  p = find(ds1 < 0 & ds2 >= 0);
  t = Range(SQ, 'ts');
  t = t(p);
  
  t_close_idx = find(diff(t) < 500);
  t_close_idx = t_close_idx(:);
  mm = s(t_close_idx) < s(t_close_idx+1);
  tz = t_close_idx + mm;
  
  
  
  t(t_close_idx) = t(tz);
  
  t_close_idx = t_close_idx + 1;
  
  tx = ones(size(t));
  tx(t_close_idx) = 0;
  
  tx = find(tx);
  t = t(tx);
  
  
  
  
  
  t_start = [];
  t_end = [];
  t_tot = [];
  
  start = Data(S);
  stop = Data(E);
  gaps = [1];
  
  for i = 1:length(Data(S))
    t1_ix(i) = min(find(t > start(i)));
    t2_ix(i) = max(find(t < stop(i)));

    pt = t(t1_ix(i):t2_ix(i));
    pt = pt(:);
    pt1 = pt(1:end-1);
    pt1 = pt1(:);
    pt2 = pt(2:end);
    pt2 = pt2(:);
    
    t_start = [t_start ; pt1];
    t_end = [t_end ; pt2];
    t_tot = [t_tot; pt];
    gaps(end+1) = gaps(end) + length(pt1);
    
  end
  
  gp1 = gaps(1:end-1);
  gp1 = gp1(:);
  gp2 = gaps(2:end);
  gp2 = gp2(:);
  
  gaps = [gp1 gp2];
  
  
  t = t_tot;

  if do_plot
    figure(1)
    clf
    PlotEEG(SQ)
    hold on
    plot(Data(S)/10000, ones(size(Data(S))), 'go'), plot(Data(E)/10000, ...
						  ones(size(Data(S))), 'ro')
    plot(t/10000, zeros(size(t)), 'rx')
    keyboard;
  end
  