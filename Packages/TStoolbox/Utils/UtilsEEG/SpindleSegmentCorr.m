function [cc, sp_sparsity, Qtot, Qbound, S, E, M] = SpindleSegmentCorr(SQ, Sp)
  
  
  
  do_plot = 0;
  [S, E, M] = findSpindles(SQ, length(Sp)/36, length(Sp)/72);
  
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
  
  clear tx
  
  start = Data(S);
  stop = Data(E);
  
  if do_plot
    figure(1)
    clf
    PlotEEG(SQ)
    hold on
    plot(Data(S)/10000, ones(size(Data(S))), 'go'), plot(Data(E)/10000, ones(size(Data(S))), 'ro')
  end
  
  for i = 1:length(Data(S))
    t1_ix(i) = min(find(t > start(i)));
    t2_ix(i) = max(find(t < stop(i)));
    
    pt = t(t1_ix(i):t2_ix(i));
  
    
    if(do_plot)
      
      plot(pt/10000, zeros(size(pt)), 'rx')
    end
    
  end
  
  if do_plot
    keyboard;
  end
  
  if length(Data(S)) > 0
  
    max_spindle = max(t2_ix-t1_ix);
  else
    max_spindle = 1;
  end
  
    
  sp_corr = zeros(length(Data(S)), max_spindle);
  sp_N = zeros(1, max_spindle);
  sp_sparsity = zeros(1, length(Data(S))* max_spindle);
  
  Qtot = [];
  Qbound = 0;
  n = 1;
  for i = 1:length(Data(S))
    
    pt = t(t1_ix(i):t2_ix(i));
    pt1 = pt(1:end-1);
    pt2 = pt(2:end);
    
    qs = zeros(length(Sp), length(pt1));
    
    for j = 1:length(Sp)
      h = Histogram_intervals(Data(Sp{j}), pt1, pt2);
      qs(j,:) = (h(:))';
    end
    
    qs = qs * diag(1./(pt2-pt1));
    
    Qtot = [Qtot qs];
    Qbound = [Qbound (Qbound(end)+length(pt))];
    
    for k = 1:length(pt1)
      sp_sparsity(n) = sum(qs(:,k) > 0) / length(Sp);     
      n = n+1;
    end
    
    cq = corrcoef(qs);
%    keyboard
    
    for j = 1:length(pt1)
      sp_corr(i,j) = sum(diag(cq, j-1));
      sp_N(j) = sp_N(j)+ length(diag(cq, j-1));
    end
    
  end
  
  cc = sum(sp_corr) ./ sp_N;
  sp_sparsity = sp_sparsity(1:n-1);
    