function [EV, EVr] = ReactRmatrices(filebase, cell_list, epochs_start, epochs_end, restrict_times)
% ReactRmatrices 
%
% computes the cell-pair correlation matrices R for sleep 1 sleep 2 and
% maze and the explained variance. The R matrices are "flattened" stored
% as column vectors (RpairC_s1, RpairC_s2, RpairC_m) corresponding to the
% correlation in cell pairs cell_I, cell_J. 
% This format is particularly convenient if data from several sessions
% are to be pooled together. 
% Inputs are the t files, and the start/end times for the epochs.
% outputs are the RpairC, cell_I, cell_J, and the EV.
% INPUTS
% filebase: stem for output fielnames  
% cell_list: name of a text file containing the names of the cells to use 
% epochs_start: vector containing the start times of sleep1, maze, sleep2
% epochs_end: vector containing the end times of sleep1, maze, sleep2
% restrict_times: cells array each element might be empty (and then
% no operation is performed) or it may contain a set of intervals, in that
% case the corresponding epoch will be restricted to the specified
% intervals. Useful to restrict maze time to actual runnign periods,
% sleep to spw periods, etc. 
  
  
  
  
  Qbin = 100 * 10; % the time bin for the Qmatrix (in 1/10000 sec)
  
  
  % loading t files: cellnames is a cell array of filenames indicating
  % which cells you want to include in the analysis, modify at will...
  
  cellnames = List2Cell(cell_list);   
  S = LoadSpikes(cellnames);
  n_cells = length(S);
  
  for i = 1:n_cells
    smin(i) = StartTime(S{i});
    smax(i) = EndTime(S{i});
  end

  
  

  
  start_s1 = epochs_start(1);
  end_s1 = epochs_end(1);
  
  start_m = epochs_start(2);
  end_m = epochs_end(2)
  
  start_s2 = epochs_start(3);
  end_s2 = epochs_end(3);

  for i = 1:length(S)
    S_s1{i} = Restrict(S{i}, start_s1, end_s1);
    
    
    S_m{i} = Restrict(S{i}, start_m, end_m);

    S_s2{i} = Restrict(S{i}, start_s2, end_s2);

  end
  
    
  
  
  
  % start calculations, the R matrices are transformed in a column
  % vector, containing only the matrix' lower triangle, and only the
  % correlations between cells not recorded on the same electrode
  % findOnSameTrodeMatrix is a function that looks at the tfiles
  % filenames and determines if they came from te same trode or not,
  % returning an output of the appropriate format for this code. 
  % the function DEPENDS on your very own convention for filenames, so go
  % look in there...
  % the outputs of this calculation are RpairC_m, RpairC_s1, RpairC_s2, (the
  % R matrices flattened in a vector, and purged of same trode stuff, and
  % cell_I, cell_J, two vectors indicating to which cell pairs the
  % correlations in the RpairCs corresponded to, useful if later you want
  % to recalculate the EVs only on subsets of cells.
  
  % notice that with the "eval" tricks you don't have to rewrite the code
  % for s1, m, s2
  
  
  sfx = { '_s1', '_m', '_s2' };
  
  [X, Y] = meshgrid(1:n_cells, 1:n_cells);
  
  idx = find((~findOnSameTrodeMatrix(cellnames)) & (X > Y) );
  
  
  
  cell_I = X(idx);
  cell_J = Y(idx);
  
  
  for i = 1:3
    sf = sfx{i};
    
    eval(['st = start' sf ';']);
    eval(['en = end' sf ';']);
    eval(['S_epoch = S' sf ';']);
    Q = MakeQfromS(S_epoch, Qbin, 'T_start', st, 'T_end', en);

    Q = tsd(Range(Q, 'ts'), full(Data(Q)));
    
    if ~isempty(restrict_times{i})
      r_st = Start(restrict_times{i}, 'ts');
      r_en = End(restrict_times{i}, 'ts');
      Q = Restrict(Q, r_st, r_en);
    end
    
    
    
    warning off
    cQ = corrcoef(Data(Q));
    warning on
    RpairC = full(cQ(idx));
    
    
    eval(['RpairC' sf ' = RpairC;']);

     
    clear Q cQ RpairC
    clear st
    clear en
    
  end
 
  
  % now compute the explained variance: ReactEV computes the partial
  % correlation coefficients, taking care of NaNs  that arise from
  % correlations from empty spike trains
%  keyboard
  [EV EVr] = ReactEV(RpairC_s1, RpairC_s2, RpairC_m);
    
  
  
  if ~isempty(filebase)
    % save everything in a file
    fname = [filebase '_' num2str(Qbin/10) '.mat'];
    save(fname, 'RpairC*', 'cell_*', 'EV*');
  end
  
  
