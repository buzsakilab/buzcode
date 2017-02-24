
fname = ['/home/fpbatta/Data/WarpRat/7027/analysis/2001-7-15_10-54-2/CSC7.dat'];

CR_tsd = ReadCR_tsd(fname);

fname = '/home/fpbatta/Data/WarpRat/7027/analysis/2001-7-15_10-54-2/old_analysis/SPWtimes0301.mat';

load(fname);

S_s1 = ts(S_s1.t);
E_s1 = ts(E_s1.t);
M_s1 = ts(M_s1.t);

S_s2 = ts(S_s2.t);
E_s2 = ts(E_s2.t);
M_s2 = ts(M_s2.t);


fname = '/home/fpbatta/Data/WarpRat/7027/analysis/2001-7-15_10-54-2/epoch_limits.mat';

load(fname);

cd /home/fpbatta/Data/WarpRat/7027/analysis/2001-7-15_10-54-2/

cn = List2Cell('cort_cells.list');

S = LoadSpikes(cn);

