function m = mean(tsa,epoch)

Returns the mean  of data in tsd object

USAGE:
m = mean(tsa) 

INPUTS:
tsa: a tsd object

OUTPUTS:
m: the rate
 
% copyright (c) 2004 Francesco P. Battaglia, 
% modified 2009 by Adrien Peyrache adrien.peyrache@gmail.com
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

m = mean(tsa.data, 1);
