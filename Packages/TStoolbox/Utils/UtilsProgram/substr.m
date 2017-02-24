function so = substr(si, s, s1)
% so = substr(si, s) substitute string 
%   
% INPUTS:
% si: the initial string  
% s: the string to be substituted
% s1: the string it will be substitued with
% OUTPUTS:
% so: the resulting string 
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  done = 0;
  
  so = '';
  ls1 = length(s1);
  ls = length(s);
  
  while ~done
    k = findstr(si, s);
    if isempty(k)
      so = [so, si];
      done = 1;
      break
    end
    so = [so, si(1:(k-1))];
    so = [so, s1];
    si = si((k+ls):end);
  end
  
  
    
    