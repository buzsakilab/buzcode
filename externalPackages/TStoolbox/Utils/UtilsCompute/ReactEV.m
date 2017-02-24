function [EV, EVr] = ReactEV(x_s1, x_s2, x_m)
%   [EV EVr] = ReactEV(x_s1, x_s2, x_m)
% 
% computes reactivation EV and reverse EV for matching columns in the
% three inputs
  
  r_m_s1 = ReactCC(x_m, x_s1);
  r_m_s2 = ReactCC(x_m, x_s2);
  r_s1_s2 = ReactCC(x_s1, x_s2);
  
  EV = ( (r_m_s2 - r_m_s1 .* r_s1_s2) ./ ...
	 sqrt((1 - (r_m_s1 .^ 2)) .* (1 - (r_s1_s2 .^ 2)))) .^ 2;

  EVr = ( (r_m_s1 - r_m_s2 .* r_s1_s2) ./ ...
	  sqrt((1 - (r_m_s2 .^ 2)) .* (1 - (r_s1_s2 .^ 2)))) .^ 2;

  
  

