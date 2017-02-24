function r =  PartialCorr(t1,t2,t3);

%  USAGE
%  
%     r = PartialCorr(t1,t2,t3) computes the partial correlation r such as 
%     r = r(t1, t2 | t3). It is the squared root of the Explain Variance
%  
%  Adrien Peyrace 2007

r12 = corrcoef([t1 t2]);
r12 = r12(1,2);

r13 = corrcoef([t1 t3]);
r13 = r13(1,2);

r23 = corrcoef([t2 t3]);
r23 = r23(1,2);

r = (r12 - r13*r23)/(sqrt(1-r13^2)*sqrt(1-r23^2));
