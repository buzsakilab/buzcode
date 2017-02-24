function [C, Clambda] = ChernoffGaussian2Cyl(mutheta1, muy1, S1, mutheta2, muy2, S2)
  
  draw = 0;
  
  A1 = inv(S1);
  A2 = inv(S2);
  
  
  lambda = linspace(0,1,21);
  
  
  
  lds1 = log(det(S1));
  lds2 = log(det(S2));
  
  mt1 = rem(mutheta1,2 * pi);
  mt2 = rem(mutheta2,2 * pi);
  
  delta1 = [mutheta1-mutheta2; muy1-muy2];

  md1 = delta1(1);
  md2 = - sign(md1) * (2 * pi - abs(md1));
  delta2 = [md2; delta1(2)];
  
  
  
  for i = 1:length(lambda)
    l = lambda(i);
    Gamma = l * A1 + (1-l) * A2;
    
    Clambda(i) = l * lds1 + (1-l) * lds2 + log(det(Gamma)); 
    
    Chi = A1 * inv(Gamma) * A2;
    
    M1 = delta1' * Chi * delta1;
    M2 = delta2' * Chi * delta2;
    
    M = min(M1, M2);
    
    Clambda(i) = Clambda(i) + l*(1-l)*M;
  end
  
  Clambda = Clambda / 2;
  C = max(Clambda);
  
  if draw
    figure(6)
    plot(lambda, Clambda);
    keyboard
  end
  
  
    
  