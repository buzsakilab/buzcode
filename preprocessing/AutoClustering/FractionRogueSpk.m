function l = FractionRogueSpk(rg,tR,tC)

dt = diff(rg);
r = sum(dt<tR);

if r==0
    l=0;
else

    nS = length(rg);
    %T = rg(end)-rg(1);
    
    %better approximation of "true" f.r. (?)
    fr = 1/median(diff(rg));
    
    %a = r*T/(2*(tR-tC)*nS);
	    a = r/(2*(tR-tC)*nS*fr);
	    
%     delta = 1-4*a;
%     l = zeros(2,1);
%     l(1) = (1+sqrt(delta))/2;
%     l(2) = (1-sqrt(delta))/2;
%     ix = l>=0 & l<=1;
%     keyboard
%     if any(ix)
%         l = l(ix);
%     else
%         l=0;
%     end
%     
    
    l = a;

end
