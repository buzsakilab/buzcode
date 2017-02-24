function icsi = icsi(rgx,rgy)

% icsi computes the common inter-spike interval (see Fee et al., J Neurosci. Meth., 1996)
% Adrien Peyrache, 2012

% rewritten 1/4/2017 by David Tingley 

temp = [rgx;rgy];
icsi = sort(temp);

% lx = length(rgx);
% ly = length(rgy);
% icsi = [];
% 
% [dummy,ix] = min([lx,ly]);
% if ix~=1
%     rgy2 = rgx;
%     rgx=rgy;
%     rgy = rgy2;
%     lx = length(rgx);
%     ly = length(rgy);
% end
%    
% tx=1;
% ty=1;
% if rgx(1)<rgy(1)
%     while rgx(tx)<rgy(ty) && tx<=lx
%         tx = tx+1;
%     end
%     if tx>1
%         tx=tx-1;
%     end
%     icsi = [icsi;rgy(ty)-rgx(tx)];
% end
% tx = tx+1;
% 
% while tx<lx && ty<ly
%     
%     while rgx(tx)>rgy(ty) && ty<ly
%         ty = ty+1;
%     end
%     if ty>1
%         ty=ty-1;
%     end
%     icsi = [icsi;rgx(tx)-rgy(ty)];
%     ty=ty+1;
% 
%     while rgx(tx)<rgy(ty) && tx<lx
%         tx = tx+1;
%     end
%     if tx>1
%         tx=tx-1;
%     end
%     %disp([tx ty])
%     icsi = [icsi;rgy(ty)-rgx(tx)];
%     tx=tx+1;
%     %keyboard
%     %h=waitbar(tx/lx);
%     %length(icsi)
% end
