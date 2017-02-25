function beta = CalcReorientAngle(vtheta, moves_indx)

%CALCREORIENTANGLE
%   Usage:  beta = CalcReorientAngle(vtheta, moves_indx)
%
% This function calculates the reorientation angle, or the change in angle 
% of progression after the fly stops. The inputs are the angle of the 
% velocity vector, vtheta, and an index of all the frames where the fly was
% moving, moves_indx.  

%Written by Dan Valente
%modified November 2007

m = 1;
before = 0;

for i = 1:length(moves_indx)-1
    if (moves_indx(i)+1 ~= moves_indx(i+1))
        beta(m) = 180*(vtheta(moves_indx(i))-before)/pi;
        if (beta(m) > 180)
            beta(m) = beta(m)-360;
        elseif (beta(m) < -180)
            beta(m) = beta(m)+360;
        end
        before = vtheta(moves_indx(i+1));
        m = m+1;
    end
end

return;