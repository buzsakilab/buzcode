function kappa = CalcCurvature(v, vtheta, indx, dt)

%CALCCURVATURE
%   Usage:  kappa = CalcCurvature(v,theta,indx,dt)
%
% This function calculates the curvature of the path that the fly took
% while walking in the frames indexed by indx.  The inputs are the speed v
% , the angle of the velocity vector vtheta,the index of walking frames,and the time resolution
% (1/sampling rate) dt. I am convinced that this function can be
% simplified, but haven't had the time or energy to find a more efficient way.

%Written by Dan Valente
%September 2007

list = [];
phi = [];
kappa = [];
indx(end+1)=-1;
for i = 1:length(indx)-1
    if (indx(i)+1 == indx(i+1))
        list = [list indx(i)];
    else
        list = [list indx(i)];
        temp = vtheta(list)';
        if (length(temp)~= 1)
            for j = 1:length(temp)-1
                d = temp(j+1)-temp(j);
                if abs(d) >= pi
                    if (temp(j+1) < 0)
                        d = (temp(j+1)+2*pi)-temp(j);
                    elseif (temp(j+1) >= 0)
                        d = (temp(j+1)-2*pi)-temp(j);
                    end
                end
                phi(j) = d/dt;
            end
            
            phi2 = phi./v(list(1:end-1))'; %transpose for FAnalyze;
            kappa = [kappa phi2];
        end
        list = [];
        phi = [];
    end
end

return;