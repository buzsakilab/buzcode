function BaslerTrackPerFrame(movname)

% 
% tspdata=load([fbasename '.tsp']);
% tspdata = tspdata(:,1);
% 
% if size(tspdata,2)>1
%   warning('tsp has data in it! Are You sure you want to continue?')
% %   keyboard
% end
% 
% nFramesTsp = size(tspdata,1);
% tspnew = [tspdata -1*ones(size(tspdata,1),6)];

readerobj  = VideoReader([movname]);
frames = read(readerobj);


whl = ApproxMedianFilter_RB_LED([movname]);
% figure(1),clf
% plot(whl(:,1),whl(:,2))
% hold on
% plot(whl(:,3),whl(:,4),'r')

nDiff = nFramesTsp - size(whl,1);
if nDiff ~=0 %in case of frame drop, interpolate the output

    twhl = (1:size(whl,1));
    ttsp = (1:size(tspdata,1));
    whl(whl==-1) = NaN;
    interpData = interp1(twhl,whl,ttsp);
    interpData(isnan(interpData)) = -1;
    tspnew(:,[2 3 6 7]) = interpData;
else
    tspnew(:,[2 3 6 7]) = whl;
end

fid = fopen([movname '.tsp'],'w');
fprintf(fid,'%i\t %f\t %f\t %f\t %f\t %f\t %f\t \n',tspnew');
fclose(fid);
%eval(['!rm ' fbasename '.ogg'])
 
      


  