%function [ptab,map] = MultimodalTest(data, max_k, Verbose)
%Silverman's method for testing succesive hypotheses
%that a dataset has 1, 2, 3,..., max_k modes.
%
% from A Sirota toolbox

function [ptab,htab] = MultimodalTest(data,max_k)

hrange = [0.001:0.001:1];
boot_iterations=200;
%max_k = 1;

m1=max(data);	%rescale based on the data.
m2=min(data);
modect.scale=(m1-m2)/200;

thr= prctile(data,[0.1 99.9]);
modect.min=thr(1);%m2-(m1-m2)/10;
modect.max=thr(2); %m1+(m1-m2)/10;

datavar = var(data);

if 0 %very slow , goes through all the h in interval 0 1
    for i=1:length(hrange)
        h = hrange(i);
        nmodes = countmodes(data,h,modect);
        map(i,:)=[h nmodes];
        %if nmodes==max_k
        %    break
        %end
    end
    for i=1:max_k;
        htab(i,:)=[i min(map(map(:,2)==i,1))];
    end
    

else
    hlow=0; hhigh=datavar; hdiff=1e-10;map=[];
    for i=1:max_k;
        %looking for the h that makes i modes - transition to i+1
        while abs(hhigh-hlow)>hdiff
            hcur = 0.5*(hlow+hhigh);
            nmodes= countmodes(data,hcur,modect);
%            nmodes_high = countmodes(data,hhigh,modect,Verbose);
            if nmodes>i
                hlow=hcur;
            else
                hhigh=hcur;
            end
        end

        htab(i,:)= [i hhigh];
    end
end

for i=1:max_k
    ptab(i)= boot(data,htab(i,2),htab(i,1),boot_iterations,modect);
    %fprintf('for null of <%d modes p=%f\n',i,ptab(i));
end
%keyboard
return


%helper functions

function p = boot(data,h0,modect_target,boot_iterations,modect)
p=0;
for i=1:boot_iterations
    if(countmodes(boot_draw(data,h0),h0,modect)<=modect_target)
        p = p+1;
    end
end
p = 1-p/boot_iterations;
return


function [nmodes,modes]=countmodes(in,h,modect)
x= modect.min:modect.scale:modect.max;
nmodes=0;
%variance=var(in);
len=length(x);
nin = length(in);
if 1
    estim = zeros(len,1);
    for i = 1:len
        estim(i)= sum(pdf('Normal',(x(i)-in(:))/(h),0,1)); %put the variance in h calculation
    end
    estim = estim/nin/h;
else
    estim = hist(in,x);
end


modes = LocalMinima(-estim,1,max(-estim));
nmodes = length(modes);

if 0
    figure(763); clf
    plot(x,estim);
    hold on
    Lines(x(modes),[],'r');
    title([num2str(nmodes) ' Modes ']);
    drawnow;
%   waitforbuttonpress;
end;
return

function d = boot_draw(data,h0)
ct = length(data);
d=data(randsample(ct,ct));
d = sqrt(1+h0.^2/var(d)).*(d+h0.*randn(ct,1));
return
