function [C, B, C2] = CtsdPETH(d, e, tmax)

  t = Range(d, 'ts');
  step = median(diff(t))/10;
  
  sz = ceil(tmax/step);
  
  if(min(size(e)) ~= 1)
    error('e must be row or column');
  end
  
  ev = findAlignment(d, e);
  
  if(size(ev, 2) ~= 1)
    ev = ev';
  end
  
  B = -sz:sz;
  dd = Data(d);
  if (size(dd, 1) ~= 1)
    dd = dd';
  end
  dd = [dd 0];  
  evtot = ev;
  C= zeros(size(B));
  C2= zeros(size(B));
  
  
  dstep = 100;
  
  for i = 1:dstep:length(evtot)
    ev = evtot(i:min(length(evtot), i+dstep-1));
    tms = repmat(ev, 1, length(B))+repmat(B, length(ev), 1);
    t1 = length(t);
    tms(find(tms < 1)) = t1+1;
    tms(find(tms > t1)) = t1 + 1;
    C = C + sum(dd(tms));
    C2 = C2 + sum(dd(tms).^2);
  end
  
  
  lt = length(evtot);
  C = C / lt;
  C2 = C2/(lt-1) - (lt/lt-1)*C.*C;
  
  C2 = sqrt(C2);

  B = B * step;
  
  
  
  
  