function h = TS2Time(t)

secs = floor(t / 10000);
d = rem(t, 10000);

mins = floor(secs/60);
secs = rem(secs, 60);

hours = floor(mins / 60);
mins = rem(mins,60);

if hours ~= 0 
  h = [num2str(hours) ':' num2str(mins, '%2d') ':' num2str(secs, '%2d') '.' num2str(d)];
else
  h = [num2str(mins, '%2d') ':' num2str(secs, '%2d') '.' num2str(d)];
end
