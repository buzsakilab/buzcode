function pos = waveForm2ProbePos(meanWaveF)

% z points downward // x leftward

x = [17 -17 15 -15 13 -13 5 -5]';
z = [-70 -50 -30 -10 10 20 40 60]';
pos = zeros(length(meanWaveF),2);

for c=1:length(meanWaveF)
    w = meanWaveF{c};
    w = sum(w'.^2);
    w = w./sum(w);
    if length(w)<8
        w = [zeros(1,8-length(w)) w];
    end
    pos(c,:) = w*[x z];
end

if 0
figure(1),clf
for c=1:length(meanWaveF)
    text(-pos(c,1)/160+0.5,-pos(c,2)/160+0.5,num2str(c+1))
    hold on
end
end