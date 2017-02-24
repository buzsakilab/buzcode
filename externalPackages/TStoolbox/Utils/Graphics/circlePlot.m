function circlePlot(M, x, y, cmap)

if length(x) ~= size(M, 1) | length(y) ~= size(M,1) | length(x) ~= length(y)
error('x and y should have the same length as dimensions of M');
end

mx = max(M);
mn = min(M);

ix = (M-mn)/(mx-mn);
colIx = ceil(63*ix+1);
colYT=mn;
for i=1:4

	colYT = [colYT mn+(mx-mn)*i/4];
end

hold on
for i=1:length(x)
	c = cmap(colIx(i),:);
	plot(x(i),y(i),'o','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerSize',4)
end

%  keyboard

c = colorbar;
set(c,'YTick',[1:16:65])
set(c,'YTickLabel',round(colYT*100)/100)
