function testspecgram(data)


%cd 'C:\Documents and Settings\Admin\Desktop\';
%data=wavread('bird109_26519_on_Aug_19_16_33.wav');

params.tapers = [3 5];
params.fpass = [100 20000];
params.Fs = 44100;
params.pad = 2;
max_time =10; % seconds per run
max_tapers = 150;
increment = 1.5;
profile on


nsamples = 1000;
if 1
slow_results = [];
    
while 1
    tic
    [S,t,f] = mtspecgramc_slow( data(1:nsamples), [0.01 0.001], params );
    time = toc;
    result = [nsamples time];
    fprintf( 'ran %d samples in %d seconds\n',nsamples, time );
    slow_results = [slow_results ;result];
    if time > max_time 
        break
    end
    nsamples = round(nsamples * increment);
end
slow_results
fig=figure();
ax=axes('XScale','log','YScale','log');
axes(ax);
h=line( 'Xdata',slow_results(:,1),'Ydata',slow_results(:,2),'Marker','*');
xlabel('number of samples')
ylabel('time');
title('Original mtspecgramc');
grid on
drawnow;
saveas(fig,'datalength_slow.png');

nsamples = 1000;
fast_results = [];
while 1
    tic
    [S,t,f] = mtspecgramc( data(1:nsamples), [0.01 0.001], params );
    time = toc;
    result = [nsamples time];
    fprintf( 'ran %d samples in %d seconds\n',nsamples, time );
    fast_results = [fast_results ;result];
    if time > max_time 
        break
    end
    nsamples = round(nsamples * increment);
end
fast_results
fig=figure();
ax=axes('XScale','log','YScale','log');
axes(ax);
h=line( 'Xdata',fast_results(:,1),'Ydata',fast_results(:,2),'Marker','*');
xlabel('number of samples')
ylabel('time');
title('Modified mtspecgramc - preallocate space');
grid on
drawnow;
saveas(fig,'datalength_fast.png');

compare = [];
n = 1;
while n <= min(length(slow_results(:,1)),length(fast_results(:,1)))
    compare_one = [slow_results(n,1) slow_results(n,2)/fast_results(n,2)];
    compare = [compare ;compare_one];
    n = n + 1; 
end
compare
fig=figure();
ax=axes('XScale','log','YScale','lin');
axes(ax);
h=line( 'Xdata',compare(:,1),'Ydata',compare(:,2),'Marker','*');
title('Preallocation slowdown/speedup of mtspecgramc');
xlabel('number of samples')
ylabel('speedup');
grid on
drawnow;
saveas(fig,'speedup.png');


end;

nsamples=10000;
results = [];
n = 1;
while 1
    tic
    params.tapers = [n (2*n-1)];
    [S,t,f] = mtspecgramc( data(1:nsamples), [0.01 0.001], params );
    time = toc;
    result = [params.tapers(2) time];
    fprintf( 'ran %d samples in %d seconds with tapers %d %d\n',nsamples, time,params.tapers(1),params.tapers(2) );
   results = [results ;result];
    if time > max_time || params.tapers(2) > max_tapers
        break
    end
    n = round(n * increment);
end
fig=figure();
ax=axes('XScale','log','YScale','log');
axes(ax);
h=line( 'Xdata',results(:,1),'Ydata',results(:,2),'Marker','*');
xlabel('tapers')
ylabel('time');

drawnow;
saveas(fig,'tapers.png');


stats = profile('info')
