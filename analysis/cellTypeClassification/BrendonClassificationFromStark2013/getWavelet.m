% [ wave, f, t, coh, phases, raw, coi, scale, period, scalef ] = getWavelet( x, Fs, fMin, fMax, nbins, graphics )
%
% x:        
% vector - wavelet analysis is applied
% multi-column matrix - wavelet analysis is applied to each column
% separately
%
% two-column matrix - treated as two time series and the csd, coherence is
% computed
% 3D array, 3rd dimension is 2 - treated as multiple samples from two time
% series, coherence etc is computed
%
% Fs - sampling frequency
% fMin, fMax, nbins - parameters for spectral analysis - min/max freqnecy
% (Hz) and number of frequency bins (Morlet is always used so this is
% frequency per se and not just scale)
% graphcis - flag 0/1
%
% does: compute the CWT of each segment, then if two signals also compute
% the cross-spectrum/coherence/phase lag. then average over all segments
% and plot
%
% output:
% wave - the PSD/CSD
% f, t - vectors
% coh - only for 2-channel input
% phases - either for the 1-channel or the phase-difference (from the CSD,
%       not the smoothed estimate)
% 
%
% call: wavelet, smoothwavelet, phaseplot, myjet, colormaps

% 08-oct-12 ES

% revisions
% 15-nov-12 (1) added single-channel phase; cross-spectrum, coherence, and
%               phase difference estimates
%           (2) local plotting of spectrogram + phasogram (single-channel)
%               or coherogram + phase differences (two-channel)

% to do: organize input / output handling better, plotting etc
%           also compute mean spectra/coherence/phase 
% also - external coherence computation for multiple segments (i.e. average
%       across trials)

%function [ wave, f, t, phases, coi, scale, raw ] = getWavelet( x, Fs, fMin, fMax, nbins )
function [ wave, f, t, coh, phases, raw, coi, scale, period, scalef ] = getWavelet( x, Fs, fMin, fMax, nbins, scaling, graphics )

if isa( x, 'int16' ), x = single( x ); end

if nargin < 6 || isempty( scaling )
    scaling = 'var';
end
if nargin < 7 || isempty( graphics )
    if nargout == 0
        graphics = 1;
    else
        graphics = 0;
    end
end

dt = 1 / Fs; 
s0 = 1 / fMax;
tMax = 1 / fMin;
dj = log2( tMax/s0 ) / nbins;

mother = 'MORLET';
%k0 = 6; 
%fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % scale->frequency

wave = [];
f = [];
t = [];
coh = [];
phases = [];
raw = [];
coi = [];
scale = [];
period = [];
scalef = [];

er = 0;
sx = size( x );
switch ndims( x )
    case 2
        switch min( sx )
            case 1
                nsignals = 1;
                nsegments = 1;
                nsamples = max( sx );
                x = x( : );
            case 2
                nsignals = 2;
                nsegments = 1;
                if size( x, 1 ) == 2
                    x = x';
                end
                y = x( :, 2 );
                x = x( :, 1 );
                nsamples = max( sx );
            otherwise
                nsignals = 1;
                nsegments = sx( 2 );
                nsamples = sx( 1 );
        end
    case 3
        if size( x, 3 ) ~= 2
            er = 1;
        else
            nsignals = 2;
            nsegments = sx( 2 );
            nsamples = sx( 1 );
            y = x( :, :, 2 );
            x = x( :, :, 1 );
        end
    otherwise
        er = 1;
end
if er 
    error( 'not supported' )
end

% parameters
if isa( scaling, 'char' ) && strcmp( scaling, 'var' )
    switch nsignals
        case 1
            scalef = var( x ); % can be a different number for each segment
        case 2
            scalef = std( x ) .* std( y );
    end
elseif isa( scaling, 'double' ) && size( scaling, 1 ) == ( nbins + 1 )
    scalef = scaling; % can be a different number for each frequency
    scaling = 'z';
else
    scalef = ones( 1, nsegments );
    scaling = 'none';
end
flipidx = ( nbins + 1 ) : -1 : 1;
t = ( 1 : nsamples )' / Fs;

% actually comptue
%[wave,period,scale,coi] = wavelet(Y,dt,pad,dj,s0,J1,mother,param);
switch nsignals
    case 1
        xw = zeros( nbins + 1, nsamples, nsegments );
        for i = 1 : nsegments
            [ xw( :, :, i ), period, scale, coi ] = wavelet( x( :, i ), dt, 1, dj, s0, nbins, mother );
        end
        xw = xw( flipidx, :, : ); % freq, time, segments
        wave = abs( xw ) .^ 2;
        phases = angle( xw );
        raw = xw;
    case 2
        xw = zeros( nbins + 1, nsamples, nsegments );
        yw = xw;
        coh = xw;
        for i = 1 : nsegments
            [ xw( :, :, i ), period, scale, coi ] = wavelet( x( :, i ), dt, 1, dj, s0, nbins, mother );
            [ yw( :, :, i ) ] = wavelet( y( :, i ), dt, 1, dj, s0, nbins, mother );
            % coherence (copied as is from wtc.m):
            sinv=1./(scale');
            X = xw( :, :, i );
            Y = yw( :, :, i );
            wxy = X .* conj( Y ); % complex, single trial
            sX=smoothwavelet(sinv(:,ones(1,nsamples)).*(abs(X).^2),dt,period,dj,scale);
            sY=smoothwavelet(sinv(:,ones(1,nsamples)).*(abs(Y).^2),dt,period,dj,scale);
            sWxy=smoothwavelet(sinv(:,ones(1,nsamples)).*wxy,dt,period,dj,scale);
            Rsq=abs(sWxy).^2./(sX.*sY);
            %phases( :, :, i ) = angle( sWxy );
            coh( :, :, i ) = Rsq( flipidx, : );
            %coh = abs( yo( :, 1, 2 ) .^ 2 ) ./ ( yo( :, 1, 1 ) .* yo( :, 2, 2 ) );
            
%             subplot( 4, 2, 1 ), [ c h ] = contourf( t, f, log2( abs( flipud( X ) ).^2 ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 3 ), [ c h ] = contourf( t, f, log2( abs( flipud( Y ) ).^2 ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 5 ), [ c h ] = contourf( t, f, log2( abs( flipud( wxy ) ).^2 ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 2 ), [ c h ] = contourf( t, f, log2( flipud( sX ) ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 4 ), [ c h ] = contourf( t, f, log2( flipud( sY ) ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 6 ), [ c h ] = contourf( t, f, log2( flipud( sWxy ) ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 7 ), [ c h ] = contourf( t, f, flipud( abs(wxy).^2./(X.*Y) ), 100 ); set( h, 'linestyle', 'none' );
%             subplot( 4, 2, 8 ), [ c h ] = contourf( t, f, flipud( abs(sWxy).^2./(sX.*sY) ), 100 ); set( h, 'linestyle', 'none' );
%             
%             x0 = (abs(X).^2);
%             fmat = sinv(:,ones(1,nsamples));
%             sX=smoothwavelet(fmat.*x0,dt,period,dj,scale);
            
        end
        
        % individual channels
        xw = xw( flipidx, :, : );
        yw = yw( flipidx, :, : );
        raw( :, :, :, 1 ) = xw;
        raw( :, :, :, 2 ) = yw;
        %xwave = abs( xw ) .^ 2;
        %ywave = abs( yw ) .^ 2;
        %xphases = angle( xw );
        %yphases = angle( yw );
        % cross spectrum
        xyw = xw .* conj( yw );
        phases = angle( xyw ); % from the CSD
        wave = abs( xyw );
        % for multiple trials, one can also compute the trial-averaged coherence and phase lag by:
        %cohTA = mean( abs( xyw ) .^ 2, 3 ) ./ ( mean( abs( xw ) .^ 2, 3 ) .* mean( abs( yw ) .^ 2, 3 ) );     
        %phasesTA = mod( atan2( mean( sin( phases ), 3 ), mean( cos( phases ), 3 ) ), 2 * pi );

end



%xw = flipud( xw );
f = 1 ./ period( flipidx );
coi = 1 ./ coi; % minimum freq to consider at each time point
%f = fliplr( 1 ./ period ); 
%f = fliplr( 1 ./ scale ); 
%t = ( 1 : length( x ) )' / Fs;

if graphics
    
    % here the scaling is by the signal variance
    
    figure
    
    if nsignals == 1
        nplots = 1; % PSD
    else
        nplots = 4; % [ PSD1 CSD; COH PSD2 ]
    end
    
    for np = 1 : nplots
    end
    
    % scale
    %scalef = mean( scalef );
    %scaleres = 0.25;
    scaleres = 10;
    scalename = '{\sigma}^2';
    pow = zeros( size( wave ) );
    switch scaling
        case { 'var', 'none' }
            for i = 1 : nsegments
                pow( :, :, i ) = log2( abs( wave( :, :, i ) / scalef( i ) ) );
            end
        case 'z'
            
            for i = 1 : length( f )
                pow( i, :, : ) = ( wave( i, :, : ) - scalef( i, 1 ) ) / scalef( i, 2 );
            end
    end
    pow = mean( pow, 3 ); % average over segments (1 signal)
    levels = min( pow(:) ) : scaleres : max( pow( : ) );
    mphases = mod( atan2( mean( sin( phases ), 3 ), mean( cos( phases ), 3 ) ), 2 * pi );
    %mp = []; for i = 1 : size( phases, 1 ), mp( :, i ) = circ_mean( squeeze( phases( i, :, : ) )' ); end
    
    % plot
    h1 = subplot( 1, 1, 1 ); %subplot( 2, 1, 1 );
    %[ c h ] = contourf( t, f, pow, levels );
    [ c h ] = contourf( t, f, pow, 100 );
    set( h, 'linestyle','none')
    %xlabel( 'Time (sec)' )
    ylabel( 'Frequency (Hz)' )
    title( sprintf( '%d segments, %d signals', nsegments, nsignals ) )

    % center color limits around log2(1)=0
    if strcmp( scaling, 'var' )
        clim=get(gca,'clim');
        clim=[-1 1]*max(clim(2),3);
        set(gca,'clim',clim)
    end
    
    % add the cone of influence
    line( t, coi, 'color', [ 0 0 0 ] );
    hold on
    tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
    hcoi=fill(tt,1./[period([end 1]) 1./coi period([1 end])],'w');
    %hcoi=fill(tt,[f([end 1]) coi f([1 end])],'w');
    set(hcoi,'alphadatamapping','direct','facealpha',.5)
    hold off
    
    set( h1, 'box', 'off', 'tickdir', 'out' )
    
    % add phase arrows (copied as is from xwt.m)
    if nsignals == 2
        Args.ArrowDensity = [30 30];
        Args.ArrowSize = 1;
        Args.ArrowHeadSize = 1;
        ad=mean(Args.ArrowDensity);
        Args.ArrowSize=Args.ArrowSize*30*.03/ad;
        Args.ArrowHeadSize=Args.ArrowHeadSize*Args.ArrowSize*220;
        phs_dt=round(length(t)/Args.ArrowDensity(1));
        tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/Args.ArrowDensity(2));
        pidx=fliplr( max(floor(phs_dp/2),1):phs_dp:length(period) );
        phaseplot(t(tidx),f(pidx),2*pi-mphases(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
    end
    
    % add colorbar
    h = colorbar;
    subplot( h )
    barylbls=rats(2.^(get(h,'ytick')'));
    %barylbls([1 end],:)=' ';
    barylbls(:,all(barylbls==' ',1))=[];
    set(h,'yticklabel',barylbls);
    title( scalename )
    set( h, 'box', 'off', 'tickdir', 'out' )

    colormap( h1, myjet ) 

    if 1
        figure, %h2 = subplot( 2, 1, 2 );
        switch nsignals
            case 1
                % plot phases separately
                [ c h ] = contourf( t, f, mphases, 10 );
                set( h, 'linestyle','none')
                xlabel( 'Time (sec)' )
                ylabel( 'Frequency (Hz)' )
                %colormap( h2, colormaps( myjet ) )
                set( gca, 'clim', [ 0 2*pi ] )

                
                % add the cone of influence
                lh = line( t, coi, 'color', [ 0 0 0 ] );
                hold on
                tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
                hcoi=fill(tt,1./[period([end 1]) 1./coi period([1 end])],'w');
                set(hcoi,'alphadatamapping','direct','facealpha',.5)
                hold off
                
                h = colorbar;
                subplot( h )
                title( 'Phase (rad)' )
                set( h, 'box', 'off', 'tickdir', 'out' )
                colormap( colormaps( myjet ) )

            case 2
                % also plot coherence (without phase arrows)
                % plot
                %h2 = subplot( 2, 1, 2 );
                [ c h ] = contourf( t, f, mean( coh, 3 ), 100 );
                %[ c h ] = contourf( t, f, cohTA, 100 );
                set( h, 'linestyle','none')
                xlabel( 'Time (sec)' )
                ylabel( 'Frequency (Hz)' )
                title( sprintf( '%d segments, %d signals', nsegments, nsignals ) )
                set( gca, 'clim', [ 0 1 ] )
                set( h1, 'box', 'off', 'tickdir', 'out' )                
                
                % add the cone of influence
                lh = line( t, coi, 'color', [ 0 0 0 ] );
                hold on
                tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
                hcoi=fill(tt,1./[period([end 1]) 1./coi period([1 end])],'w');
                set(hcoi,'alphadatamapping','direct','facealpha',.5)
                hold off
                
                % add phase plots (only for high coherence values inside the coi)
                aaa=2*pi-mphases;
                aaa(mean( coh, 3 )<.5)=NaN; 
                aaa( bsxfun( @lt, f' * ones( 1, nsamples ), coi ) ) = NaN;
                phaseplot(t(tidx),f(pidx),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
                
                % add colorbar
                h = colorbar;
                subplot( h )
                title( 'Coherence' )
                set( h, 'box', 'off', 'tickdir', 'out' )
                colormap( myjet )
        end
        
    end
    
    %figure, xwt( [t x( :, 1 )],[t y( :, 1 ) ], 'Pad', 1, 'Dj', dj, 'S0', s0, 'J1', nbins );
    %figure, [ Rsq,aWxy ] = waveletCoherence( [t x( :, 1 )],[t y( :, 1 ) ], 'Pad', 1, 'Dj', dj, 'S0', s0, 'J1', nbins, 'mcc', 0, 'MakeFigure', 1 ); title( 'COH' )
    %figure, [ Wxy ] = xwt( [t x( :, i )],[t y( :, i ) ], 'Pad', 1, 'Dj', log2( fMax/fMin ) / nbins, 'S0', 1/fMax, 'MaxScale', 1/fMin, 'MakeFigure', 1 ); title( 'CSD' )

end

%t = [ 0 : 1 : length( x ) - 1 ]' / Fs;
%[ rawwave, period, scale, coi sig95 ] = wt( [ t x ], 'Pad', 1, 'dj', dj, 's0', s0, 'j1', nbins, 'mother', mother, 'MakeFigure', 1 ); 

return

% notes:
% (1) it is clear what the frequencies are (nbins, log-spaced between fMin and
% fMax), but the amp is unclear to me           DONE
% (2) should use a similar approach to filter the signal at various
% frequency ranges and calculate the spiking rate frequency/phase maps

% to do:
% (1) adjust the scale properly                 DONE
% (2) get the multi-segment version working; basically use the same call to
% wavelet.m, but concatenate the segments (with intervening portions), then
% call once with the first segment to get the proper coi.       DONE
% (3) get the plotting of phases working on this diagram        DONE
% (4) 2-signal version: compute the coherence as in Torrence and Compo
% (smoothing in time- and frequency-domains each of the spectra and the
% cross-spectrum) and add the phase                             DONE

% OK now (15nov12), but not super elegant:
% should partition into two function - getWavelet and plotWavelet,
% the first should only compute, the second should plot with options:
% single signal: just power,power + phase plots.
%       e.g. [ 1 1 ] would make two plots, whereas [ 1 0 ] just on
% two signals: power/phase for each signal; csd/coh/phase for the joint
%       i.e. if only power, plots 2x2 - [ s1 s12_csd; s12_coh; s2 ]
%       if also phase, cuts each plot into two and adds the phase

getWavelet( [ x0( :, 1 ) y0( :, 1 ) ], Fs, Fmin, Fmax, nBins );
getWavelet( [ x0( :, 1 ) ], Fs, Fmin, Fmax, nBins );
getWavelet( [ x0( :, 1 : 10 ) ], Fs, Fmin, Fmax, nBins );
xy = [];
xy( :, :, 1 ) = x0( :, 1 : 10 );
xy( :, :, 2 ) = y0( :, 1 : 10 );
getWavelet( xy, Fs, Fmin, Fmax, nBins );

