%
% This is a calling routine to test & check out the power spectrum &
% spectrogram routines for unequal segment lengths. In addition, use it 
% to compare with Chronux routines when segments are of equal length. 
%
clear all;

if 0
    dir = 'G:\ravi\Chrowser\Pass~ Tioga_0e927741-9673-46e5-9050-ca1d7541bf22\';
    xfile = 'Pass~ Tioga_0e927741-9673-46e5-9050-ca1d7541bf22'
    %dir = 'G:\ravi\Chrowser\sample~ data_8ef647e3-e5ea-43a6-8c69-fb848b8db7c2\';
    %xfile = 'sample~ data_8ef647e3-e5ea-43a6-8c69-fb848b8db7c2'
else
    dir = 'Z:\xltekRawData\Wallis~ Terry_c3f44891-afa7-4fa7-a643-55c772a05241\'
    xfile = 'Wallis~ Terry_c3f44891-afa7-4fa7-a643-55c772a05241'
end

% Get header info
% Channels are labelled from C1 through C127 and '' 
% total of 128 channels
hdr = eegMex( dir, xfile);


gram = 1 ; % 0=spectra, 1=coherence
chronux = 0 ; % 0=no comparison with Chronux; 1=compare with chronux

%nSamples = 4210; 
nChannels = 2; 
nSegments = 1;
movingwin = [25, 25];

%
% Spectral Parameters
%
params.fpass = [ 0 0.5 ];
params.pad = 2;
params.err = [2 0.05];  % err(1)=0 is no err estimates; err(1)=1 is asymptotic estimates; err(1)=2 is Jacknife
params.trialave = 1;
params.Fs = 1;

%
% Tapers issues
%
halfBandWidth = 2.5; 
kCap = 2*halfBandWidth - 1;
%params.tapers = [ halfBandWidth, kCap ];
params.tapers = [ halfBandWidth, 2 ];

%
% Basic checks on inputs
%
if (gram==1) && (nChannels < 2), error( 'Coherence requires at least 2 channels...' ); end 
%if (nSegments==1) && (params.err(1)==2), error( 'Jacknife requires more than 1 segment'); end

%
% Generate segments endpoints randomly
% myrandint is a 3rd party routine (from matlab site)
%
% Randomly generated segment end points
sMarkers = reshape( sort( myrandint( 2*nSegments, 1, [ ceil(hdr.nSamples/500) : ceil(hdr.nSamples/50) ], 'noreplace' ) )', 2, nSegments )';
%sMarkers = [ ceil(hdr.nSamples/80), ceil(hdr.nSamples/65) ];

%
% Randomly select a few channels
%
if ~chronux
    %chIndices = sort( myrandint( nChannels, 1, [ round(hdr.nChans/4) : round(3*hdr.nChans/4) ], 'noreplace' ) );
    chIndices = [ 3 : 3+nChannels-1 ];
else
    %chIndices = [ 10 : 10+nChannels-1 ];
    chIndices = [ 3, 7 ];
end

%
% Randomly generate the time series
%
fulldata = eegMex( dir, xfile, chIndices, [ 1 hdr.nSamples/50 1 ] );
mDiscardBits = 0;
conversionFactor = ( 8711 / (2^21 - 0.5) ) * 2^mDiscardBits;
fulldata{:} = fulldata{:} * conversionFactor;

%
% Create a data matrix with all the segments aligned one after another
%
totalSegmentLength = sum( sMarkers(:,2) - sMarkers(:,1) + 1 );
data = zeros( totalSegmentLength, length(chIndices) ); % preallocate to ensure contiguous memory
newMarkers(1,1) = 1; 
newMarkers(1,2) = sMarkers(1,2) - sMarkers(1,1) + 1;
data( newMarkers(1,1):newMarkers(1,2), : ) = detrend( fulldata{1}( sMarkers(1,1):sMarkers(1,2), :) );
for sg = 2:size( sMarkers, 1 )
    newMarkers(sg,1) = newMarkers(sg-1,2) + 1; 
    newMarkers(sg,2) = newMarkers(sg,1) + sMarkers(sg,2) - sMarkers(sg,1);
    data( newMarkers(sg,1):newMarkers(sg,2), : ) = detrend( fulldata{1}( sMarkers(sg,1):sMarkers(sg,2), :) );
end

% To ensure that we check results from array indices beyond 1
if nChannels > 1
    ix = sort( myrandint( 1, 2, [1:length(chIndices)], 'noreplace' ) ); % Arbitrarily pick two indices from selected channels for testing results
    i1=ix(1); i2=ix(2);
    % iC = m + (n-1)*(n-2)/2, for elements of the the coherence matrix, Cmn
    iC = ix(1) + (ix(2)-1)*(ix(2)-2)/2;
else
    ix = sort( myrandint( 1, 1, [1:length(chIndices)], 'noreplace' ) ); % Arbitrarily pick 1 indices from selected channels for testing results
    i1=ix(1);
end

%
% Power spectrum/spectrogram/coherence/coherogram
%
if gram==0
    [ S, f, Serr ] = avgSpectrum( data, movingwin, params, newMarkers );
    figure; plot( f, 10*log10( S(:,i1) ), 'k', f, 10*log10( Serr(2,:,i1) ), 'g--', f, 10*log10( Serr(1,:,i1)), 'g--' ); title('Avg. Routine:: Spectrum');
    %figure; plot( f, 10*log10( S(:,i1) )); title('Avg. Routine:: Spectrum');
elseif gram==1
    [Cmn,Phimn,Smn,Smm,f,ConfC,PhiStd,Cerr] = avgCoherence( data, movingwin, params, newMarkers );
    %  C(i,j) = Cmn(:,k) where k = j + (1/2)*(i-1)*(i-2)
    figure; plot( f, Cmn(:,iC), 'k', f, Cerr(2,:,iC), 'g--', f, Cerr(1,:,iC), 'g--' ); 
    title('Avg. Routine:: Coherence'); ylim([0 1])
    %figure; plot( f, 10*log10( Cmn(:,iC) ) ); title('Avg. Routine:: Coherence-Magnitude');
    %figure; plot( f, phimn(:,iC) ); title('Avg. Routine:: Coherence-Phase');
    disp( ['Confidence level for C (confC) at (1-p) level: ' num2str( ConfC(iC)) ] );
end



%
% Use to check against Chronux: only for equal length segments
%
if chronux

    win = floor( newMarkers(1,2) / movingwin(1) );
    newMarkers(1,2) = newMarkers(1,2) - mod( newMarkers(1,2), win );
    cdata = data( [1:newMarkers(1,2)], i1 );
    cdata = detrend( reshape( cdata, [ newMarkers(1,2)/win, win ] ) );
    cdata2 = data( [1:newMarkers(1,2)], i2 );
    cdata2 = detrend( reshape( cdata2, [ newMarkers(1,2)/win, win ] ) );
    params.trialave = 1;
    if gram==0
        [ cS, cf, cSerr ] = mtspectrumc( cdata, params );
        figure; plotvector( cS, cf, 'l', cSerr );
        %figure; plot( cf, 10*log10( cS )); title('Chronux:: Spectrum');
        figure; plot( cf, 10*log10(cSerr(1,:)), cf, 10*log10(cSerr(2,:)) ); title('Chronux Error-Bar Computations');
        figure; plot( cf, 10*log10( cS ) - 10*log10( S(:,i1) )); title('Error in Spectrum = |New Routines - Chronux|');
        figure; plot( cf, 10*log10(cSerr(1,:)) - 10*log10(Serr(1,:,i1)), cf, 10*log10(cSerr(2,:)) - 10*log10(Serr(2,:,i1)) );title('Error in Error-Bar Computations = |New Routines - Chronux| ');
    elseif gram==1
        
        [cC,cphi,cS12,cS1,cS2,cf,cconfC,cphistd,cCerr]=coherencyc( cdata, cdata2, params );
        
        %figure; plotvector( cC(:,1), cf, 'n', cCerr );
        figure; plot( cf, cC(:,iC), 'k', cf, cCerr(2,:,iC), 'g--', cf, cCerr(1,:,iC), 'g--' );
        title('Chronux:: Coherence'); ylim([0 1])
        %figure; plot( cf, 10*log10( cC(:,1) ) ); title('Chronux:: Coherence-Magnitude');
        figure; plot( cf, 10*log10( cC(:,1) ) - 10*log10( Cmn(:,iC) ) ); title('Error in Coherence = |New Routines - Chronux|');
        % Phase may give a problem of 2pi difference... look into it.
        figure; plot( cf, cphi(:,1) -  Phimn(:,iC) ); title('Error in Phase = |New Routines - Chronux|');
        %
        % Note the remaining quantities do not really need to checked since
        % coherence = cross-spectrum/power spectra* power spectra, ie C = S12/(S1*S2)
        % so unlikely that S12, S1, S2 are incorrect if C is correct.
        if 1
            figure; plot( cf, 10*log10( cS1(:,1) ) - 10*log10( Smm(:,ix(1)) ) ); title('Error in Power Spectrogram-1 = |New Routines - Chronux|');
            figure; plot( cf, 10*log10( cS2(:,1) ) - 10*log10( Smm(:,ix(2)) ) ); title('Error in Power Spectrogram-2 = |New Routines - Chronux|');
        end
        %
        % Error-Bars & Confidence Levels
        disp( ['Confidence levelfor C (confC) at (1-p) level: ' num2str( cconfC) ' (Chronux)' ] );
        disp( ['Error in confidence level, confC: ' num2str( ConfC(iC) - cconfC ) ] );     
        %figure; plot( cf, cphistd(:,1), f, phistd(:,iC) ); title('Phase-Error-Bar Computations');
        figure; plot( cf, cphistd(:,1) -  PhiStd(:,iC) ); title('Error in PhiStd-1');
        figure; plot( cf, cphistd(:,1) -  PhiStd(:,iC) ); title('Error in PhiStd-2');
        figure; plot( cf, abs(cCerr(1,:,1) -  Cerr(1,:,iC)), cf, abs(cCerr(2,:,1) -  Cerr(2,:,iC))  ); title('Error in Abs(Coherence)-Error-Bar Computations = |New Routines - Chronux|');
    end
end

    



