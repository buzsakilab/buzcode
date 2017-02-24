%
% This is a calling routine to test & check out the power spectrum &
% spectrogram routines for unequal segment lengths. In addition, use it 
% to compare with Chronux routines when segments are of equal length. 
%
clear all;

gram = 0 ; % 0=spectra, 1=coherence
chronux = 0; % 1=compare, 0=don't
line = 0;

nSamples = 4210; 
nChannels = 1; 
nSegments = 1;
movingwin = [50, 50];

%
% Spectral Parameters
%
params.fpass = [ 0 0.5 ];
params.pad = 2;
params.err = [1 0.85];  % err(1)=0 is no err estimates; err(1)=1 is asymptotic estimates; err(1)=2 is Jacknife
params.trialave = 1;
params.Fs = 1;

%
% Tapers issues
%
halfBandWidth = 2.5; 
kCap = 2*halfBandWidth - 1;
%params.tapers = [ halfBandWidth, kCap ];
params.tapers = [ halfBandWidth, 1 ];

%
% Basic checks on inputs
%
if (gram==1) && (nChannels < 2), error( 'Coherence requires at least 2 channels...' ); end 
if (nSegments==1) && (params.err(1)==2), error( 'Jacknife requires more than 1 segment'); end

%
% Generate segments endpoints randomly
% myrandint is a 3rd party routine (from matlab site)
%
if ~chronux
    % Randomly generated segment end points
    sMarkers = reshape( sort( myrandint( 2*nSegments, 1, [ 1 : nSamples ], 'noreplace' ) )', 2, nSegments )';
else
    % Equal length segments (to compare with Chronux)
    sMarkers = ones( nSegments, 2 ); 
    sMarkers( :, 2 ) = round( nSamples/2 ); 
end

%
% Randomly generate the time series
%
fulldata = randn( nSamples, nChannels );  
if line % add line harmonic for testing purposes 
    f1 = 0.45; a1 = 0.20;
    f2 = 0.25; a2 = 0.15;
    for c=1:size(fulldata,2)
        mx = max(fulldata(:,c));
        fulldata(:,c) = fulldata(:,c) + a1*mx*sin( f1*2*pi*[1:size(fulldata,1)]') + a2*mx*sin( f2*2*pi*[1:size(fulldata,1)]') ;
    end
end

%
% Randomly select a few channels
%
chIndices = sort( myrandint( ceil(nChannels/1.5), 1, [ 1 : nChannels ], 'noreplace' ) );
%chIndices = [ 1 : nChannels ]; 

%
% Create a data matrix with all the segments aligned one after another
%
totalSegmentLength = sum( sMarkers(:,2) - sMarkers(:,1) + 1 );
data = zeros( totalSegmentLength, length(chIndices) ); % preallocate to ensure contiguous memory
newMarkers(1,1) = 1; 
newMarkers(1,2) = sMarkers(1,2) - sMarkers(1,1) + 1;
data( newMarkers(1,1):newMarkers(1,2), : ) = fulldata( sMarkers(1,1):sMarkers(1,2), chIndices(:));
for sg = 2:size( sMarkers, 1 )
    newMarkers(sg,1) = newMarkers(sg-1,2) + 1; 
    newMarkers(sg,2) = newMarkers(sg,1) + sMarkers(sg,2) - sMarkers(sg,1);
    data( newMarkers(sg,1):newMarkers(sg,2), : ) = fulldata( sMarkers(sg,1):sMarkers(sg,2), chIndices(:));
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
    [ S, f, Serr ] = ueSpectrogram( data, movingwin, params, newMarkers );
    figure; plotvector( S(:,i1), f, 'l', squeeze(Serr(:,:,i1)) );
    %figure; plot( f, 10*log10( S(:,i1) )); title('Avg. Routine:: Spectrum');
elseif gram==1
    [Cmn,Phimn,Smn,Smm,f,I,ConfC,PhiStd,Cerr] = ueCoherence( data, params, newMarkers );
    %  C(i,j) = Cmn(:,k) where k = j + (1/2)*(i-1)*(i-2)
    figure; plot( f, 10*log10( Cmn(:,iC) ) ); title('Avg. Routine:: Coherence-Magnitude');
    %figure; plot( f, phimn(:,iC) ); title('Avg. Routine:: Coherence-Phase');
    disp( ['Confidence levelfor C (confC) at (1-p) level: ' num2str( ConfC(iC)) ] );
end

%
% Use to check against Chronux: only for equal length segments
%
if chronux
    cdata = repmat( fulldata( sMarkers(1,1):sMarkers(1,2), chIndices(i1) ), 1, nSegments ); % round(nSamples/2) x nSegments
    params.trialave = 1;
    if gram==0
        [ cS, cf, cSerr ] = mtspectrumc( cdata, params );
        figure; plot( cf, 10*log10( cS )); title('Chronux:: Spectrum');
        figure; plot( cf, 10*log10(cSerr(1,:)), cf, 10*log10(cSerr(2,:)) ); title('Chronux Error-Bar Computations');
        figure; plot( cf, 10*log10( cS ) - 10*log10( S(:,i1) )); title('Error in Spectrum = |New Routines - Chronux|');
        figure; plot( cf, 10*log10(cSerr(1,:)) - 10*log10(Serr(1,:,i1)), cf, 10*log10(cSerr(2,:)) - 10*log10(Serr(2,:,i1)) );title('Error in Error-Bar Computations = |New Routines - Chronux| ');
    elseif gram==1
        cdata2 = repmat( fulldata( sMarkers(1,1):sMarkers(1,2), chIndices(i2) ), 1, nSegments ); % round(nSamples/2) x nSegments
        params.trialave = 1;
        
        [cC,cphi,cS12,cS1,cS2,cf,cconfC,cphistd,cCerr]=coherencyc( cdata, cdata2, params );
        
        figure; plot( cf, 10*log10( cC(:,1) ) ); title('Chronux:: Coherence-Magnitude');
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
        figure; plot( cf, cphistd(1,:,1) -  PhiStd(1,:,iC) ); title('Error in PhiStd-1');
        figure; plot( cf, cphistd(2,:,1) -  PhiStd(2,:,iC) ); title('Error in PhiStd-2');
        figure; plot( cf, abs(cCerr(1,:,1) -  Cerr(1,:,iC)), cf, abs(cCerr(2,:,1) -  Cerr(2,:,iC))  ); title('Error in Abs(Coherence)-Error-Bar Computations = |New Routines - Chronux|');
    end
end

    



