%extracts spikes from raw signal
%applies band-pass filtering to raw signal before doing so.
%
%params.nrNoiseTraces: 0 if no noise should be estimated
%               >0 : # of noise traces to be used to estimate autocorr of
%               noise, returned in variable autocorr
%
%
%params.detectionMethod: 1 -> from power signal, 2 threshold positive, 3 threshold negative, 4 threshold abs, 5 wavelet
%params.detectionParams: depends on detectionMethod. 
%       if detectionmethod==1, detectionParams.kernelSize
%       if detectionmethod==4, detectionParams.scaleRanges (the range of scales (2 values))
%                              detectionParams.waveletName (which wavelet to use)
%
%params.peakAlignMethod: 1-> find peak, 2->none, 3->peak of power signal, 4->peak of MTEO signal.
%
%
%orig: urut/2004
%revised: urut/feb07
%
function [rawMean, filteredSignal, rawTraceSpikes,spikeWaveforms, spikeTimestamps, runStd2, upperlim, noiseTraces] = extractSpikes(rawSignal, Hd, params )
defineSortingConstants;

%set default parameters if they aren't set
detectionMethod=METHOD_EXTRACTION_POWER;  %default
if isfield(params,'detectionMethod')
    detectionMethod=params.detectionMethod;
else
    params.detectionMethod=detectionMethod;
end
if ~isfield(params,'peakAlignMethod')
   params.peakAlignMethod=1;  %default is normal findPeak method if is not set 
end

%calculate running mean raw signal
rawMean = runningAverage( rawSignal, 50);
filteredSignal = filterSignal( Hd,  rawSignal);


%
%calculate the to-be-thresholded signal, depending on the method used
%
searchInds=[];
runStd2=[];
switch(detectionMethod)
    case METHOD_EXTRACTION_POWER   %power method
        %calculate local energy
        kernelSize=25;

        if isfield(params,'detectionParams')
            kernelSize=params.detectionParams.kernelSize;
        end
        
        runStd2 = runningStd(filteredSignal, kernelSize);  %5=0.2ms. use ~ 1ms , because thats a typical spike width
        d0 = size(rawSignal,1) - size(runStd2,1);
        end0 = runStd2(end);
        runStd2(end:end+d0)=end0;

        %---STD
        upperlimFixed = mean( runStd2 ) + params.extractionThreshold * std(runStd2);    %extractionThreshold default is 5
        upperlim=ones(length(runStd2),1)*upperlimFixed;
    case METHOD_EXTRACTION_AMPP   %amplitude thresholding method (positive)
        runStd2=filteredSignal;
        upperlimFixed = params.extractionThreshold * median(abs(filteredSignal)/0.6745);
        upperlim=ones(length(runStd2),1)*upperlimFixed;

    case METHOD_EXTRACTION_AMPN   %amplitude thresholding method (negative)
        runStd2=-1*filteredSignal;
        upperlimFixed = params.extractionThreshold * median(abs(filteredSignal)/0.6745);
        upperlim=ones(length(runStd2),1)*upperlimFixed;
    case METHOD_EXTRACTION_AMPA   %amplitude thresholding negative+positive
        runStd2=abs(filteredSignal);
        upperlimFixed = params.extractionThreshold * median(runStd2/0.6745);
        upperlim=ones(length(runStd2),1)*upperlimFixed;
        
    case METHOD_EXTRACTION_WDM  %WDM method.
        runStd2=[];
        upperlim=[];
        
        L=params.extractionThreshold;
        
        plotFlag=0;
        commentFlag=1;
        
        %defaults
        scalesRange = [0.5 1.0]; %in ms
    	waveletName='haar'; 
        
        if isfield(params,'detectionParams')
            scalesRange=params.detectionParams.scalesRange;
            waveletName=params.detectionParams.waveletName;
        end
        
        timepointsWE = detect_spikes_wavelet( filteredSignal, params.samplingFreq/1000, scalesRange, 'reset', L, waveletName, plotFlag, commentFlag);
        searchInds=timepointsWE; %directly supply timepoints to extraction method; detection method wont threshold again in this case.        
    otherwise
        error('unknown detection method ');
end


%compute the signal for peak finding, if required
peakFindSignal=[];
switch( params.peakAlignMethod )
    case METHOD_PEAKFIND_POWER
        kernelSize=25;
        peakFindSignal = runningStd(filteredSignal, kernelSize);
        peakFindSignal = [powerSignal; zeros(kernelSize-1,1) ];          
    case METHOD_PEAKFIND_MTEO
        MTEOScales=[1 2 3]; %how many scales to use to compute the MTEO signal.
        peakFindSignal = MTEO( filteredSignal', MTEOScales);        
end

if params.prewhiten
          disp(['Prewhiten the signal - nr datapoints ' num2str(length(filteredSignal))]);	
          x= filteredSignal(1:400000);
          [a,e]=lpc(x,5);
          a=real(a);
          e=sqrt(e);
          filteredSignal2 = filter(a, e, filteredSignal );
          %filteredSignal=filteredSignal2;
end

if ~isfield(params,'nrNoiseTraces')
    params.nrNoiseTraces=0;
end

%tic
[rawTraceSpikes,spikeWaveforms, spikeTimestamps, noiseTraces] = detectSpikes(filteredSignal, rawMean, runStd2, upperlim, params, searchInds, peakFindSignal );
%toc
