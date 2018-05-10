function snr = calcSnr(output,target)
%CALCSNR Calculate the separation SNR
% 
%   SNR = IOSR.BSS.CALCSNR(OUTPUT,TARGET) calculate the separation
%   signal-to-noise ratio in dB for a re-synthesised output compared to the
%   ideal target. OUTPUT and TARGET should be vectors of equal length.

%   Copyright 2016 University of Surrey.

    % check input
    if ~isvector(output)
        error('iosr:calcSnr:invalidOutput','''output'' must be a vector')
    end
    if ~isvector(target)
        error('iosr:calcSnr:invalidTarget','''target'' must be a vector')
    end
    if numel(output)~=numel(target)
        error('iosr:calcSnr:inputLengths','inputs must be the same length')
    end
    if size(output,1)==1
        output = output';
    end
    if size(target,1)==1
        target = target';
    end

    % remove delay caused by convolution
    cc = xcorr(output.^2,target.^2);

    delay = find(cc==max(cc))-length(output);

    if delay > 0
            target = [zeros(delay,1); target];
            output = [output; zeros(delay,1)];
    elseif delay < 0
            delay = -delay;
            output = [zeros(delay,1); output];
            target = [target; zeros(delay,1)];
    end

    % account for arbitrary gain
    if sum(abs(output(:))) > 0
        G = output\target;
        output = output.*G;
    end
    
    snr = 10*log10(sum(target.^2)/sum((output-target).^2));

end
