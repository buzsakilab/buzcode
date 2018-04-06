function g = loudWeight(f,phon)
%LOUDWEIGHT Calculate loudness weighting coefficients
% 
%   G = IOSR.AUDITORY.LOUDWEIGHT(F) returns loudness-weighting linear-
%   magnitude coefficients for frequencies F (Hz). The function is based on
%   the loudness at 65 phons defined in ISO 226:2003. The coefficients are
%   scaled such that the G = 1 when F = 1000.
% 
%   G = IOSR.AUDITORY.LOUDWEIGHT(F,PHON) returns loudness weighting
%   coefficients at the loudness level PHON. PHON should be a scalar and,
%   according to the standard, is only valid at all frequencies such that
%   20<=PHON<80 (although the function will return extrapolated
%   coefficients outside of this range).
% 
%   The input f may be an array of any size. The outputs will be the same
%   size as f, with coefficients calculated for each element.
% 
%   G = IOSR.AUDITORY.LOUDWEIGHT(F,METHOD) returns loudness weighting
%   coefficients for a variety of methods. Specifying METHOD as 'A', 'C',
%   or 'Z' selects frequency weighting curves defined in IEC 61672-1:2013;
%   'ISO-226' selects a loudness weighting curve derived from ISO 226:2003
%   at 65 phons; 'B' selects a frequency weighting curve defined in IEC
%   60651:1979; 'D' selects a weighting curve defined in IEC 537:1976.
% 
%   See also IOSR.AUDITORY.ISO226, IOSR.AUDITORY.DUPWEIGHT.

%   Copyright 2016 University of Surrey.

    %% Check input

    assert(all(f(:)>=0), 'iosr:loudWeight:invalidF', 'f should be greater than or equal to zero!')

    if nargin<2
        phon = 65;
    end
    
    if ischar(phon)
        method = phon;
        if strcmpi(method,'iso-226')
            phon = 65;
        end
    elseif isnumeric(phon)
        method = 'iso-226';
        assert(isscalar(phon), 'iosr:loudWeight:invalidPhon', 'phon must be a scalar.')
    else
        error('iosr:loudWeight:invalidArg','Second argument should be a scalar or char array.')
    end
    
    %% coefficients
    
    fr = 1000;
    fL = 10^1.5;
    fH = 10^3.9;
    fA = 10^2.45;
    D = sqrt(0.5);
    b = (1/(1-D)) * ((fr^2) + (((fL^2)*(fH^2))/(fr^2)) - (D*((fL^2)+(fH^2))));
    c = (fL^2)*(fH^2);
    f1 = sqrt((-b - sqrt((b^2) - (4*c))) / (2));
    f4 = sqrt((-b + sqrt((b^2) - (4*c))) / (2));
    f2 = ((3 - sqrt(5)) / 2) * fA;
    f3 = ((3 + sqrt(5)) / 2) * fA;

    %% Calculate weighting coefficients
    
    switch lower(method)
        case 'a'
            g = a_weighting(f);
        case 'b'
            g = b_weighting(f);
        case 'c'
            g = c_weighting(f);
        case 'd'
            g = d_weighting(f);
        case 'z'
            g = z_weighting(f);
        case 'iso-226'
            % calculate weighting coefficients
            gdB = iosr.auditory.iso226(phon, 1000) - ...
                squeeze(iosr.auditory.iso226(phon, f));
            g = 10.^(gdB./20);
        otherwise
            error('iosr:loudWeight:unknownMethod','Unknown method.')
    end
    
    function w = a_weighting(f)
    %A_WEIGHTING return A-weighting magnitude coefficients
        function w = calculate(f)
            w = ((f4^2).*(f.^4))./...
                ( ((f.^2)+(f1^2)) .* sqrt((f.^2)+(f2^2)) .* sqrt((f.^2)+(f3^2)) .* ((f.^2)+(f4^2)) );
        end
        w = calculate(f)./calculate(1000);
    end

    function w = b_weighting(f)
    %B_WEIGHTING return B-weighting magnitude coefficients
        w = ((12200^2).*(f.^3))./...
            (((f.^2)+(20.6^2)).*sqrt((f.^2)+(158.5^2)).*((f.^2)+(12200^2)));
    end

    function w = c_weighting(f)
    %C_WEIGHTING return C-weighting magnitude coefficients
        function w = calculate(f)
            w = ((f4^2).*(f.^2)) ./...
                ( ((f.^2)+(f1^2)) .* ((f.^2)+(f4^2)) );
        end
        w = calculate(f)./calculate(1000);
    end

    function w = d_weighting(f)
    %D_WEIGHTING return D-weighting magnitude coefficients
        hf = (((1037918.48-(f.^2)).^2)+(1080768.16.*(f.^2)))./...
            (((9837328-(f.^2)).^2)+(11723776.*(f.^2)));
        w = (f./(6.8966888496476*(10^(-5)))).*sqrt(hf./(((f.^2)+79919.29).*((f.^2)+1345600)));
    end

    function w = z_weighting(f)
    %Z_WEIGHTING return Z-weighting magnitude coefficients
        w = ones(size(f));
    end

end
