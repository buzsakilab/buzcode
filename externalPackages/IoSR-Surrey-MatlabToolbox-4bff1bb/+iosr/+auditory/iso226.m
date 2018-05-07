function [spl, f, params] = iso226(phon,fq,sq)
%ISO226 ISO 226:2003 Normal equal-loudness-level contours
%   
%   [SPL,F] = IOSR.AUDITORY.ISO226(PHON) returns the sound pressure level
%   (SPL) (dB) of pure tone frequencies F (Hz) at the loudness level(s)
%   PHON. The values are calculated according to ISO 226:2003 using the
%   reference frequencies specified in the standard. According to the
%   standard, PHON is only valid at all frequencies if 20<=PHON<80
%   (although the function will return SPL values outside of this range).
%   If PHON is 0, the threshold of hearing is returned.
% 
%   PHON may be an array of any size; SPL and F will be of size
%   [1,29,M,N,P,...] where M,N,P,... are the dimensions of PHON.
%
%   [SPL,F] = IOSR.AUDITORY.ISO226(PHON,FQ) returns the SPL of the pure
%   tone frequencies in FQ at the specified loudness level(s). For
%   non-standard frequencies, the SPL is calculated by interpolating the
%   parameters used in its calculation. According to the standard, FQ is
%   only valid between 20 Hz and 12.5 kHz; the function will extrapolate
%   SPL values above 12.5 kHz by mirroring 20 Hz values at 20 kHz.
% 
%   FQ may be an array of any size; SPL and F will be of size
%   [Q,R,S,...,M,N,P,...] where Q,R,S,... are the dimensions of FQ.
%   
%   ... = IOSR.AUDITORY.ISO226(PHON,FQ,SQ) specifies whether singleton
%   dimensions will be removed from the output. With sq=false, singleton
%   dimensions will be retained (default), else they will be removed.
% 
%   ... = IOSR.AUDITORY.ISO226(PHON,[],SQ) uses the standard reference
%   frequencies for SPL calculations.
% 
%   [SPL,F,PARAMS] = IOSR.AUDITORY.ISO226(...) returns the reference
%   parameters used to calculate the normal equal-loudness-level contours.
%   PARAMAS is a structure with the following fields:
%       'f'         : the reference frequencies,
%       'alpha_f'   : the exponent of loudness perception,
%       'L_U'       : magnitude of the linear transfer function normalized
%                     at 1000 Hz, and
%       'T_f'       : the threshold of hearing.
% 
%   Example
% 
%       % Plot equal-loudness contours between 20 and 80 phon
% 
%       % Calculate SPLs
%       phons = 20:10:80;
%       [spl,f] = iosr.auditory.iso226(phons,[],true);
% 
%       % plot
%       figure; semilogx(f,spl)
%       set(gca,'xlim',[min(f(:)) max(f(:))])
%       legend(num2str(phons'),'location','southwest');
%       title('Equal loudness contours for different loudness levels (in phons)')
%       xlabel('Frequency [Hz]')
%       ylabel('SPL [dB]')
% 
%   See also IOSR.AUDITORY.LOUDWEIGHT.

%   Copyright 2016 University of Surrey.

    %% Check input

    if any(phon > 80)
        warning('iosr:iso226:phonRange','SPL values may not be accurate for loudness levels above 80 phon.')
    elseif any(phon(phon~=0) < 20)
        warning('iosr:iso226:phonRange','SPL values may not be accurate for loudness levels below 20 phon.')
    end

    if nargin>1
        if ~isempty(fq)
            if any(fq(:) < 20 | fq(:) > 4000) && any(phon > 90)
                warning('iosr:iso226:frequencyRange','ISO 226:2003 is valid for 20?4000 Hz only up to 90 phon. SPL values may be inaccurate.')
            elseif any(fq(:) < 5000 | fq(:) > 12500) && any(phon > 80)
                warning('iosr:iso226:frequencyRange','ISO 226:2003 is valid for 5000?12500 Hz only up to 80 phon. SPL values may be inaccurate.')
            elseif any(fq(:)>12500)
                warning('iosr:iso226:frequencyRange','ISO 226:2003 defines loudness levels up to 12.5 kHz. SPL values for frequencies above 12.5 kHz may be inaccurate.')
            end
            assert(all(fq(:)>=0), 'iosr:iso226:invalidFrequencies', 'Frequencies must be greater than or equal to 0 Hz.')
        end
    else
        fq = [];
    end

    if nargin<3
        sq = false;
    else
        assert(islogical(sq), 'iosr:iso226:invalidSq', 'sq must be logical.')
    end

    %% References

    % reference frequencies
    params.f = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 ...
        500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 ...
        6300 8000 10000 12500];

    % exponent of loudness perception
    params.alpha_f = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 ...
        0.367 0.349 0.330 0.315 0.301 0.288 0.276 0.267 0.259...
        0.253 0.250 0.246 0.244 0.243 0.243 0.243 0.242 0.242...
        0.245 0.254 0.271 0.301];

    % magnitude of linear transfer function normalized at 1 kHz
    params.L_U = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 ...
        -6.2 -4.5 -3.1 -2.0 -1.1 -0.4 0.0 0.3 0.5 0.0 -2.7 ...
        -4.1 -1.0 1.7 2.5 1.2 -2.1 -7.1 -11.2 -10.7 -3.1];

    % threshold of hearing
    params.T_f = [78.5 68.7 59.5 51.1 44.0 37.5 31.5 26.5 22.1 17.9...
        14.4 11.4 8.6 6.2 4.4 3.0 2.2 2.4 3.5 1.7 -1.3 -4.2...
        -6.0 -5.4 -1.5 6.0 12.6 13.9 12.3];

    %% Calculate

    % determine frequency range
    if isempty(fq)
        f = params.f;
    else
        f = fq;
    end

    % output size
    out_dims = [size(f) size(phon)];

    % independent outputs
    f_squeeze = zeros(numel(f),numel(phon));
    spl_squeeze = zeros(numel(f),numel(phon));

    % iterate through phons
    for p = 1:numel(phon)
        % frequencies for phon level
        f_squeeze(:,p) = f(:);
        % interpolate reference parameters
        if nargin>1
            if any(f_squeeze(:,p) > 12500)
                % extrapolate - mirror 20Hz behaviour at 20kHz
                f_r_extrap = [params.f 20000];
                alpha_f_r_extrap = [params.alpha_f params.alpha_f(1)];
                L_U_r_extrap = [params.L_U params.L_U(1)];
                T_f_r_extrap = [params.T_f params.T_f(1)];
            else
                f_r_extrap = params.f;
                alpha_f_r_extrap = params.alpha_f;
                L_U_r_extrap = params.L_U;
                T_f_r_extrap = params.T_f;
            end
            % interpolate parameters
            alpha_f = interp1(f_r_extrap, alpha_f_r_extrap, f_squeeze(:,p)', 'spline', 'extrap');
            L_U = interp1(f_r_extrap, L_U_r_extrap, f_squeeze(:,p)', 'spline', 'extrap');
            T_f = interp1(f_r_extrap, T_f_r_extrap, f_squeeze(:,p)', 'spline', 'extrap');
        else
            alpha_f = params.alpha_f;
            L_U = params.L_U;
            T_f = params.T_f;
        end
        % calculate SPL
        A_f = 0.00447 * ((10^(0.025*phon(p)))-1.15) + ...
            ((0.4*(10.^(((T_f+L_U)./10)-9))).^alpha_f);
        if phon(p) > 0
            spl_squeeze(:,p) = ((10./alpha_f).*log10(A_f)) - L_U + 94;
        else
            spl_squeeze(:,p) = T_f;
        end
    end

    % reshape outputs
    f = reshape(f_squeeze,out_dims);
    spl = reshape(spl_squeeze,out_dims);

    % remove singleton dimensions if requested
    if sq
        f = squeeze(f);
        spl = squeeze(spl);
    end

end
