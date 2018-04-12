classdef mixture < iosr.dsp.audio
%MIXTURE Class of sound source separation mixture.
% 
%   iosr.bss.mixture objects contain information about a mixture of sound
%   sources.
% 
%   IOSR.BSS.MIXTURE is a subclass of IOSR.DSP.AUDIO.
% 
%   IOSR.BSS.MIXTURE properties:
%       azi_sep         - The azimuthal separation of the widest sources
%                         (read-only)
%       decomp          - The result of the time-frequency decomposition
%                         (mixture) (read-only)
%       decomp_i        - The result of the time-frequency decomposition
%                         (interferer) (read-only)
%       decomp_t        - The result of the time-frequency decomposition
%                         (target) (read-only)
%       decomposition   - Set the time-frequency decomposition. The options
%                         are:
%                             'stft'        : short-time fourier transform
%                                             (default)
%                             'gammatone'   : gammatone filterbank
%       elevation       - The median elevation of the mixture (read-only)
%       filename_t      - Name of the target audio file (based on the
%                         mixture filename) (read-only)
%       filename_i      - Name of the interferer audio file (based on the
%                         Mixture filename) (read-only)
%       gammatone       - Settings for the gammatone filterbank
%                         decomposition. The property is a structure
%                         containing the following fields:
%                             'cfs'     : the centre frequencies of the
%                                         gammatone filterbank
%                             'frame'   : the frame length (in samples)
%       sofa_path       - Path to a SOFA file containing spatialisation data
%       ibm             - The ideal binary mask
%       irm             - The ideal ratio mask
%       int_fns         - A char array containing the filenames of all of
%                         the interfering sources (read-only)
%       interferers     - An array of interferer sources of type
%                         iosr.bss.source
%       signal_t        - The sampled data (target) (read-only)
%       signal_i        - The sampled data (interferer) (read-only)
%       stft            - Settings for the STFT decomposition. The property
%                         is a structure containing the following fields:
%                             'hop'     : the hop size of the STFT
%                             'win'     : the STFT window
%                         See iosr.dsp.stft for more information
%       target          - The target source of type iosr.bss.source
%       tfPower         - The (frame-based) time-frequency power (mixture)
%                         (read-only)
%       tfPower_i       - The (frame-based) time-frequency power
%                         (interferer) (read-only)
%       tfPower_t       - The (frame-based) time-frequency power (target)
%                         (read-only)
%       tir             - The target-to-interferer ratio (target or
%                         interferers are attenuated in order that their
%                         RMS amplitudes have this ratio)
%       wdo             - The w-disjoint orthogonality (read-only)
%       wdo_lw          - The loudness-weighted w-disjoint orthogonality
%                         (read-only)
% 
%   IOSR.BSS.MIXTURE methods:
%       mixture         - Create the mixture
%       clearCache      - Clears the internally cached data, reducing
%                         file size.
%       copy            - Create an independent copy of the mixture, its
%                         sources, and any rendered files
%       applyMask       - Apply a time-frequency mask
%       deleteFiles     - Delete the mixture audio files
%       sound_t         - Replay the target
%       sound_i         - Replay the interferer
%       write           - Save the mixture, target, and interferers to audio files
%     Static methods:
%       maskCentroids   - Calculate the centroids of a time-frequency mask
%       mixWdo          - WDO of a mixture
%       mixWdo_lw       - Loudness-weighted WDO of a mixture
%       mixWdo_Stokes   - WDO calculated using Stokes's method
% 
%   Note that target and interferer properties may be modified as
%   MIXTURE.TARGET.PROPERTY_NAME and MIXTURE.INTERFERERS(N).PROPERTY_NAME,
%   except for the sampling frequency FS, which cannot be overridden. This
%   ensures that the sampling frequencies are identical for the mixture and
%   its sources.
% 
%   See also IOSR.DSP.AUDIO, IOSR.BSS.SOURCE, SOFALOAD, IOSR.DSP.STFT,
%   IOSR.AUDITORY.GAMMATONEFAST.

%   Copyright 2016 University of Surrey.

%TODO: Add more decompositions.
    
    properties (AbortSet)
        decomposition = 'stft'  % The time-frequency decomposition
        gammatone = struct      % Settings for the gammatone filterbank decomposition
        sofa_path               % Path to a SOFA file containing spatialisation data
        stft = struct           % Settings for the STFT decomposition
        tir = 0                 % The target to interferer ratio
        interferers             % An array of interferer sources of type iosr.bss.source
        target                  % The target source of type iosr.bss.source
    end
    
    properties (Dependent, SetAccess = protected)
        signal          % The sampled data (read-only)
    end
        
    properties (Dependent, SetAccess = private)
        azi_sep     % The azimuthal separation of the widest sources (read-only)
        decomp      % The result of the time-frequency decomposition (mixture) (read-only)
        decomp_i    % The result of the time-frequency decomposition (interferer) (read-only)
        decomp_t    % The result of the time-frequency decomposition (target) (read-only)
        elevation   % The median elevation of the mixture (read-only)
        filename_t  % Name of the target audio file (read-only)
        filename_i  % Name of the interferer audio file (read-only)
        ibm         % Ideal binary mask (read-only)
        irm         % Ideal ratio mask (read-only)
        int_fns     % Filenames of all of the interfering sources (read-only)
        numchans    % Number of audio channels in the mixture (read-only)
        signal_t    % Return target (read-only)
        signal_i    % Return interferer (read-only)
        tfPower     % The time-frequency power (mixture) (read-only)
        tfPower_i   % The time-frequency power (interferer) (read-only)
        tfPower_t   % The time-frequency power (target) (read-only)
        wdo         % Return the w-disjoint orthogonality metric (read-only)
        wdo_lw      % Return the loudness-weighted w-disjoint orthogonality metric (read-only)
        wdo_stokes  % Return the w-disjoint orthogonality metric using Stokes's method (read-only)
    end
    
    properties (Access = private)
        cached_decomp               % the cached mixture decomposition
        decomp_cached = false       % whether the mixture decomposition is cached
        cached_decomp_i             % the cached interferer decomposition
        decomp_i_cached = false     % whether the interferer decomposition is cached
        cached_decomp_t             % the cached target decomposition
        decomp_t_cached = false     % whether the target decomposition is cached
        cached_tfPower              % the cached mixture TF power
        tfPower_cached = false      % whether the mixture TF power is cached
        cached_tfPower_i            % the cached interferer TF power
        tfPower_i_cached = false    % whether the interferer TF power is cached
        cached_tfPower_t            % the cached target TF power
        tfPower_t_cached = false    % whether the target TF power is cached
    end
    
    methods
        
        % constructor
        function obj = mixture(target,interferers,varargin)
        %MIXTURE Create a mixture
        % 
        %   OBJ = IOSR.BSS.MIXTURE(TARGET,INTERFERER) creates a mixture by
        %   summing together the target source and the interferer
        %   source(s). The sources are mixed together such that their RMS
        %   amplitudes are equal (target-to-interferer ratio is 0dB). The
        %   OBJ sampling rate is equal to the target sampling rate.
        %   Information about the sources' spatial location is ignored.
        %   
        %   OBJ = IOSR.BSS.MIXTURE(...,'PARAMETER',VALUE) allows additional
        %   options to be specified. The options are ({} indicate
        %   defaults):
        %   
        %       'decomposition' : {'stft'} | 'gammatone'
        %       'filename'      : {[]} | str
        %           A filename used when writing the file with the write()
        %           method. The filename may also be set when calling the
        %           write() method. Filenames for the target and interferer
        %           are determined automatically, by append the filename
        %           with '_target' and '_interferer' respectively.
        %       'fs'            : {obj.target.fs} | scalar
        %           The sampling frequency of the mixture. All HRTFs and/or
        %           sources will be resampled to this frequency each time
        %           the signal is requested.
        %       'gammatone'
        %           Settings for the gammatone filterbank
        %           decomposition. The property is a structure
        %           containing the following fields:
        %               'cfs'     : the centre frequencies of the gammatone
        %                           filterbank (default is 64 channels
        %                           equally spaced on the ERB scale between
        %                           20Hz and the Nyquist limit
        %               'frame'   : the frame length (in samples) (the
        %                           default is 20ms)
        %       'sofa_path'     : {[]} | str
        %           A path to a SOFA file containing HRTFs that are
        %           convolved with sources in order to generate the
        %           mixture.
        %       'stft'
        %           Settings for the STFT decomposition. The property
        %           is a structure containing the following fields:
        %               'hop'     : the hop size of the STFT (the default
        %                           is 512)
        %               'win'     : the STFT window (the default is 1024)
        %       'tir'           : {0} | scalar
        %           The RMS ratio of the target and interfer sources.
        %           Interferer sources are individually set to this level
        %           prior to their summation.
        %
        %   OBJ = IOSR.BSS.MIXTURE creates an empty mixture object with
        %   empty target and interferer sources.
        %         
        %   To speed up time-frequency operations, the  IOSR.BSS.MIXTURE
        %   class caches time-frequency data internally, which can lead to
        %   large file sizes if the object is saved. Use the CLEARCACHE()
        %   method to empty the internal cache.
        %
        %   Note that this is a handle class, as is IOSR.BSS.SOURCE. Target
        %   and interferer(s) are hence passed by reference. Use the COPY()
        %   method to create an independent copy of the mixture and its
        %   sources.
            
            if nargin > 0
                
                assert(nargin>1, 'iosr:mixture:nargin', 'Not enough input arguments')
        
                propNames = {'filename','fs','sofa_path','tir','decomposition','stft','gammatone'};

                % set sources
                obj.target = target;
                obj.interferers = interferers;

                % defaults
                obj.fs = obj.target.fs;
                obj.sofa_path = [];
                obj.tir = 0;
                obj.rendered = false;
                
                % default decomposition settings
                obj.stft = struct('win',1024,'hop',512);
                obj.gammatone = struct(...
                    'cfs',iosr.auditory.makeErbCFs(20, obj.fs/2, 64),...
                    'frame', round(0.01*obj.fs)...
                );

                % read parameter/value inputs
                if nargin > 1 % if parameters are specified
                    obj.set_properties(varargin,propNames)
                end

                % set sample rate to SOFA file if set and fs not specified
                if ~isempty(obj.sofa_path) && all(~strcmp('fs',varargin))
                    SOFAobj = SOFAload(obj.sofa_path); % load SOFA object
                    obj.fs = SOFAobj.Data.SamplingRate;
                end

                % ensure fs for sources matches instance
                obj.target.fs = obj.fs;
                for n = 1:numel(obj.interferers)
                    obj.interferers(n).fs = obj.fs;
                end
                
                % set parent
                obj.target.parent = obj;
                for n = 1:numel(obj.interferers)
                    obj.interferers(n).parent = obj;
                end
                
            end
            
        end
        
        function z = applyMask(obj,m)
        %APPLYMASK Apply a time-frequency mask to the mixture
        % 
        %   Z = IOSR.BSS.MIXTURE.APPLYMASK(M) applies the time-frequency
        %   mask M to the mixture. The mixture is transformed according to
        %   the objects decomposition settings. The transform is then
        %   multiplied with M. Lastly, the result is then inverse
        %   transformed (or reynthesised). The time-domain output Z is
        %   returned.
            
            switch lower(obj.decomposition)
                case 'stft'
                    s = obj.decomp;
                    % correct dims
                    d = obj.stft_shifted_dims(s);
                    s = ipermute(s,d);
                    m = ipermute(m,d);
                    for c = 1:obj.numchans
                        % apply mask
                        z(:,c) = iosr.bss.applyMask( ...
                            s(:,:,c), ...
                            m(:,:,min(c,size(m,3))), ...
                            obj.stft.win, ...
                            obj.stft.hop, ...
                            obj.fs ...
                        ); %#ok<AGROW>
                    end
                case 'gammatone'
                    x = obj.signal;
                    for c = 1:obj.numchans
                        z(:,c) =  iosr.bss.resynthesise( ...
                            x(:,c), ...
                            obj.fs, ...
                            iosr.bss.cfs2fcs(obj.gammatone.cfs, obj.fs), ...
                            m(:,:,c), ...
                            'frame_length', obj.gammatone.frame, ...
                            'filter', 'sinc', ...
                            'kernel', gausswin(obj.gammatone.frame) ...
                        ); %#ok<AGROW>
                    end
            end
            
            z = obj.setlength(z,length(obj.signal));
            
        end
        
        function sound_t(obj)
        %SOUND_T Replay the target signal
        %
        %   IOSR.BSS.MIXTURE.SOUND_T() replays the target signal.
            
            obj.replay(obj.signal_t)
            
        end
        
        function sound_i(obj)
        %SOUND_I Replay the interferer signal
        %
        %   IOSR.BSS.MIXTURE.SOUND_I() replays the interferer signal.
            
            obj.replay(obj.signal_i)
            
        end
        
        function write(obj,filename)
        %WRITE Save the mixture to an audio file
        %
        %   IOSR.BSS.MIXTURE.WRITE() writes the mixture to an audio file
        %   specified by MIXTURE.FILENAME, and also writes the target and
        %   interferer signals to automatically-determined file names
        %   (MIXTURE.FILENAME_T and MIXTURE.FILENAME_I respectively).
        %
        %   IOSR.BSS.MIXTURE.WRITE(FILENAME) uses the specified FILENAME
        %   and updates MIXTURE.FILENAME.
            
            if exist('filename','var')==1
                fn = filename;
            else
                fn = obj.filename;
            end
        
            if obj.rendered && exist(fn, 'file') == 2
                if strcmp(obj.filename,fn)
                    % nothing to do
                else
                    copyfile(obj.filename,fn);
                    copyfile(obj.filename_t,obj.make_target_filename(fn));
                    copyfile(obj.filename_i,obj.make_interferer_filename(fn));
                end
            else
                % check filename is valid
                obj.filename = fn;
                assert(ischar(obj.filename) && ~isempty(obj.filename),'iosr:mixture:invalidFilename','FILENAME must be a non-empty char array. Set filename as MIXTURE.FILENAME or MIXTURE.WRITE(FILENAME).')

                % normalize
                M = obj.signal;
                T = obj.signal_t;
                I = obj.signal_i;
                [~,gain] = obj.normalize([M T I]);
                M = M.*gain;
                T = T.*gain;
                I = I.*gain;

                % ensure path
                obj.ensure_path(obj.filename)

                % write audio files
                audiowrite(obj.filename,M,obj.fs);
                audiowrite(obj.filename_t,T,obj.fs);
                audiowrite(obj.filename_i,I,obj.fs);

                % set flag to indicate mixture has been rendered to an audio file
                obj.rendered = true;
                
            end
            
        end
        
        function deleteFiles(obj)
        %DELETEFILES Delete the mixture audio files
        
            if exist(obj.filename, 'file') == 2
                delete(obj.filename);
            end
            if exist(obj.filename_t, 'file') == 2
                delete(obj.filename_t);
            end
            if exist(obj.filename_i, 'file') == 2
                delete(obj.filename_i);
            end
            obj.rendered = false;
            
        end
        
        function clearCache(obj)
        %CLEARCACHE Clear the internal cache
            
            obj.updateCachedFalse();
            obj.cached_decomp = [];
            obj.cached_decomp_i = [];
            obj.cached_decomp_t = [];
            obj.cached_tfPower = [];
            obj.cached_tfPower_i = [];
            obj.cached_tfPower_t = [];
        end

        % set/validate properties
        
        % set decomposition
        function set.decomposition(obj,val)
            assert(any(strcmpi(val,{'stft','gammatone'})), 'iosr:mixture:invalidDecomposition', '''decomposition'' must be ''stft'' or ''gammatone''');
            obj.updateCachedFalse();
            obj.decomposition = val;
        end
        
        % set gammatone settings
        function set.gammatone(obj,val)
            assert(isstruct(val), 'iosr:mixture:invalidGammatoneStruct', '''gammatone'' must be a struct')
            assert(all(ismember(fieldnames(val)',{'cfs','frame'})), 'iosr:mixture:invalidGammatoneFields', 'gammatone structure must contain fields ''cfs'' and ''frame''')
            obj.updateCachedFalse;
            obj.gammatone = val;
        end
        
        % set stft settings
        function set.stft(obj,val)
            assert(isstruct(val), 'iosr:mixture:invalidStftStruct', '''stft'' must be a struct')
            assert(all(ismember(fieldnames(val)',{'win','hop'})), 'iosr:mixture:invalidStftFields', 'stft structure must contain fields ''win'' and ''hop''')
            obj.updateCachedFalse;
            obj.stft = val;
        end
        
        % set tir
        function set.tir(obj,val)
            obj.tir = val;
            obj.property_changed('tir',val);
        end
        
        % validate sofa_path
        function set.sofa_path(obj,val)
            assert(ischar(val) || isempty(val), 'iosr:mixture:invalidSofaPath', 'sofa_path must be a char array or an empty array')
            if ~isempty(val)
                assert(exist(val,'file')==2, 'iosr:mixture:invalidSofaPath', 'SOFA file does not exist')
            end
            obj.sofa_path = val;
            obj.property_changed('sofa_path',val);
        end
        
        % validate interferers
        function set.interferers(obj,val)
            assert(isa(val,'iosr.bss.source'), 'iosr:mixture:invalidInterferers', 'INTERFERERS must be of type source')
            obj.interferers = val;
            obj.property_changed('interferers',val);
        end
        
        % validate target
        function set.target(obj,val)
            assert(isa(val,'iosr.bss.source') && numel(val)==1, 'iosr:mixture:invalidTarget', 'TARGET must be a scalar of type source')
            obj.target = val;
            obj.property_changed('target',val);
        end
        
        % dependent properties
        
        % get decomposition result (mixture)
        function s = get.decomp(obj)
            if obj.decomp_cached
                s = obj.cached_decomp;
            else
                s = obj.decompose(obj.signal);
                obj.cached_decomp = s;
                obj.decomp_cached = true;
            end
        end
        
        % get decomposition result (interferer)
        function s = get.decomp_i(obj)
            if obj.decomp_i_cached
                s = obj.cached_decomp_i;
            else
                s = obj.decompose(obj.signal_i);
                obj.cached_decomp_i = s;
                obj.decomp_i_cached = true;
            end
        end
        
        % get decomposition result (target)
        function s = get.decomp_t(obj)
            if obj.decomp_t_cached
                s = obj.cached_decomp_t;
            else
                s = obj.decompose(obj.signal_t);
                obj.cached_decomp_t = s;
                obj.decomp_t_cached = true;
            end
        end
        
        % get azimuthal separation
        function s = get.azi_sep(obj)
            [s,~] = get_loc(obj);
        end
        
        % get median elevation
        function e = get.elevation(obj)
            [~,e] = get_loc(obj);
        end
        
        % target filename
        function fn = get.filename_t(obj)    
            fn = obj.make_target_filename(obj.filename);
        end
        
        % interferer filename
        function fn = get.filename_i(obj)
            fn = obj.make_interferer_filename(obj.filename);
        end
        
        % ibm
        function m = get.ibm(obj)
            switch lower(obj.decomposition)
                case 'stft'
                    for c = 1:obj.numchans
                        [~,m(:,:,c)] = iosr.bss.idealMasks( ...
                            obj.decomp_t(:,:,c), ...
                            obj.decomp_i(:,:,c) ...
                        ); %#ok<AGROW>
                    end
                case 'gammatone'
                    T = obj.tfPower_t;
                    I = obj.tfPower_i;
                    m = +(T>I);
            end
        end
        
        % irm
        function m = get.irm(obj)
            switch lower(obj.decomposition)
                case 'stft'
                    for c = 1:obj.numchans
                        m(:,:,c) = iosr.bss.idealMasks( ...
                            obj.decomp_t(:,:,c), ...
                            obj.decomp_i(:,:,c) ...
                        ); %#ok<AGROW>
                    end
                case 'gammatone'
                    T = obj.tfPower_t;
                    I = obj.tfPower_i;
                    m = T./(T+I);
            end
            m(isnan(m) | isinf(m)) = 0;
        end
        
        % filenames of all interferers
        function fns = get.int_fns(obj)
            for n = 1:length(obj.interferers)
                if n==1
                    fns = obj.interferers(n).filename;
                else
                    fns = [fns ', ' obj.interferers(n).filename]; %#ok<AGROW>
                end
            end
        end
        
        % number of audio channels
        function n = get.numchans(obj)
            if obj.hrtf_is_set()
                s = SOFAload(obj.sofa_path);
                n = size(s.Data.IR,2);
            else
                iN = [obj.interferers.numchans];
                n = max([obj.target.numchans; iN(:)]);
            end
        end
        
        % return target signal
        function signal_t = get.signal_t(obj)
            if obj.rendered && exist(obj.filename_t,'file')==2 % don't bother calculating
                signal_t = audioread(obj.filename_t);
            else % calculate
                signal_t = return_source(obj,obj.target);
                if obj.tir<0
                    % attenuate according to TIR
                    Trms = iosr.dsp.rms(signal_t(:));
                    Irms = iosr.dsp.rms(obj.signal_i(:));
                    signal_t = signal_t./(Trms/Irms); % match to interferer
                    signal_t = signal_t.*(10^(obj.tir/20)); % attenuate
                end
            end
        end
        
        % return interferer signal
        function signal_i = get.signal_i(obj)
            if obj.rendered && exist(obj.filename_i,'file')==2 % don't bother calculating
                signal_i = audioread(obj.filename_i);
            else % calculate
                signal_i = zeros(1, obj.numchans); % initialise
                maxlength = 0; % initialise
                for n = 1:numel(obj.interferers) % step through each interferer
                    % source signal (ensure 2-channel)
                    s = obj.return_source(obj.interferers(n));
                    % ensure signals are same length, or zero-pad
                    maxlength = max([length(s) maxlength]);
                    s = obj.setlength(s,maxlength);
                    signal_i = obj.setlength(signal_i,maxlength);
                    % add source to interferer
                    signal_i = signal_i + s;
                end
                if obj.tir>=0
                    % attenuate according to TIR
                    Trms = iosr.dsp.rms(obj.signal_t(:));
                    Irms = iosr.dsp.rms(signal_i(:));
                    signal_i = signal_i./(Irms/Trms); % match to interferer
                    signal_i = signal_i./(10^(obj.tir/20)); % attenuate
                end
            end
        end
        
        % return mixture signal
        function signal = get.signal(obj)
            if obj.rendered && exist(obj.filename,'file')==2 % don't bother calculating
                signal = audioread(obj.filename);
            else % calculate
                % return target and interferers
                T = obj.signal_t;
                I = obj.signal_i;
                % mix, ensuring equal length
                maxlength = max([length(T) length(I)]);
                signal = obj.setlength(T,maxlength) + obj.setlength(I,maxlength);
            end
        end
        
        % return time-frequency power (mixture)
        function val = get.tfPower(obj)
            if obj.tfPower_cached
                val = obj.cached_tfPower;
            else
                val = obj.tf_power(obj.decomp);
                obj.cached_tfPower = val;
                obj.tfPower_cached = true;
            end
        end
        
        % return time-frequency power (interferer)
        function val = get.tfPower_i(obj)
            if obj.tfPower_i_cached
                val = obj.cached_tfPower_i;
            else
                val = obj.tf_power(obj.decomp_i);
                obj.cached_tfPower_i = val;
                obj.tfPower_i_cached = true;
            end
        end
        
        % return time-frequency power (target)
        function val = get.tfPower_t(obj)
            if obj.tfPower_t_cached
                val = obj.cached_tfPower_t;
            else
                val = obj.tf_power(obj.decomp_t);
                obj.cached_tfPower_t = val;
                obj.tfPower_t_cached = true;
            end
        end
        
        % return wdo
        function w = get.wdo(obj)
            w = obj.mixWdo(abs(obj.decomp_t), abs(obj.decomp_i));
        end
        
        % return loudness-weighted wdo
        function w = get.wdo_lw(obj)
            [S,f] = obj.decompose(obj.signal_t);
            w = obj.mixWdo_lw(abs(S), abs(obj.decomp_i), f);
        end
        
        % return Stokes's wdo
        function w = get.wdo_stokes(obj)
            w = obj.mixWdo_Stokes(obj.irm);
        end
        
    end
    
    methods (Static, Access = public)
       
        function [C, Ct, Cf] = maskCentroids(m, fs, decomposition, nfft, hop)
        %MASKCENTROIDS Return various time-frequency mask centroids
        %
        %   C = IOSR.BSS.MIXTURE.MASKCENTROIDS(M,FS,NFFT,HOP) returns the
        %   time-frequency mask centroid C calculated from mask M, sample
        %   rate FS, FFT-length NFFT, and and hop size HOP. The mask
        %   centroid is a unitless measure of the centroid of a
        %   time-frequency mask. The mask centroid is a measure of the
        %   frequency content of a time-frequency mask.
        %   
        %   [C,Ct,Cf] = IOSR.BSS.MIXTURE.MASKCENTROIDS(...) returns the
        %   time-related spectral centroid (in Hz) Ct and the
        %   frequency-related spectral centroid (in seconds) Cf.
            
            assert(ischar(decomposition), 'iosr:mixture:invalidDecomposition', '''decomposition'' must be a char array')
        
            numchans = size(m,3);
            frame_frequency = fs/hop;
            switch lower(decomposition)
                case 'stft'
                    assert(isscalar(nfft), 'iosr:mixture:invalidNfft', '''nfft'' must be a scalar')
                    quefrency = fs/nfft;
                case 'gammatone'
                    assert(isvector(nfft), 'iosr:mixture:invalidNfft', '''cfs'' must be a scalar')
                    cfs = nfft;
                    quefrency = fs/min(diff(iosr.auditory.erbRate2hz(cfs)));
                otherwise
                    error('iosr:mixture:unknownDecomp','Unknown decomposition.')
            end
            
            f_dct = repmat((((0:size(m,1)-1)./size(m,1)).*frame_frequency)',1,size(m,2));
            c_dct = repmat(((0:size(m,2)-1)./size(m,2)).*quefrency,size(m,1),1);
            
            C = zeros(1, numchans);
            Ct = zeros(1, numchans);
            Cf = zeros(1, numchans);
            for c = 1:numchans
                X = abs(dct(m(:,:,c)));
                sumX = sum(sum(X));
                C(c) = sum(sum(X.*f_dct.*c_dct))./sumX;
                Ct(c) = sum(sum(X.*f_dct))./sumX;
                Cf(c) = sum(sum(X.*c_dct))./sumX;
            end
            
        end
        
        function w = mixWdo(st, si)
        %MIXWDO Calculate the mixture w-disjoint orthogonality
        %
        %   W = IOSR.BSS.MIXTURE.MIXWDO(St,Si) calculates the w-disjoint
        %   orthogonality of a target and interferer(s) from their
        %   time-frequency magnitudes St and Si. The w-disjoint
        %   orthogonality is estimated by the two-dimensional
        %   cross-correlation coefficient of the source magnitudes.
            
            numchans = size(st,3);
            w = zeros(1,numchans);
            for c = 1:numchans
                w(c) = corr2(st(:,:,c),si(:,:,c));
            end
            
        end
        
        function w = mixWdo_lw(st, si, f)
        %MIXWDO_LW Calculate the loudness-weighted w-disjoint orthogonality
        %
        %   W = IOSR.BSS.MIXTURE.MIXWDO_LW(St,Si,F) calculates the
        %   loudness-weighted w-disjoint orthogonality of a target and
        %   interferer(s) from their time-frequency magnitudes St and Si.
        %   The weighting is calculated from the ISO-226 65-phon equal
        %   loudness contour using the frequency vector F containing the
        %   frequency of each bin of St and Si. The magnitudes are
        %   multiplied with the loudness-weighting factor prior to
        %   calculating the WDO.
        %   
        %   See also IOSR.AUDITORY.LOUDWEIGHT.
        
            lw = iosr.auditory.loudWeight(f(:), 65)';
            lw = repmat(lw,size(st,1),1,size(st,3));
            w = iosr.bss.mixture.mixWdo(st.*lw,si.*lw);
            
        end
        
        function w = mixWdo_Stokes(irm)
        %MIXWDO_LW Calculate the w-disjoint orthogonality using Stokes's method
        %
        %   W = IOSR.BSS.MIXTURE.MIXWDO_STOKES(IRM) calculates the
        %   w-disjoint orthogonality using Stokes's method [1]. The method
        %   utilizes the ideal ratio mask (IRM), which is 0.5 when the
        %   target and interferer have equal power, and tends towards 0 or
        %   1 when either source dominates.
        %   
        %   References
        %
        %   [1] Stokes, Tobias W. (2015) Improving the perceptual quality
        %       of single-channel blind audio source separation. Doctoral
        %       thesis, University of Surrey.
        
            h = hist(irm(:),11)./numel(irm);
            hw = [0:0.2:1 0.8:-0.2:0];
            w = sum(h.*hw);
            
        end
        
    end
        
    methods (Static, Access = private)
        
        function y = setlength(x,signal_length)
        %SETLENGTH Crop or zero-pad signal to specified length
            
            d = size(x);
            if size(x,1)>signal_length % need to crop
                subsidx = [{1:signal_length} repmat({':'},1,ndims(x)-1)];
                y = x(subsidx{:});
            elseif size(x,1)<signal_length % need to zero-pad
                y = [x; zeros([signal_length-size(x,1),d(2:end)])];
            else % do nothing
                y = x;
            end
            
        end
        
        function fn = append_filename(filename,append)
        %APPEND_FILENAME Append strings to MIXTURE.FILENAME
            
            if isempty(filename) % do nothing if filename is empty
                fn = [];
            else % append
                [filepath,name,ext] = fileparts(filename); % break up filename
                newname = [name append ext]; % append
                if isempty(filepath) % only filename specified
                    fn = newname;
                else % path specified
                    fn = [filepath filesep newname];
                end
            end
            
        end
        
        function d = stft_shifted_dims(s)
        %STFT_SHIFTED_DIMS Dim order for STFT, swapping time and bins
        
            dim_order = 1:ndims(s);
            if ndims(s) > 2 %#ok<ISMAT>
                d = dim_order([2 1 dim_order(3:end)]);
            else
                d = [2 1];
            end
        end
        
    end
        
    methods (Access = private)
        
        function c = hrtf_is_set(obj)
        %HRTF_IS_SET Determine whether HRTFs are specified
            
            c = ~isempty(obj.sofa_path);
            
        end
        
        function s = return_source(obj,src)
        %RETURN_SOURCE Return signal from specified source
            
            % ensure correct channel count
            if ~src.precomposed && obj.hrtf_is_set()
                src.numchans = 1; % signal will be spatialised
            else
                src.numchans = obj.numchans; % signal will be mixed directly
            end
            x = src.signal;
            % return/convolve signal
            if ~src.precomposed && obj.hrtf_is_set() % convolve
                s = obj.spat(obj.sofa_path,x,src.azimuth,src.elevation,obj.fs);
            else % return directly
                s = x; 
            end
            
        end
        
        function [s,em] = get_loc(obj)
        %GET_LOC Return location information
            
            % get angles from the sources
            a = zeros(numel(obj.interferers)+1,1);
            e = a;
            for n = 1:numel(obj.interferers)
                a(n) = obj.interferers(n).azimuth;
                e(n) = obj.interferers(n).elevation;
            end
            a(end) = obj.target.azimuth;
            e(end) = obj.target.elevation;
            
            % get separation
            if any(a<0) % assume angles are -179:180
                s = abs(max(a)-min(a));
            else % assume angles are 0:359
                s = (360-max(a))+min(a);
            end
            s = mod(s,180);
            
            % get elevation
            em = median(e);
            
        end
        
        function fn = make_target_filename(obj,filename)
        %MAKE_TARGET_FILENAME Return the automated target filename
        
            fn = obj.append_filename(filename,'_target');
            
        end
        
        function fn = make_interferer_filename(obj,filename)
        %MAKE_INTERFERER_FILENAME Return the automated interferer filename
        
            fn = obj.append_filename(filename,'_interferer');
            
        end
        
        function [s,f,t] = decompose(obj,x)
        %DECOMPOSE Transform a signal into the time-frequency domain
            
            switch lower(obj.decomposition)
                case 'stft'
                    for c = 1:obj.numchans
                        [s(:,:,c),f,t] = iosr.dsp.stft(x(:,c), obj.stft.win, obj.stft.hop, obj.fs); %#ok<AGROW>
                    end
                    % put time down column
                    s = permute(s,obj.stft_shifted_dims(s));
                    s = obj.setlength(s, obj.get_frame_count(length(obj.signal)));
                case 'gammatone'
                    f = obj.gammatone.cfs;
                    t = (0:length(x)-1)./obj.fs;
                    for c = 1:obj.numchans
                        s(:,:,c) = iosr.auditory.gammatoneFast(x(:,c), obj.gammatone.cfs, obj.fs); %#ok<AGROW>
                    end
            end
            
        end
        
        function y = tf_power(obj,x)
        %TF_POWER Calculate TF power.
        
            switch lower(obj.decomposition)
                case 'gammatone'
                    frame_count = obj.get_frame_count(size(x,1));
                    y = zeros(frame_count,length(obj.gammatone.cfs),obj.numchans);
                    for m = 1:frame_count
                        samples = (m-1)*obj.gammatone.frame+1:m*obj.gammatone.frame;
                        for N = 1:size(x,3)
                            y(m,:,N) = sum(x(samples,:,N).^2);
                        end
                    end
                case 'stft'
                    y = abs(x).^2;
            end
            
        end
        
        function frame_count = get_frame_count(obj, signal_length)
        %GET_FRAME_COUNT Calculate number of T-F frames.

            assert(isscalar(signal_length), 'iosr:mixture:invalidSignalLength', 'SIGNAL_LENGTH must be a scalar.')

            switch lower(obj.decomposition)
                case 'gammatone'
                    frame_count = floor(signal_length/obj.gammatone.frame);
                case 'stft'
                    frame_count = fix((signal_length-(obj.nfft-obj.stft.hop))/obj.stft.hop);
            end

        end
        
        function n = nfft(obj)
        %NFFT return the fft length
        
            if isscalar(obj.stft.win)
                n = obj.stft.win;
            else
                n = length(obj.stft.win);
            end
        end
        
        function updateCachedFalse(obj)
            obj.decomp_cached = false;
            obj.decomp_i_cached = false;
            obj.decomp_t_cached = false;
            obj.tfPower_cached = false;
            obj.tfPower_i_cached = false;
            obj.tfPower_t_cached = false;
        end
        
    end
    
    methods(Access = protected)
        
        function property_changed(obj,name,val)
        %PROPERTY_CHANGED handle property changes
        
            obj.rendered = false;
            obj.updateCachedFalse();
            
            switch lower(name)
                case 'fs'
                    obj.target.fs = val;
                    for n = 1:numel(obj.interferers)
                        obj.interferers(n).fs = val;
                    end
            end
            
        end
    
        function cpObj = copyElement(obj)
        %COPYELEMENT Overload copy method with additional functionality
        
            % Make a shallow copy of all properties
            cpObj = copyElement@iosr.dsp.audio(obj);
            
            % Changed rendered file name
            cpObj.filename = obj.append_filename(cpObj.filename,'_copy');
            
            % copy files
            if obj.rendered
                if exist(obj.filename,'file')==2
                    copyfile(obj.filename,cpObj.filename)
                end
                if exist(obj.filename_i,'file')==2
                    copyfile(obj.filename_i,cpObj.filename_i)
                end
                if exist(obj.filename_t,'file')==2
                    copyfile(obj.filename_t,cpObj.filename_t)
                end
            end
        
            % Make a deep copy of target and interferers
            cpObj.target = copy(obj.target);
            cpObj.interferers = copy(obj.interferers);
            
        end
        
    end
    
end
