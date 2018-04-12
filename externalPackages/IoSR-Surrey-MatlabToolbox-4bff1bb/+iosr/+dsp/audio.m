classdef (Abstract) audio < matlab.mixin.Copyable
%AUDIO Abstract superclass providing audio-related properties and methods.
% 
%   The IOSR.DSP.AUDIO superclass is an abstract class that defines
%   audio-related properties and methods. As an abstract class, it cannot
%   be instantiated.
% 
%   IOSR.DSP.AUDIO is a subclass of MATLAB.MIXIN.COPYABLE, which is a
%   subclass of the HANDLE class. Hence subclasses of AUDIO are handle
%   classes.
% 
%   IOSR.DSP.AUDIO properties:
%       filename        - Name of the corresponding audio file
%       fs              - Sampling frequency (Hz)
%       signal          - The sampled data (read-only)
% 
%   IOSR.DSP.AUDIO methods:
%       sound           - Replay the audio signal
%     Static methods:
%       up_down_mix     - Remix a source to a designated channel count
%       normalize       - Normalize an audio signal
%       spat            - Spatialise a sound using a SOFA file
% 
%   See also IOSR.BSS.MIXTURE, IOSR.BSS.SOURCE, HANDLE,
%            MATLAB.MIXIN.COPYABLE.

%   Copyright 2016 University of Surrey.
    
    properties (AbortSet)
        filename            % Name of the audio file
        fs = 16000          % Sampling frequency (Hz)
    end
    
    properties (Abstract, Dependent, SetAccess = protected)
        signal              % The sampled data (read-only)
    end
    
    properties (Access = protected)
        rendered = false    % Flag indicates whether files match object
    end
    
    methods % public methods
        
        function obj = audio(filename,fs)
        %AUDIO Create the audio object
            
            if nargin > 0
                obj.filename = filename;
            end
            if nargin > 1
                obj.fs = fs;
            end
            
        end
        
        function sound(obj)
        %SOUND Replay the audio signal
        %
        %   IOSR.DSP.AUDIO.SOUND() replays the audio signal.
            
            obj.replay(obj.signal)
            
        end
        
        % validate properties
        
        % validate filename
        function set.filename(obj,val)
            assert(ischar(val) || isempty(val), 'iosr:audio:invalidFile', 'FILENAME must be a char array or an empty array')
            obj.filename = val;
            obj.property_changed('filename',val);
        end
        
        % validate fs
        function set.fs(obj,val)
            assert(isscalar(val), 'iosr:audio:invalidFs', 'FS must be a scalar')
            assert(val > 0, 'iosr:audio:invalidFs', 'FS must be greater than 0')
            obj.fs = val;
            obj.property_changed('fs',val);
        end
        
    end
    
    methods (Access = protected)
        
        function property_changed(obj,name,val) %#ok<INUSD>
        %PROPERTY_CHANGED handle property changes
            
            obj.rendered = false;
        end
        
        function replay(obj,x)
        %REPLAY Replay any audio signal
            
            sound(obj.normalize(x),obj.fs) % replay
            
        end
       
        function set_properties(obj,vgin,propNames)
        %SET_PROPERTIES Parse varargin input to set object properties
        % 
        %   IOSR.DSP.AUDIO.SET_PROPERTIES(VGIN,PROPNAMES) searches through
        %   the cell array VGIN (a cell array of parameter/value pairs) for
        %   properties specified in the cell array PROPNAMES. Any matching
        %   properties are written to OBJ.

            % count arguments
            nArgs = length(vgin);
            if round(nArgs/2)~=nArgs/2
               error('iosr:audio:nameValuePair',[class(obj) ' needs propertyName/propertyValue pairs'])
            end
            % overwrite defults
            for pair = reshape(vgin,2,[]) % pair is {propName;propValue}
               IX = strcmpi(pair{1},propNames); % find match parameter names
               if any(IX)
                  % do the overwrite
                  obj.(propNames{IX}) = pair{2};
               else
                  error('iosr:audio:unknownOption','%s is not a recognized parameter name',pair{1})
               end
            end

        end
        
        function cpObj = copyElement(obj)
        %COPYELEMENT Overload copy method with additional functionality
            
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
        end
        
    end
    
    methods (Static)
        
        function [y,gain] = normalize(x)
        %NORMALIZE Normalize an audio signal
        %
        %   Y = IOSR.DSP.AUDIO.NORMALIZE(X) normalises X such that
        %   ABS(Y(:))<1.
        % 
        %   [Y,GAIN] = IOSR.DSP.AUDIO.NORMALIZE(X) returns the scalar gain
        %   GAIN applied to X in order to normalize it.
            
            gain = (1-((2*16)/(2^16)))/max(abs(x(:)));
            y = x.*gain;
            
        end
        
        function y = spat(SOFAfile,x,azimuth,elevation,fs)
        %SPAT Spatialise a sound using a SOFA file
        % 
        %   IOSR.DSP.AUDIO.SPAT(SOFAFILE,X,AZIMUTH,ELEVATION) spatialises
        %   the signal X using spatialisation data contained in the SOFA
        %   file SOFAFILE at the corresponding AZIMUTH and ELEVATION. The
        %   sampling frequency is that of the SOFA data. X should be a
        %   column vector; if it is not, it will be down-mixed to a column
        %   vector.
        %   
        %   IOSR.DSP.AUDIO.SPAT(...,FS) uses the sampling frequency FS for
        %   the spatialisation. If the SOFA data are not at the sampling
        %   rate FS, the data are resampled.
            
            assert(size(x,1)>1, 'iosr:audio:invalidX', 'X should be a column vector. Matrices will be downmixed to a column vector');
        
            if size(x,2)>1 % mix to mono
                x = audio.up_down_mix(x,1);
            end
            SOFAobj = SOFAload(SOFAfile); % load SOFA object
            if ~exist('fs','var')==1 % set default FS
                fs = SOFAobj.Data.SamplingRate;
            end
            
            % convert navigational coordinates to spherical coordinates
            if azimuth<0
                [azimuth,elevation] = nav2sph(azimuth,elevation);
            end
            
            % find index of corresponding IR is SOFA obj
            dist = (SOFAobj.SourcePosition(:,1)-azimuth).^2 + ...
                (SOFAobj.SourcePosition(:,2)-elevation).^2;
            [~,idx] = min(dist);
            
            % return IR
            IR = squeeze(SOFAobj.Data.IR(idx,:,:))';
            
            % resample?
            if SOFAobj.Data.SamplingRate~=fs
                disp(['Resampling audio to ' num2str(fs) ' Hz.'])
                IR = resample(IR,fs,SOFAobj.Data.SamplingRate);
            end
            
            % convolve
            y = zeros(length(x)+size(IR,1)-1,size(IR,2));
            for c = 1:size(IR,2)
                y(:,c) = iosr.dsp.convFft(x,IR(:,c));
            end
            
            % crop
            y = y(1:length(x),:);
            
        end
        
        function y = up_down_mix(x,N)
        %UP_DOWN_MIX Remix a source to a designated channel count
        % 
        %   Y = IOSR.DSP.AUDIO.UP_DOWN_MIX(X,N) mixes the signal X to have
        %   the specified channel count N. Data for each channel in X
        %   should be stored down its columns. For input channel count M,
        %   the up-/down-mix coefficients are calculated as
        %   ONES(M,N)./(M*N).
            
            M = size(x,2);
            y = x*(ones(M,N)./(M*N));
            
        end
        
    end
        
    methods (Static, Access = protected)
        
        function ensure_path(filename)
        %ENSURE_PATH Ensure path exists for specified filename
            
            % get path
            filepath = fileparts(filename);
            if ~isempty(filepath)
                % ensure path exits
                if exist(filepath,'dir')~=7
                    mkdir(filepath)
                end
            end
            
        end
        
    end
    
end
