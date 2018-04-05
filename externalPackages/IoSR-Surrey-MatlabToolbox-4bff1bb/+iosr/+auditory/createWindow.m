function OutWin = createWindow(varargin)
% Create a Hann or exp. window with specified onsets/offsets
%
% OutWin = iosr.auditory.createWindow(WinLen, Type, Duration)
%
% Other (non-symmetric) uses include:
%
% OutWin = iosr.auditory.createWindow(WinLen, ...
%   Hann, ...
%   OnDuration, ...
%   Hann, ...
%   OffDuration,)
% OutWin = iosr.auditory.createWindow(WinLen, ...
%   Hann, ...
%   OnDuration, ...
%   Exp, ...
%   OffDuration, ...
%   OffSlope)
% OutWin = iosr.auditory.createWindow(WinLen, ...
%   Exp, ...
%   OnDuration, ...
%   OnSlope, ...   
%   Hann, ...
%   OffDuration)
% OutWin = iosr.auditory.createWindow(WinLen, ...
%   Exp, ...
%   OnDuration, ...
%   OnSlope, ...   
%   Exp, ...
%   OffDuration, ...
%   OffSlope)
%
% OutWin = iosr.auditory.createWindow(type, samples)
% OutWin = iosr.auditory.createWindow(Ontype, Onsamples, ...
%     Offtype, Offsamples)
%
% In the first mode of operation, the onset and offset are assumed to be
% identical (symmetric window). In other modes, the OnRamp and OffRamp
% durations (and slopes if necessary) are specified independently
%
% OutWin: the created window
%
% WinLen: the total duration (in samples) of the desired window
%
% Type (inc. Ontype and Offtype):
%   
%   - 'Hann': a raised cosine ramp, following value should always be the
%   Duration of the ramp
%   
%   - 'Exp': an exponential ramp, following value should be always be the
%   Duration of the ramp, and then the Slope. e.g Slope > 3 very shallow,
%   slope < 0.1 very steep
%
% Duration: number of samples over which the OnRamp and OffRamp are applied
%
% Slope: the steepness of the exponential ramp
% 
% Example case: 100 samples onset Hann, followed by 50 unramped samples
% (i.e. value = 1), followed by 500 samples of shallow exponential offset
% ramp
%
%   OutWin = iosr.auditory.createWindow(650, 'Hann', 100, 'Exp', 500, 3)

%   Copyright 2016 University of Surrey.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Create Window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% test input types
assert(isnumeric(varargin{1}),['1st input should be window length' ...
'specified in samples']);
assert(size(varargin, 2)>2,'not enough input arguments');
assert(size(varargin, 2)<8,'too many input arguments');

%% extract inputs
nInputs = (size(varargin, 2));
WinLen = varargin{1};
Slope = [0 0];
switch nInputs

    case 3 % Hann or Rec symmetric window
        assert(~strcmp(varargin{2},('Exp')),'Must specify slope of exponential');
        assert(strcmp(varargin{2},('Hann')) || strcmp(varargin{2},('Rec')),'Input type not specified');
        RampOn = varargin{2};
        RampOff = varargin{2};
        RampOnDur = varargin{3};
        RampOffDur = varargin{3};
    case 4 % Exp symmetric window
        assert(strcmp(varargin{2},('Exp')),'wrong number of input arguments or ramp type');
        assert(varargin{4}~=0,'exponential slope may not be 0');
        RampOn = 'Exp';
        RampOff = 'Exp';
        RampOnDur = varargin{3};
        RampOffDur = varargin{3};
        Slope(1) = varargin{4};
    case 5 % Hann or Rec ramp on and ramp off
        assert(strcmp(varargin{2},('Hann')) || strcmp(varargin{2},('Rec')),'Input type not specified');
        assert(strcmp(varargin{4},('Hann')) || strcmp(varargin{4},('Rec')),'Input type not specified');
        RampOn = varargin{2};
        RampOff = varargin{4};
        RampOnDur = varargin{3};
        RampOffDur = varargin{5};
    case 6 % (Hann or Rec) + Exp, or Exp + (Hann or Rec)  
        if (strcmp(varargin{2},'Hann')||strcmp(varargin{2},'Rec'))
            RampOn = varargin{2};
            RampOnDur = varargin{3};
            RampOff = 'Exp';
            RampOffDur = varargin{5};
            Slope(2) = varargin{6};
        elseif strcmp(varargin{2},'Exp')
            RampOn = 'Exp';
            RampOnDur = varargin{3};
            assert(isnumeric(varargin{4}),'exponential slope should be specified as arg 4');
            assert(varargin{4}~=0,'exponential slope may not be 0');
            Slope(1) = varargin{4};
            RampOff = varargin{5};
            RampOffDur = varargin{6};
        else
            error('incorrect input arguments');
        end
case 7 % Ramp on = Exp, Ramp off = Exp
        assert(strcmp(varargin{2},('Exp')) && strcmp(varargin{2},('Exp')),'both ramps should be Exp for 7 input arguments');
        RampOn = 'Exp';
        RampOff = 'Exp';
        RampOnDur = varargin{3};
        Slope(1) = varargin{4};
        RampOffDur = varargin{6};
        Slope(2) = varargin{7};
    otherwise
        error('incorrect number of inputs');   
end
assert(WinLen>=((RampOnDur)+(RampOffDur)),'ramp lengths greater than signal duration!');

%% Test Hann gives correct answer for dummy example
assert( 0.01 > ( max(abs( AppWin(100,'Hann',10,'Hann',10) - ...
    [ 0; 0.0245; 0.0955; 0.2061; 0.3455; 0.5000; ...
    0.6545; 0.7939; 0.9045; 0.9755; ...
    ones(80,1); ...
    0.9755; 0.9045; 0.7939; 0.6545; 0.5000; ...
    0.3455; 0.2061; 0.0955; 0.0245; 0 ]))), ...
    'Hann not functioning corectly');

%% Do Work
OutWin = AppWin(WinLen,RampOn,RampOnDur,RampOff,RampOffDur,Slope);

end

% ----------------------------------------------------------
% DO AppWin
% ----------------------------------------------------------
function OutWin = AppWin(WinLen,RampOn,RampOnDur,RampOff,RampOffDur,Slope)

% generate rectangular window
OutWin = ones(WinLen,1);

%% Create the ramp on
switch RampOn
    case 'Hann'
        OutWin(1:RampOnDur) = (1-cos(pi/RampOnDur*[0:RampOnDur-1])) / 2;
        
    case 'Exp'
        OutWin(1:RampOnDur) = exp(([1:RampOnDur]./(Slope(1)*RampOnDur))-(1/Slope(1)));
        
    case 'Rec'
        
    otherwise
        error('ramp on unspecified');
end

%% Create the ramp off
switch RampOff
    case 'Hann'
        OutWin(end-RampOffDur+1:end) = (1+cos(pi/RampOffDur*[1:RampOffDur]))/2;
        
    case 'Exp'
        OutWin(end-RampOffDur+1:end) = exp((fliplr([1:RampOffDur])./(Slope(2)*RampOffDur))-(1/Slope(2)));
    
    case 'Rec'
            
    otherwise
        error('ramp off unspecified');
end

end
