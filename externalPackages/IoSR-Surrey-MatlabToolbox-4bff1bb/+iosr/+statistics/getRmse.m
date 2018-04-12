function RMSE = getRmse(X, y, varargin)
%GETRMSE Calculate the root-mean-square error between input data
%   
%   RMSE = IOSR.STATISTICS.GETRMSE(X, Y) calculates the RMSE of the inputs
%   data where
%
%       RMSE = sqrt( 1/N-D * SUM(err^2) )
%       err = abs(X-Y)
%
%   The input data X and Y can be vectors or matrices; D=1. X and Y must
%   have equal size.
%
%   RMSE = IOSR.STATISTICS.GETRMSE(X, Y, D) allows the degrees of freedom D
%   to be specified. The default is D=1 (you may wish to set it to 0, or
%   the degrees of freedom).
%
%   RMSEEPS = IOSR.STATISTICS.GETRMSE(X, Y, D, EPSILON) calculates
%   epsilon-insensitive RMSE (RMSE*). The parameter EPSILON is a threshold
%   for err values: err values less than the respective EPSILON value are
%   set to 0; err values greater than EPSILON have EPSILON subtracted from
%   them. This allows for the calculation of error 'after' a certain effect
%   (such as subjective error). RMSEEPS will ALWAYS be lower than RMSE. The
%   value of EPSILON will generally be set to half the 95% confidence
%   interval or similar

%   Copyright 2016 University of Surrey.

    % test input types
    assert(isnumeric(X)&isnumeric(y),'Input data must be numeric!');
    assert(length(size(X))==length(size(y)),'Input data X and y must be the same size');
    assert(all(size(X)==size(y)),'Input data X and y must be the same size');

    % extract inputs
    nInputs = (numel(varargin) + 2);
    d = -1;
    epsilon = -1;
    epsFlag = 0;

    if nInputs == 2
        if numel(X)>1    
            d = 1;
        else
            error('you should not calculate RMSE for a singular')
        end
    elseif nInputs == 3
        if ( varargin{1} < numel(X) )
            d = varargin{1};
        else
            error('d cannot be larger than or equal to N');
        end
    elseif nInputs == 4
        if ( varargin{1} < numel(X) )
            d = varargin{1};
        else
            error('d cannot be larger than or equal to N');
        end
        if ( varargin{2} >= 0)
            epsilon = varargin{2};
        else
            error('epsilon should not be negative!');
        end
        epsFlag = 1;
    elseif (nInputs<2) || (nInputs>4)
        error('should be >1 and <5 input arguments')    
    end

    % Test RMSE == 0 when input vectors identical
    assert(Calc(1:10,1:10,1,0,0)==0,'RMSE of x=Y for equal vectors does not compute correctly');
    assert(Calc(ones(3,3,3),ones(3,3,3),1,0,0)==0,'RMSE of x=Y for equal matrices does not compute correctly');

    % Test RMSE calculated properly for dummy data
    assert(Calc([1 2 3; 1 2 3; 1 2 3],[3 2 1; 3 2 1; 3 2 1],1,0,0)==(sqrt(3)),'RMSE of x=Y for dummy data does not compute correctly');

    % Test RMSEeps calculated properly for dummy data
    assert(single(Calc([1 2 3],[1.4 2.4 3.4],1,0.3,1))==single(sqrt(0.015)),'RMSE epsilon insensitive for dummy data does not compute correctly');
    assert(single(Calc([1 2 3],[1.4 2.4 3.4],1,[0.3 0.4 0.4],1))==single(sqrt(0.005)),'RMSE epsilon insensitive for dummy data does not compute correctly');

    RMSE = Calc(X,y,d,epsilon,epsFlag);

end

function RMSE = Calc(X,y,d,epsilon,epsFlag)
%CALC Calculate RMSE & RMSE*

    % Calculate RMSE
    err = abs(X-y);
    if epsFlag
         err = max(err-epsilon,0);
    end
    RMSE = sqrt((1/(numel(err)-d)) * sum(err(:).^2));

end
