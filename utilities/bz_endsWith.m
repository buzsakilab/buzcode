function tf = bz_endsWith(s, pattern, varargin)
%ENDSWITH True if text ends with pattern.
%   TF = endsWith(S,PATTERN) returns true if any element of string array S 
%   ends with PATTERN. TF is the same size as S.
%
%   S can be a string array, a character vector, or a cell array of
%   character vectors. So can PATTERN. PATTERN and S need not be the same
%   size. If PATTERN is nonscalar, endsWith returns true if S ends with any
%   element of PATTERN.
%
%   TF = endsWith(S,PATTERN,'IgnoreCase',IGNORE) ignores case when searching 
%   for PATTERN at the end of S if IGNORE is true. The default value of IGNORE 
%   is false.
%
%   Examples
%       S = string('data.tar.gz');
%       P = string('gz');
%       endsWith(S,P)                   returns  1
%
%       S = string({'abstracts.docx','data.tar.gz'});
%       P = 'docx';         
%       endsWith(S,P)                   returns  [1 0]
%
%       S = string('abstracts.docx');
%       P = {'docx','tar.gz'};
%       endsWith(S,P)                   returns  1
%
%       S = string({'DATA.TAR.GZ','SUMMARY.PPT'});
%       P = string('ppt');
%       endsWith(S,P,'IgnoreCase',true) returns  [0 1]
%
%   See also startsWith, contains.
 
%   Copyright 2015-2016 The MathWorks, Inc.
 
    narginchk(2, inf);
 
    if ~ischar(s) && ~iscellstr(s) && ~isstring(s)
        firstInput = getString(message('MATLAB:string:FirstInput'));
        error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    end
 
    try
        stringS = string(s);
        tf = stringS.endsWith(pattern, varargin{:});
    catch E
        throw(E)
    end
end

