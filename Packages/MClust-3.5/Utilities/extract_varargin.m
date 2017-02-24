% Extract_varargin

%   NOT A FUNCTION -- this allows it to access the current workspace
%
%   expects varargin to consist of sequences of 'variable', value
%   sets variable to value for each pair.
%   changes the current workspace!

% ADR 1998
% version L4.0
% status: PROMOTED

for iV = 1:2:length(varargin)
  if size(varargin{iV},2) == 1
    varargin{iV} = varargin{iV}';
  end
  eval([varargin{iV}, ' = ', 'varargin{iV+1};']);
end