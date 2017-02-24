function WriteHeader(fp, varargin)

% WriteHeader(fp, H1, H2, H3, ...)
%
% INPUTS
%    fp = file pointer
%    H1, H2, H3, ... = lines (strings) to write out as header
% 
% OUTPUTS: none
%  
% Writes NSMA header
%
% ADR 1998
% version U3.1
% status PROMOTED

global MCLUST_VERSION

% v 3.1 now accepts cell arrays as well as strings
fprintf(fp, '%%%%BEGINHEADER\n');
fprintf(fp, '%% Program: matlab\n');
fprintf(fp, '%% MClust version: %s\n', MCLUST_VERSION);
fprintf(fp, '%% Date: %s\n', datestr(now));
fprintf(fp, '%% Directory: %s\n', pwd);

if ~isempty(getenv('HOST'))
   fprintf(fp, [ '%% Hostname: ', getenv('HOST'), '\n']);

end

if ~isempty(getenv('USER'))
   fprintf(fp, [ '%% User: ', getenv('USER'), '\n']);

end

for iH = 1:length(varargin)

   if isa(varargin{iH}, 'cell')

      for jH = 1:length(varargin{iH})

         fprintf(fp, '%% %s\n', varargin{iH}{jH});

      end

   elseif isa(varargin{iH}, 'char')
      fprintf(fp, '%% %s\n', varargin{iH});

   else

      error('Unknown input type.');

   end
end
fprintf(fp, '%%%%ENDHEADER\n');
