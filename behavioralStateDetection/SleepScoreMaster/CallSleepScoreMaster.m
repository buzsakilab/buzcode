function CallSleepScoreMaster(basepath,basename,varargin)
% Lets you call SleepScoring by simply sitting in the path of the data to
% be sleepscored.  Also lets you call it from elsewhere with basepath and
% basename inputs.  Does not allow for SleepScoreMaster's multi-folder
% usage.  Does pass any additional arguments along to SleepScoreMaster...
% see optional inputs to SleepScoreMaster
%
% Brendon Watson 2016


if ~exist('basepath','var')
    basepath = cd;
    [~,basename] = fileparts(cd);
end

% if ~exist(fullfile(basepath,[basename '_SleepScore.mat']))
%     upperpath = fileparts(basepath);
    s = strfind(basepath,filesep);
    s(s==length(basepath)) = [];
    upperpath = basepath(1:s(end)-1);
    
    if isempty(varargin)
        SleepScoreMaster(upperpath,basename)
    else
        arginstring = '';
        for a = 1:length(varargin);
            eval(['v' num2str(a) ' = varargin{a};'])
            arginstring = strcat(arginstring,',v',num2str(a));
        end

        eval(['SleepScoreMaster(upperpath,basename' arginstring ');'])
    end
% else
%     disp([basename ' already SleepScored.  Rename _SleepScore.mat file if you want this to execute.'])
% end