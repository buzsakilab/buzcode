% this function recompiles all .c scripts to create the apropriate .mex* files on your local machine
%NOTE: current directory must be the buzcode base directory
    
%try
%cd('externalPackages/chronux_2_12/locfit/Source')
%compile % compiles chronux
%cd('../../../../')
%catch
%    warning('CHRONUX DIDN''T COMPILE. sad.')
%end


cd('analysis/spikes/correlation/')
mex -O CCGHeart.c
cd('../../..')