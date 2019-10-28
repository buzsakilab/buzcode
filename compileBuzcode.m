% this function recompiles all .c scripts to create the apropriate .mex* files on your local machine
%NOTE: current directory must be the buzcode base directory

buzpath = pwd;

addpath(genpath('externalPackages'))

compilefma % compiles FMAToolbox

try
cd(['externalPackages' filesep 'FilterM' filesep ''])
catch
    display('Please navigate to your local buzcode main directory')
end

if isunix | ismac
    mex('CFLAGS="\$CFLAGS -std=c99"', 'FilterX.c') % the above line fails with newer compilers but this works
    cd(['..' filesep '..'])
elseif ispc
    mex -O FilterX.c
    cd(['..' filesep '..'])  
end
try
cd(['externalPackages' filesep 'chronux_2_12' filesep 'locfit' filesep 'Source'])
compile % compiles chronux
cd(['..' filesep '..' filesep '..' filesep '..' filesep ''])
catch
    warning('CHRONUX DIDN''T COMPILE. sad.')
    cd(buzpath);
end

cd(['externalPackages' filesep 'xmltree-2.0' filesep '@xmltree' filesep 'private' filesep ''])
mex -O xml_findstr.c
cd(['..' filesep '..' filesep '..' filesep '..'])

cd(['analysis' filesep 'spikes' filesep 'correlation' filesep ''])
mex -O CCGHeart.c
cd(['..' filesep '..' filesep '..'])

% below is incomplete
% below adds buzcode to the matlab path and saves


% below adds the ' filesep 'buzcode' filesep 'generalComputation' filesep 'scripts' filesep ' folder to the system path
if isunix


end