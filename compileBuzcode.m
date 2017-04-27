% this function recompiles all .c scripts to create the apropriate .mex* files on your local machine
%NOTE: current directory must be the buzcode base directory

addpath(genpath('externalPackages'))

compilefma % compiles FMAToolbox

try
cd('externalPackages/FilterM/')
catch
    display('Please navigate to your local buzcode main directory')
end

% mex -O FilterX.c % compiles FilterM (faster filtering than filtfilt)
mex('CFLAGS="\$CFLAGS -std=c99"', 'FilterX.c') % the above line fails with newer compilers but this works
cd('../..')
    

cd('externalPackages/chronux_2_12/locfit/Source')
compile % compiles chronux
cd('../../../../')

cd('externalPackages/xmltree-2.0/@xmltree/private/')
mex -O xml_findstr.c
cd('../../../..')

cd('analysis/spikes/correlation/')
mex -O CCGHeart.c
cd('../../..')

% below is incomplete
% below adds buzcode to the matlab path and saves


% below adds the /buzcode/generalComputation/scripts/ folder to the system path
if isunix


end