% this function recompiles all .c scripts to create the apropriate .mex* files on your local machine

addpath(genpath('externalPackages'))

compilefma % compiles FMAToolbox

cd('externalPackages/FilterM/')
% mex -O FilterX.c % compiles FilterM (faster filtering than filtfilt)
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

