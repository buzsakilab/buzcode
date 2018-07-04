% this script is a brief tutorial for mapping the firing rates of neurons
% onto the behavioral positions occupied by an animal
%
% This script assumes that you are currently in the
% /buzcode/tutorials/exampleDataStructs/20170505_396um_0um_merge folder, 
% it will not work otherwise!



sessionInfo = bz_getSessionInfo;
lfpChan = sessionInfo.ca1;
spikes = bz_GetSpikes('noprompts',true);
try
    lfp = bz_GetLFP(lfpChan);
catch
    error('couldnt load LFP, have you run the "download_DATA" script yet?')
end

behavior = bz_LoadBehavior(pwd,'track');

[firingMaps] = bz_firingMap1D(spikes,behavior,5,'savemat',false);
[phaseMaps] = bz_phaseMap1D(spikes,behavior,lfp,5,'savemat',false);




