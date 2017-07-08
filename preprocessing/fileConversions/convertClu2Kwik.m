function [] = convertClu2Kwik(varargin)
% USAGE
% [] = restoreOriginalKwik(kwikfile)
%
% INPUTS
% 
%   kwikfile - filename of .kwik file in current working directory 
% 
% OUTPUT
%
% this function restores the original klustakwik results to the /main branch of the .kwik file
% David Tingley, 2017


p = inputParser;
addParameter(p,'tkwik','',@isstr)
addParameter(p,'clu','',@isvector)
parse(p,varargin{:})
tkwik = p.Results.tkwik;
clu = p.Results.clu;

disp(tkwik)

if isempty(tkwik)  % remove this? we shouldn't assume there is a single .kwik in the cwd..
    try
    tkwik = dir('*.kwik');
    tkwik = tkwik.name;
    catch
    end
end
cluster_groups = unique(clu);

if ~isempty(tkwik)
    elec = split(tkwik,'.');
    elec = split(elec{1},'_sh');
    elec = str2num(elec{end});  % assuming you break the .kwik files up by seperate shanks like I do..

    kwikinfo_original = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/original']);
    kwikinfo_main = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/main']);
    % first, let's clear out the /main branch
    fid = H5F.open(kwikinfo_main.Filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    for gg = 1:length(kwikinfo_main.Groups)
        H5L.delete(fid,[kwikinfo_main.Groups(gg).Name],'H5P_DEFAULT')
    end
    H5F.close(fid)
    % now re-populate from /original branch..
    fid = H5F.open(kwikinfo_original.Filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    for g = 1:length(cluster_groups)
%         clust = split(kwikinfo_original.Groups(g).Name,'/');
%         clust = clust{end};
        clust = cluster_groups(g);
        try
        H5L.copy(fid,kwikinfo_original.Groups(g).Name,...
                        fid,[kwikinfo_main.Name '/' num2str(clust)],...
                        'H5P_DEFAULT','H5P_DEFAULT')
        catch
            disp([ clust ' cluster group already exists'])
        end
    end
    for g= 1:length(cluster_groups)
        if g > 2
        h5writeatt(tkwik,['/channel_groups/' num2str(elec) '/clusters/main/' ...
            num2str(cluster_groups(g))],'cluster_group',2) 
        end
    end
    H5F.close(fid)
    % above we restore cluster groups, but not spike ID
    % So now let's restore cluster ID's
    try
        clu_orig = h5read(tkwik,['/channel_groups/' num2str(elec) '/spikes/clusters/original']);
    catch
        error(['could not restore cluster IDs for ' tkiwk])
    end
    h5write(tkwik,['/channel_groups/' num2str(elec) '/spikes/clusters/main'],uint32(clu));
    disp(['added clu file cluster IDs for ' tkwik])
end
