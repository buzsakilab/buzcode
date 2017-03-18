function [] = restoreOriginalKwik()

tkwik = dir('*.kwik');
if ~isempty(tkwik)
    tkwik = tkwik.name;
    elec = split(tkwik,'.');
    elec = split(elec{1},'_sh');
    elec = str2num(elec{end});

    kwikinfo_original = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/original']);
    kwikinfo_main = h5info(tkwik,['/channel_groups/' num2str(elec) '/clusters/main']);

    fid = H5F.open(kwikinfo_main.Filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    for gg = 1:length(kwikinfo_main.Groups)
        H5L.delete(fid,[kwikinfo_main.Groups(gg).Name],'H5P_DEFAULT')
    end
    H5F.close(fid)

    fid = H5F.open(kwikinfo_original.Filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    for g = 1:length(kwikinfo_original.Groups)
        clust = split(kwikinfo_original.Groups(g).Name,'/');
        clust = clust{end};
        try
        H5L.copy(fid,kwikinfo_original.Groups(g).Name,...
                        fid,[kwikinfo_main.Name '/' clust],...
                        'H5P_DEFAULT','H5P_DEFAULT')
        catch
            disp([ clust ' cluster group already exists'])
        end
    end
    H5F.close(fid)
    % above we restore cluster groups, but not spike ID
    basename = split(tkwik,'_sh');
    basename = basename{1};
    restoreClu(pwd,basename,elec); % restores cluster ID's
end
