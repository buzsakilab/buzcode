 function nwb = addElectrodeInfo (xml,nwb)
    % Adds electrode info in: nwb.general_extracellular_ephys_electrodes
    % Konstantinos Nasiotis 2019
    
    %%
    nShanks = length(xml.spikeDetection.channelGroups.group);

    groups = xml.spikeDetection.channelGroups.group; % Use this for simplicity

    all_shank_channels = cell(nShanks,1); % This will hold the channel numbers that belong in each shank            

    % Initialize variables
    x                 = [];
    y                 = [];
    z                 = [];
    imp               = [];
    location          = [];
    shank             = [];
    group_name        = [];
    group_object_view = [];
    filtering         = [];
    shank_channel     = [];
    amp_channel_id    = [];

    device_name = 'implant';
    nwb.general_devices.set(device_name, types.core.Device());
    device_link = types.untyped.SoftLink(['/general/devices/' device_name]);

    for iGroup = 1:nShanks
        for iChannel = 1:length(groups{iGroup}.channels.channel)
            all_shank_channels{iGroup} = [all_shank_channels{iGroup} str2double(groups{iGroup}.channels.channel{iChannel}.Text)];
            shank_channel     = [shank_channel; iChannel-1];
            amp_channel_id    = [amp_channel_id; str2double(groups{iGroup}.channels.channel{iChannel}.Text)];
            shank             = [shank; iGroup];

            if nShanks > 9 && iGroup<10
                group_name = [group_name; 'shank' num2str(iGroup) ' '];
            else
                group_name = [group_name; 'shank' num2str(iGroup)];
            end

            group_object_view = [group_object_view; types.untyped.ObjectView(['/general/extracellular_ephys/' ['shank' num2str(iGroup)]])];

            if ~isfield(groups{iGroup}.channels.channel{iChannel},'position')
                x = [x; NaN];
                y = [y; NaN];
                z = [z; NaN];
            end
            if ~isfield(groups{iGroup}.channels.channel{iChannel},'imp')
                imp = [imp; NaN];
            end  
            if ~isfield(groups{iGroup}.channels.channel{iChannel},'location')
                location{end+1,1} = 'unknown';
            end  
            if ~isfield(groups{iGroup}.channels.channel{iChannel},'filtering')
                filtering = [filtering; NaN];
            end      

        end
        nwb.general_extracellular_ephys.set(['shank' num2str(iGroup)], ...
            types.core.ElectrodeGroup( ...
            'description', ['electrode group for shank' num2str(iGroup)], ...
            'location', 'unknown', ...
            'device', device_link));

    end

    variables = {'x'; 'y'; 'z'; 'imp'; 'location'; 'filtering'; 'group'; 'group_name'; 'shank'; 'shank_channel'; 'amp_channel'};

    % In order to insert string to a table, they need to be converted to a cell
    % first (e.g. location(iElectrode))
    for iElectrode = 1:length(x)
        if iElectrode == 1
            tbl = table(x(iElectrode),y(iElectrode),z(iElectrode),imp(iElectrode),{location{iElectrode}},filtering(iElectrode),group_object_view(iElectrode),{group_name(iElectrode,:)},shank(iElectrode),shank_channel(iElectrode),amp_channel_id(iElectrode),...
                       'VariableNames', variables);
        else
            tbl = [tbl; {x(iElectrode),y(iElectrode),z(iElectrode),imp(iElectrode),{location{iElectrode}},filtering(iElectrode),group_object_view(iElectrode),{group_name(iElectrode,:)},shank(iElectrode),shank_channel(iElectrode),amp_channel_id(iElectrode)}];
        end
    end

    % add the |DynamicTable| object to the NWB file in
    % /general/extracellular_ephys/electrodes
    electrode_table = util.table2nwb(tbl, 'metadata about extracellular electrodes');
    nwb.general_extracellular_ephys_electrodes = electrode_table;

    disp('Electrode info added..')

end