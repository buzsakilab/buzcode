%  A set of defaults used for FD file creation and other processes in
%  multiple locations of the MClust code
%
%
 
% ncst 03 Dec 03

global MClust_max_records_to_load MClust_ChannelValidity MClust_TTdn MClust_TTfn MClust_TTdn

if ~isempty(MClust_max_records_to_load)
	record_block_size = MClust_max_records_to_load;
else
	record_block_size = 80000;  % maximum number of spikes to load into memory at one time
end

template_matching = 0; % used to remove noise spikes which are not "spike-like,", template-matching is not a currently supported function

NormalizeFDYN = 'no'; % 'yes' if you want to normalize the feature data files to mean = 0, std = 1