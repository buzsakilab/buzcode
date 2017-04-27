function binnedmatrix = toy_simulation(Network_opts,Assembly_opts)

% creating matrix
binnedmatrix=poissrnd(Network_opts.meanspikebin,[Network_opts.nneurons Network_opts.nbins]); 

allassemblyneurons = unique([Assembly_opts.assembly_neurons{:}]);
nonassemblyneurons = setdiff(1:Network_opts.nneurons,allassemblyneurons); 
nassemblies = length(Assembly_opts.assembly_neurons);

% add activations of all assemblies in random bins
for assemblyindex=1:nassemblies
    
    assemblysize = length(Assembly_opts.assembly_neurons{assemblyindex});
    
    % drawing activations
    activation_bins = randi([1,Network_opts.nbins],1,Assembly_opts.number_of_activations);
    
    % introducting activations
    binnedmatrix(Assembly_opts.assembly_neurons{assemblyindex},activation_bins) = ...
        poissrnd(Assembly_opts.meanspikerate_activations,...
        assemblysize,Assembly_opts.number_of_activations);
    
end
