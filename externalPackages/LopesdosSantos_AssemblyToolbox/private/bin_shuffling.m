function control_max_eig = bin_shuffling(SpikeCount,number_of_surrogates)

control_max_eig=zeros(1,number_of_surrogates);
for surr_idx = 1:number_of_surrogates
    ControlSpikeCount=zeros(size(SpikeCount));
    for neuron_idx = 1:size(SpikeCount,1)
        ControlSpikeCount(neuron_idx,:) = SpikeCount(neuron_idx,randperm(length(SpikeCount)));
    end
    CorrelationMatrix = corr(ControlSpikeCount');
    [~,eigenvalues]=eig(CorrelationMatrix);
    control_max_eig(surr_idx)=max(reshape(eigenvalues,1,[]));
end