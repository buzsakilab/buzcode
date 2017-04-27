clear all % just clearing workspace

% define mean firing rate
Network_opts.meanspikebin = 1;

% define number of neurons
Network_opts.nneurons = 20 ;

% define number of bins
Network_opts.nbins = 10000;

% above define assembly membership

% line below sets neurons 1,2,3 and 4 to be in assembly 1
Assembly_opts.assembly_neurons{1} = [1 2 3 4]; 
% line below sets neurons 5,6 and 7 to be in assembly 2
Assembly_opts.assembly_neurons{2} = [5 6 7]; 

% defines number of activation bins
Assembly_opts.number_of_activations = 300;

% defines mean rate in activation bins
Assembly_opts.meanspikerate_activations = 3;

% running the function
Activitymatrix = toy_simulation(Network_opts,Assembly_opts);

%%

correlationmat = corr(Activitymatrix');
figure(1),clf
imagesc(correlationmat)

%%

opts.threshold.permutations_percentile = 95;
opts.threshold.number_of_permutations = 20;
opts.threshold.method = 'circularshift';
opts.AssemblyTemplate.method = 'PCA';
AssemblyTemplates = assembly_patterns(Activitymatrix,opts);
 
%%
figure(3),clf
subplot(211)
stem(AssemblyTemplates(:,1))
subplot(212)
stem(AssemblyTemplates(:,2))

%%

clear all % just clearing workspace

Network_opts.meanspikebin = 1;
Network_opts.nneurons = 20;
Network_opts.nbins = 10000;

Assembly_opts.assembly_neurons{1} = [1 2 3 4 5]; 
Assembly_opts.assembly_neurons{2} = [4 5 6 7 8]; 
Assembly_opts.number_of_activations = 300;
Assembly_opts.meanspikerate_activations = 3;
Activitymatrix = toy_simulation(Network_opts,Assembly_opts);

opts.threshold.method = 'MarcenkoPastur';
opts.AssemblyTemplate.method = 'PCA';
AssemblyTemplates = assembly_patterns(Activitymatrix,opts);

figure(2),clf
subplot(211)
stem(AssemblyTemplates(:,1))
subplot(212)
stem(AssemblyTemplates(:,2))


%%
opts.threshold.method = 'MarcenkoPastur';
opts.AssemblyTemplate.method = 'ICA';
opts.AssemblyTemplate.number_of_iterations = 200;
AssemblyTemplates = assembly_patterns(Activitymatrix,opts);


figure(3),clf
subplot(211)
stem(AssemblyTemplates(:,1))
subplot(212)
stem(AssemblyTemplates(:,2))

%%

Activities = assembly_activity(AssemblyTemplates,Activitymatrix);

figure(4),clf
subplot(211)
imagesc(Activitymatrix)
xlim([0 100])
subplot(212)
plot(Activities')
xlim([0 100])




    