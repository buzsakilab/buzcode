function [time_projection] = assembly_activity(AssemblyTemplates,SpikeCount)

% Activities = assembly_activity(Patterns,Activitymatrix): computes the
% time course of the activity of assembly patterns defined in Patterns in
% the spike matrix Activitymatrix with single bin resolution.
% 
% Description of inputs: 	
%   Activitymatrix: spike matrix. Rows represent neurons, columns represent
%   time bins. Patterns: assembly patterns. Columns denote assembly # and
%   rows neuron #.
% 
% Description of output: 	
%   Activities: Time course of the activity of assemblies. Rows represent
%   assemblies and columns represent time bins. 
%
% This framework is described in: Lopes-dos-Santos V, Ribeiro S, Tort ABL 
% (2013) Detecting cell assemblies in large neuronal populations, Journal
% of Neuroscience Methods.
%
% Please send bug reports to vitor@neuro.ufrn.br (V�tor)

zSpikeCount = zscore(SpikeCount')';

time_projection=zeros(size(AssemblyTemplates,2),size(zSpikeCount,2));
for assembly_idx = 1:size(AssemblyTemplates,2)
    
    % computing projector
    ASSEMBLYPROJECTOR=AssemblyTemplates(:,assembly_idx)*AssemblyTemplates(:,assembly_idx)';
    ASSEMBLYPROJECTOR=squeeze(ASSEMBLYPROJECTOR)-diag(diag(squeeze(ASSEMBLYPROJECTOR)));
    
    % computing activity time course
    for ntime=1:size(zSpikeCount,2)
        
        time_projection(assembly_idx,ntime)=(zSpikeCount(:,ntime)'*ASSEMBLYPROJECTOR*zSpikeCount(:,ntime));
        
    end
    
end