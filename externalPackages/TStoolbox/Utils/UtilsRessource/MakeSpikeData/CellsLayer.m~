function A = truc(A)

% read a csv file given by histologist to put in DB the kayers
% Strucutre : 1:PrL, 2:PrL/IL, 3:IL/ 
% Layer : 1:superficial, 2 : middle, 3 : deep.


A = getResource(A,'SpikeData');

A = registerResource(A, 'Layer', 'numeric', {[], []}, ...
    'layer', ['cells layer : 1:superficial, 2 : middle, 3 : deep.']);

A = registerResource(A, 'Structure', 'numeric', {[], []}, ...
    'structure', ['structures, 1:PrL, 2:PrL/IL, 3:IL']);


d = dlmread([parent_dir(A) filesep 'CellsLayer.csv']);
days = d(:,1);
Structure = d(:,4);
Layer = d(:,5);

[p,dataset,e] = fileparts(current_dir(A));
t = find(days == str2num(dataset));

if length(t)==length(S);

	layer = Layer(t);
	structure = Structure(t);

else

	keyboard
	layer = zeros(length(S),1);
	structure = zeros(length(S),1);

end

A = saveAllResources(A);