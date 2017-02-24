function A = truc(A)

% Parameters
spkWidthLim = 0.42; %in ms, obtained from EM_GM on the whole dataset

A = getResource(A,'PWidth');

A = registerResource(A, 'CellType', 'numeric', {[], []}, ...
    'cellType', ['IN=1, PC=2, unknown=0']);


cellType = zeros(size(pWidth));
cellType(pWidth>0) = 1;
cellType(pWidth>spkWidthLim) = 2;

A = saveAllResources(A);