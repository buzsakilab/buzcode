function A = makeSpikeData(A) 

maxTT = 6;

A = registerResource(A, 'CellNames', 'cell', {[], 1}, ...
    'cellnames', ...
    'name of cells ', 'mfile');

A = registerResource(A, 'SpikeData', 'tsdArray', {[], 1}, ...
    'S', ...
    'spike trains', 'mfile');

A = registerResource(A, 'Tetrode', 'numeric', {[], 1}, ...
    'TT', ...
    'tetrode from whcih the cell   was recorded', 'mfile');

[dummy,dsetPrefipx,dummy1] = fileparts(current_dir(A));

S = {};
cellnames = {};
TT = [];

for t = 1:maxTT 
    cluFile = [current_dir(A) filesep dsetPrefipx '.clu.' num2str(t)];
    timFile = [current_dir(A) filesep dsetPrefipx '.tim.' num2str(t)];
    
    if exist(cluFile, 'file')
        clu = load(cluFile);
        cellIdx = 2:(clu(1)-1);
        clu = clu(2:end);
        tim = load(timFile);
        tim = tim(:,2);
        for i = cellIdx
            cn = ['TT' num2str(t) '_c' num2str(i)];
            cellnames = [ cellnames ; { cn } ];
            TT = [TT ; t];
            sp = ts(tim(find(clu==i)));
            S = [ S ; { sp } ];
        end

    end
end



S = tsdArray(S);

A = saveAllResources(A);
