dataDir = '/media/sdb6/Data';

cd(dataDir);

A = Analysis(pwd);

datasets = List2Cell('datasets.list');

!son2Mat.py --events 