function [ddtData] = get_all_from_ddt(filename)

ddtData  =[];
[dirstr, name, ext] = fileparts(filename);
ddtData.FileName = [name ext];

[ddtData.raw.nch, ddtData.raw.npoints, ddtData.raw.freq, ddtData.raw.d] = ddt(filename);
[ddtData.v.nch, ddtData.v.npoints, ddtData.v.freq, ddtData.v.d] = ddt_v(filename);