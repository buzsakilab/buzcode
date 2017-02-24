function SaveRes(FileName,times)

BufSize = inf;

formatstring = '%d';
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');

fprintf(outputfile,formatstring,round(times));

fclose(outputfile);
