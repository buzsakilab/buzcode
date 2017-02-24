function Make_TSToolbox()

objDirs = dir('@*');
dirNames = {'@intervalSet';'@tsd'};
parent_dir = pwd;

for ii=1:length(dirNames)
   fprintf('In %s \n',dirNames{ii})
   cd([dirNames{ii} filesep 'private'])
   cfiles = dir('*.c');
   for jj=1:length(cfiles)
       fprintf('Compiling %s \n',cfiles(jj).name)
       mex(cfiles(jj).name);
   end
   cd(parent_dir)
end
