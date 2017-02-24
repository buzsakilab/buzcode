dataDir = '/media/sdc6/Data';
cd(dataDir)
A = Analysis(dataDir);
datasets = List2Cell([ dataDir filesep 'datasets_noYMPpb.list' ] );


for day=51:length(datasets)

	fprintf(['Day : ' datasets{day} ' -- day # ' num2str(day) '\n']);
	answer = 'No';
	cd(datasets{day});
	while answer(1)~='Y'	
		YMProcessPosFile
		answer = questDlg('Satisfied?');
	end;
	cd('../..');

end