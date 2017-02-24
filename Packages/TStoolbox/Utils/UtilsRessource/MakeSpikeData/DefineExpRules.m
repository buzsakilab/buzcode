function A = DefineExpRules(A)

% user is asked to define the different tasks of each session. It is saved in trialrules.txt. It is then saved in the database in the function DefCorrectError with the correctError vector.


A = getResource(A,'StartTrial');
dset = current_dataset(A);
fileRules = [current_dir(A) filesep 'trialRules.txt']

doitagain = 'Y';

if exist(fileRules)
	doitagain = questdlg('File already exist, do it again?');
end

if doitagain=='Y'
	
	nbTrials = length(Data(startTrial{1}))
	
	trialRules = zeros(nbTrials,1);
	shiftTrial = 1;
	badAnswer = 1;
	
	while badAnswer
		
		goAhead = 1;
		trialRules = zeros(nbTrials,1);
		shiftTrial = 1;
			
		while goAhead
			isShift = [];	
	
			rule = inputdlg(['Startegy ' num2str(length(shiftTrial)) ' (1=right, 2=light,3=left and 4=dark) for ' dset]);
			rule=str2num(rule{1});
			while (length(isShift)==0)
				isShift = questdlg('Is there a strategy shift?');
			end;
			
			if isShift(1)=='N'
				for ix=shiftTrial(end):nbTrials
					trialRules(ix) = rule;
				end;
				goAhead = 0;
	
			elseif isShift(1)=='Y'
				shift = inputdlg('First trial of next strategy :');
				shift = str2num(shift{1})
				if shift<nbTrials
					shiftTrial = [shiftTrial shift];	
					for ix=shiftTrial(end-1):(shiftTrial(end)-1)
						trialRules(ix) = rule;
					end;
				else 
					display('error! shift > nbTrials')
					goAhead = 0;
				end;
			else
				diplay('Hmmm let''s start again')
				goAhead = 0;
			end;
		
		end;
	
		display('shifts :');
		display(shiftTrial);
		trialRules = trialRules
		isFinish = questdlg('Is this correct?')
		if isFinish(1)=='Y'
			badAnswer = 0;
		end;	
	
	end;
	
	save(fileRules,'trialRules','-ascii');

end