function  [result] = verify_mexplex(mexPlexData, dataDir)
result = 1;

for i=1:length(mexPlexData.plxs)
    plx = mexPlexData.plxs{i};
    fpath = fullfile(dataDir, plx.FileName);
    plxnew = get_all_from_plx(fpath);
    [r1 r2 er] = comp_struct(plx, plxnew,'old','new');
    if(max(size(er)) > 0) result = 0; end
end
ddtData = mexPlexData.ddt;
fpath = fullfile(dataDir, ddtData.FileName);
ddtNew = get_all_from_ddt(fpath);
[r1 r2 er] = comp_struct(ddtData, ddtNew,'old','new');
if(max(size(er)) > 0) result = 0; end

bad_par_test = plx_test_bad_parameters();
if bad_par_test.allTestsPassed == 0
	result = 0;
end
file_test = plx_test_file_variations();
if file_test.allTestsPassed == 0
	result = 0;
end

if result == 0
    disp('VERIFICATION FAILED');
else
    disp('VERIFICATION TESTS PASSED');
end
	