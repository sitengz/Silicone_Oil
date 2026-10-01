function tests = test_film_surface_tension
% Run: results = runtests('test_film_surface_tension.m'); assertSuccess(results)
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = tempname; mkdir(root);
testCase.TestData.root = root;
makeCase(root,'test9',false);
makeCase(root,'test11',true);
end

function testNineColumnStatistics(testCase)
r = analyze_film_surface_tension('test9',testCase.TestData.root, ...
    fullfile(testCase.TestData.root,'results'));
verifyEqual(testCase,height(r.windows),16);
verifyEqual(testCase,r.summary.gamma_mN_m,14.5,'AbsTol',1e-10);
verifyEqual(testCase,r.summary.block_sem_mN_m,2.5,'AbsTol',1e-10);
verifyEqual(testCase,r.summary.block_count,2);
verifyEqual(testCase,r.summary.log_complete,true);
verifyEqual(testCase,r.summary.wall_diagnostics_prod,"absent");
verifyTrue(testCase,isfile(fullfile(r.output_dir,'test9_production_mean_errorbar.tif')));
end

function testSelectionAndTail(testCase)
r = analyze_film_surface_tension('test11',testCase.TestData.root, ...
    fullfile(testCase.TestData.root,'selected'), ...
    'ProductionStartNs',1,'ProductionEndNs',9,'BlockNs',3);
verifyEqual(testCase,r.summary.gamma_mN_m,13.5,'AbsTol',1e-10);
verifyEqual(testCase,r.summary.block_sem_mN_m,1.5,'AbsTol',1e-10);
verifyEqual(testCase,r.summary.excluded_partial_tail_ns,2);
verifyEqual(testCase,r.summary.wall_contacts_selected,0);
end

function testMissingFile(testCase)
verifyError(testCase,@() analyze_film_surface_tension('missing', ...
    testCase.TestData.root),'SurfaceTension:MissingFile');
end

function testBadCaseName(testCase)
verifyError(testCase,@() analyze_film_surface_tension('../escape', ...
    testCase.TestData.root),'SurfaceTension:CaseName');
end

function testDuplicateTimeRejected(testCase)
root = testCase.TestData.root;
makeCase(root,'duplicate',true);
path = fullfile(root,'energy.duplicate.film.dat');
text = fileread(path);
lines = splitlines(string(text));
lines(3) = lines(2);
fid = fopen(path,'wt'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',join(lines,newline));
clear cleanup;
verifyError(testCase,@() analyze_film_surface_tension('duplicate',root), ...
    'SurfaceTension:Clock');
end

function testBatchContinues(testCase)
[summary, failures] = analyze_film_surface_tension_batch( ...
    ["test11","missing"],testCase.TestData.root, ...
    fullfile(testCase.TestData.root,'batch'));
verifyEqual(testCase,height(summary),1);
verifyEqual(testCase,height(failures),1);
verifyEqual(testCase,failures.case_name,"missing");
verifyEqual(testCase,summary.gamma_mN_m,14.5,'AbsTol',1e-10);
end

function testAllFailuresWriteEmptySummary(testCase)
root = testCase.TestData.root;
[summary, failures] = analyze_film_surface_tension_batch("missing",root, ...
    fullfile(root,'all_failed'));
verifyEqual(testCase,height(summary),0);
verifyEqual(testCase,height(failures),1);
verifyTrue(testCase,isfile(fullfile(root,'all_failed','surface_tension_summary.csv')));
end

function makeCase(root,name,walls)
for phase = {'film_eq','film'}
    path = fullfile(root,['energy.' name '.' phase{1} '.dat']);
    fid = fopen(path,'wt'); cleanup = onCleanup(@() fclose(fid));
    header = '# time_fs temp_K pe_kcal_per_mol pxx_atm pyy_atm pzz_atm lx_A ly_A lz_A';
    if walls, header = [header ' wall_lo_force wall_hi_force']; end
    fprintf(fid,'%s\n',header);
    for k = 0:10
        gamma = 10;
        if strcmp(phase{1},'film'), gamma = gamma+k; end
        row = [k*1e6,300,-1000,2,2,2+gamma/0.506625,100,100,100];
        if walls, row = [row 0 0]; end %#ok<AGROW>
        fprintf(fid,'%.16g',row(1));
        fprintf(fid,' %.16g',row(2:end)); fprintf(fid,'\n');
    end
    clear cleanup;
end
fid = fopen(fullfile(root,['out.' name '.film']),'wt');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'Synthetic fixture\nTotal wall time: 0:00:01\n');
end
