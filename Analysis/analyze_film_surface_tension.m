function result = analyze_film_surface_tension(caseName, dataDir, outputDir, varargin)
%ANALYZE_FILM_SURFACE_TENSION Standard film pressure-tensor analysis.
% result = analyze_film_surface_tension('N4_PDI1', dataDir, outputDir)
% Optional name/value settings: WindowNs (5), ShiftNs (1), BlockNs (5),
% ProductionStartNs (0), ProductionEndNs (Inf). Selection times are relative
% to the production clock, not the concatenated equilibration/production axis.
% Error bars are block SEM, not SD or the spread of overlapping windows.
% No toolbox beyond MATLAB is required. Inputs are never modified.

caseName = char(string(caseName));
if isempty(regexp(caseName, '^[A-Za-z0-9][A-Za-z0-9_.-]*$', 'once'))
    error('SurfaceTension:CaseName', 'Use a case name without path separators.');
end
if nargin < 2 || isempty(dataDir), dataDir = pwd; end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
positiveScalar = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0;
p = inputParser;
addParameter(p, 'WindowNs', 5, positiveScalar);
addParameter(p, 'ShiftNs', 1, positiveScalar);
addParameter(p, 'BlockNs', 5, positiveScalar);
addParameter(p, 'ProductionStartNs', 0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ProductionEndNs', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && ~isnan(x) && x > 0);
parse(p, varargin{:});
settings = p.Results;
eqPath = fullfile(dataDir, ['energy.' caseName '.film_eq.dat']);
prodPath = fullfile(dataDir, ['energy.' caseName '.film.dat']);
logPath = fullfile(dataDir, ['out.' caseName '.film']);
eq = readEnergy(eqPath);
prod = readEnergy(prodPath);
eqTime = eq.time_fs * 1e-6;
prodTime = prod.time_fs * 1e-6;
if eqTime(1) ~= 0 || prodTime(1) ~= 0
    error('SurfaceTension:Clock', 'Both energy files must start at time zero.');
end
eqEnd = eqTime(end);
prodEnd = prodTime(end);
boundaryNames = {'pxx_atm','pyy_atm','pzz_atm','lx_A','ly_A','lz_A'};
if any(abs(eq{end,boundaryNames} - prod{1,boundaryNames}) > 1e-6)
    error('SurfaceTension:Boundary', ...
        'Stage boundary snapshots differ; check that the files belong to the same run.');
end
selectedStart = settings.ProductionStartNs;
selectedEnd = min(settings.ProductionEndNs, prodEnd);
if selectedEnd <= selectedStart
    error('SurfaceTension:Interval', 'The requested production interval is unavailable.');
end
if isfinite(settings.ProductionEndNs) && settings.ProductionEndNs > prodEnd
    warning('SurfaceTension:ShortRun', 'Requested end exceeds available production; using %g ns.', prodEnd);
end

% Two interfaces; global box pressure paired with the same global box Lz.
% LAMMPS real units: atm * angstrom -> mN/m.
eqGamma = calculateGamma(eq);
prodGamma = calculateGamma(prod);
time = [eqTime(1:end-1); eqEnd + prodTime];
gamma = [eqGamma(1:end-1); prodGamma];
windows = windowMeans(time, gamma, 0, time(end), settings.WindowNs, settings.ShiftNs);
if isempty(windows)
    error('SurfaceTension:ShortRun', 'No complete sliding window is available.');
end
windows.stage = repmat("straddles boundary", height(windows), 1);
windows.stage(windows.end_ns <= eqEnd) = "equilibration";
windows.stage(windows.start_ns >= eqEnd) = "production";
blocks = windowMeans(prodTime, prodGamma, selectedStart, selectedEnd, ...
    settings.BlockNs, settings.BlockNs);
blocks.combined_start_ns = blocks.start_ns + eqEnd;
blocks.combined_end_ns = blocks.end_ns + eqEnd;
stats = blockStatistics(blocks.gamma_mN_m);
if stats.count < 2
    warning('SurfaceTension:TooFewBlocks', 'Fewer than two complete blocks: SEM and CI are unavailable.');
end
selected = prodTime >= selectedStart & prodTime < selectedEnd;
allProduction = prodTime < prodEnd;
lateStart = max(selectedStart, selectedEnd - 20);
lateBlocks = windowMeans(prodTime, prodGamma, lateStart, selectedEnd, ...
    settings.BlockNs, settings.BlockNs);
lateStats = blockStatistics(lateBlocks.gamma_mN_m);
nHalf = floor(height(blocks) / 2);
firstHalf = NaN; lastHalf = NaN; slope = NaN;
if nHalf > 0
    firstHalf = mean(blocks.gamma_mN_m(1:nHalf));
    lastHalf = mean(blocks.gamma_mN_m(end-nHalf+1:end));
end
if height(blocks) >= 2
    coefficients = polyfit(blocks.center_ns, blocks.gamma_mN_m, 1);
    slope = coefficients(1);
end
[eqContacts, eqWallStatus] = wallContacts(eq);
[prodContacts, prodWallStatus] = wallContacts(prod);
selectedContacts = NaN;
if prodWallStatus == "available"
    selectedContacts = nnz(selected & ...
        (prod.wall_lo_force ~= 0 | prod.wall_hi_force ~= 0));
end
if selectedContacts > 0
    warning('SurfaceTension:WallContact', ...
        '%d selected production samples have wall contact; this is not an unconfined-film estimate.', selectedContacts);
end
logPresent = isfile(logPath); logComplete = false; logError = false; logNewer = false;
if logPresent
    logText = fileread(logPath);
    logComplete = contains(logText, 'Total wall time:');
    logError = ~isempty(regexp(logText, '(?m)^\s*ERROR:', 'once'));
    logInfo = dir(logPath); eqInfo = dir(eqPath); prodInfo = dir(prodPath);
    logNewer = ~logComplete && ...
        logInfo.datenum > max(eqInfo.datenum, prodInfo.datenum) + 1/86400;
end
if ~logComplete || logError
    warning('SurfaceTension:Log', ...
        'The matching log is missing, incomplete, or contains ERROR; verify job status before calling this final.');
end
if logNewer && ~logComplete
    warning('SurfaceTension:FileAge', ...
        'The incomplete log is newer than both energy files; these may be saved data from an earlier run.');
end
if prodWallStatus == "absent"
    warning('SurfaceTension:WallsUnknown', ...
        'Production has no wall-force columns; check the input/log or film geometry manually.');
end

summary = table(string(caseName), eqEnd, prodEnd, selectedStart, selectedEnd, ...
    settings.WindowNs, settings.ShiftNs, settings.BlockNs, ...
    mean(prodGamma(allProduction)), mean(prodGamma(selected)), ...
    stats.mean, stats.sem, stats.sd, stats.count, stats.ci(1), stats.ci(2), ...
    selectedEnd-selectedStart-stats.count*settings.BlockNs, ...
    lateStart, lateStats.mean, lateStats.sem, lateStats.count, ...
    firstHalf, lastHalf, lastHalf-firstHalf, slope, ...
    eqContacts, prodContacts, selectedContacts, eqWallStatus, prodWallStatus, ...
    logPresent, logComplete, logError, logNewer, string(eqPath), string(prodPath), string(logPath), ...
    'VariableNames', {'case_name','equil_duration_ns','production_duration_ns', ...
    'selected_start_ns','selected_end_ns','window_ns','shift_ns','block_ns', ...
    'all_production_mean_mN_m','selected_sample_mean_mN_m','gamma_mN_m', ...
    'block_sem_mN_m','block_sd_mN_m','block_count','ci95_low_mN_m','ci95_high_mN_m', ...
    'excluded_partial_tail_ns','late_start_ns','late_mean_mN_m','late_sem_mN_m', ...
    'late_block_count','first_half_mean_mN_m','last_half_mean_mN_m', ...
    'half_difference_mN_m','block_slope_mN_m_per_ns','wall_contacts_eq', ...
    'wall_contacts_prod','wall_contacts_selected','wall_diagnostics_eq', ...
    'wall_diagnostics_prod','log_present','log_complete','log_has_error', ...
    'incomplete_log_newer_than_energy','eq_file','production_file','log_file'});
summary.review_status = "REVIEW trend, block-size sensitivity, and film geometry before accepting";
caseDir = fullfile(outputDir, caseName);
if ~isfolder(caseDir), mkdir(caseDir); end
stem = fullfile(caseDir, caseName);
writetable(windows, [stem '_windows.csv']);
writetable(blocks, [stem '_blocks.csv']);
writetable(lateBlocks, [stem '_late_blocks.csv']);
writetable(summary, [stem '_summary.csv']);
plotTrend(windows, [0 time(end)], eqEnd, [stem '_sliding_surface_tension']);
prodWindows = windows(windows.start_ns >= eqEnd, :);
if ~isempty(prodWindows)
    plotTrend(prodWindows, [eqEnd time(end)], [], [stem '_production_surface_tension']);
end
fig = newFigure(); cleanup = onCleanup(@() close(fig));
ax = axes(fig);
if isfinite(stats.sem)
    errorbar(ax, 1, stats.mean, stats.sem, 'o', 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410], 'CapSize', 12);
else
    plot(ax, 1, stats.mean, 'o', 'MarkerSize', 8);
end
styleAxes(ax);
set(ax, 'XTick', 1, 'XTickLabel', {caseName}, 'TickLabelInterpreter', 'none');
xlim(ax, [0.5 1.5]);
saveFigure(fig, [stem '_production_mean_errorbar']);
result.summary = summary;
result.windows = windows;
result.blocks = blocks;
result.late_blocks = lateBlocks;
result.output_dir = caseDir;
fprintf('\n%s: %.4f +/- %.4f mN/m (block SEM; %d x %g ns blocks).\n', ...
    caseName, stats.mean, stats.sem, stats.count, settings.BlockNs);
fprintf('Selected production: %g-%g ns; 95%% CI [%.4f, %.4f].\n', ...
    selectedStart, selectedEnd, stats.ci);
fprintf('First/last half block means: %.4f / %.4f; slope %.5f mN/m/ns.\n', ...
    firstHalf, lastHalf, slope);
fprintf('Last up-to-20 ns: %.4f +/- %.4f mN/m; log complete=%d, ERROR=%d.\n', ...
    lateStats.mean, lateStats.sem, logComplete, logError);
fprintf('Saved %d overlapping windows and all outputs to %s\n', height(windows), caseDir);
end

function gamma = calculateGamma(data)
gamma = 0.00506625 .* data.lz_A .* ...
    (data.pzz_atm - 0.5 .* (data.pxx_atm + data.pyy_atm));
end

function data = readEnergy(path)
if ~isfile(path), error('SurfaceTension:MissingFile', 'File not found: %s', path); end
lines = regexp(fileread(path), '\r\n|\n|\r', 'split');
required = {'time_fs','pxx_atm','pyy_atm','pzz_atm','lx_A','ly_A','lz_A'};
names = {}; rows = cell(numel(lines), 1); count = 0;
for k = 1:numel(lines)
    line = strtrim(lines{k});
    if isempty(line), continue; end
    if startsWith(line, '#')
        candidate = strsplit(strtrim(line(2:end)));
        if all(ismember(required, candidate))
            if ~isempty(names) && ~isequal(names, candidate)
                error('SurfaceTension:Header', 'Header changed within %s.', path);
            end
            names = candidate;
        end
        continue;
    end
    values = str2double(strsplit(line));
    if isempty(names) || numel(values) ~= numel(names) || any(~isfinite(values))
        error('SurfaceTension:Data', 'Invalid data/header at line %d of %s.', k, path);
    end
    count = count+1; rows{count} = values;
end
if count < 2, error('SurfaceTension:Data', 'Too few samples in %s.', path); end
data = array2table(vertcat(rows{1:count}), 'VariableNames', names);
dt = diff(data.time_fs);
if any(dt <= 0)
    error('SurfaceTension:Clock', 'Duplicate/decreasing times in %s; separate restarted runs.', path);
end
if any(abs(dt - median(dt)) > max(1e-6, 1e-9*median(dt)))
    error('SurfaceTension:Sampling', 'Irregular sampling or missing rows in %s.', path);
end
if any(data.lx_A <= 0 | data.ly_A <= 0 | data.lz_A <= 0)
    error('SurfaceTension:Box', 'Invalid box lengths in %s.', path);
end
end

function windows = windowMeans(time, gamma, beginNs, endNs, width, shift)
count = max(0, floor((endNs-beginNs-width)/shift + 1e-9)+1);
start_ns = beginNs + (0:count-1)'*shift;
end_ns = start_ns+width; center_ns = start_ns+width/2;
gamma_mN_m = NaN(count,1); samples = zeros(count,1);
for k = 1:count
    take = time >= start_ns(k) & time < end_ns(k);
    samples(k) = nnz(take);
    if samples(k) < 2
        error('SurfaceTension:Sampling', 'A window contains fewer than two samples.');
    end
    gamma_mN_m(k) = mean(gamma(take));
end
windows = table(start_ns,end_ns,center_ns,gamma_mN_m,samples);
end

function stats = blockStatistics(values)
stats.count = numel(values); stats.mean = mean(values);
stats.sem = NaN; stats.sd = NaN; stats.ci = [NaN NaN];
if stats.count >= 2
    stats.sd = std(values,0); stats.sem = stats.sd/sqrt(stats.count);
    df = stats.count-1;
    q = betaincinv(0.05,df/2,0.5);
    critical = sqrt(df*(1/q-1));
    stats.ci = stats.mean + [-1 1]*critical*stats.sem;
end
end

function [count,status] = wallContacts(data)
status = "absent"; count = NaN;
if all(ismember({'wall_lo_force','wall_hi_force'}, data.Properties.VariableNames))
    status = "available";
    count = nnz(data.wall_lo_force ~= 0 | data.wall_hi_force ~= 0);
end
end

function fig = newFigure()
fig = figure('Color','w','Visible','off','Units','pixels', ...
    'Position',[100 100 800 600],'PaperPositionMode','auto');
end

function styleAxes(ax)
ylabel(ax, '$\gamma\;(\mathrm{mN}\,\mathrm{m}^{-1})$', 'Interpreter','latex');
set(ax,'FontName','Times New Roman','FontSize',14,'LineWidth',1, ...
    'Box','on','TickDir','in');
end

function saveFigure(fig, stem)
savefig(fig,[stem '.fig']);
print(fig,[stem '.tif'],'-dtiff','-r300');
print(fig,[stem '.png'],'-dpng','-r150');
end

function plotTrend(windows,limits,boundary,stem)
fig = newFigure(); cleanup = onCleanup(@() close(fig));
ax = axes(fig);
plot(ax,windows.center_ns,windows.gamma_mN_m,'-', ...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5);
hold(ax,'on');
if ~isempty(boundary), xline(ax,boundary,'--k','HandleVisibility','off'); end
styleAxes(ax); xlim(ax,limits); xlabel(ax,'Time (ns)');
saveFigure(fig,stem);
end
