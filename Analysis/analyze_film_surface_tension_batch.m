function [summary, failures] = analyze_film_surface_tension_batch(caseNames, dataDir, outputDir, varargin)
%ANALYZE_FILM_SURFACE_TENSION_BATCH Analyze a list, or discover energy pairs.
% [summary, failures] = analyze_film_surface_tension_batch([], dataDir, outputDir)
% Discovers energy.*.film.dat and processes each case independently. A failed
% case is recorded rather than stopping the remaining cases. Other arguments
% are forwarded to analyze_film_surface_tension. Existing results for the
% same case/output directory are replaced; source files remain untouched.
if nargin < 1, caseNames = []; end
if nargin < 2 || isempty(dataDir), dataDir = pwd; end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if isempty(caseNames)
    files = dir(fullfile(dataDir,'energy.*.film.dat'));
    caseNames = strings(numel(files),1);
    for k = 1:numel(files)
        caseNames(k) = extractBefore(extractAfter(string(files(k).name),'energy.'),'.film.dat');
    end
else
    caseNames = string(caseNames);
    caseNames = caseNames(:);
end
caseNames = unique(caseNames,'stable');
if isempty(caseNames)
    error('SurfaceTension:NoCases','No production energy files/case names were supplied.');
end
if ~isfolder(outputDir), mkdir(outputDir); end
summary = table();
failures = table(strings(0,1),strings(0,1), ...
    'VariableNames',{'case_name','error_message'});
for k = 1:numel(caseNames)
    try
        result = analyze_film_surface_tension(caseNames(k),dataDir,outputDir,varargin{:});
        if isempty(summary)
            summary = result.summary;
        else
            summary = [summary; result.summary]; %#ok<AGROW>
        end
    catch err
        failures = [failures; {caseNames(k),string(err.message)}]; %#ok<AGROW>
        warning('SurfaceTension:CaseFailed','%s failed: %s',caseNames(k),err.message);
    end
end
writetable(failures,fullfile(outputDir,'surface_tension_failures.csv'));
if isempty(summary)
    summary = table(strings(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',{'case_name','gamma_mN_m','block_sem_mN_m'});
end
writetable(summary,fullfile(outputDir,'surface_tension_summary.csv'));
if ~isempty(summary)
    fig = figure('Color','w','Visible','off','Units','pixels', ...
        'Position',[100 100 800 600],'PaperPositionMode','auto');
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig);
    errorbar(ax,1:height(summary),summary.gamma_mN_m,summary.block_sem_mN_m, ...
        'o','LineWidth',1.5,'MarkerSize',7, ...
        'MarkerFaceColor',[0 0.4470 0.7410],'CapSize',12);
    set(ax,'XTick',1:height(summary),'XTickLabel',cellstr(summary.case_name), ...
        'TickLabelInterpreter','none','FontName','Times New Roman','FontSize',14, ...
        'LineWidth',1,'Box','on','TickDir','in');
    xtickangle(ax,45); xlim(ax,[0.5 height(summary)+0.5]);
    ylabel(ax,'$\gamma\;(\mathrm{mN}\,\mathrm{m}^{-1})$','Interpreter','latex');
    stem = fullfile(outputDir,'surface_tension_production_comparison');
    savefig(fig,[stem '.fig']);
    print(fig,[stem '.tif'],'-dtiff','-r300');
    print(fig,[stem '.png'],'-dpng','-r150');
end
fprintf('\nBatch: %d analyzed, %d failed. Inspect review flags before accepting values.\n', ...
    height(summary),height(failures));
end
