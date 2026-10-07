function figures = plot_pdms_series()
%PLOT_PDMS_SERIES Overlay six SZ fits and their observed chain counts.
%
% Run from MATLAB:
%   addpath('simulations');
%   plot_pdms_series
%
% Reads the same model.conf chain_count rows consumed by the C++ generator.
% Saves three editable .fig figures and 300-dpi .tif images in
% simulations/figures. Solid lines use table_N*.csv SZ parameters; scatter
% points are the exact integer counts in the corresponding model.conf files.

base = fileparts(mfilename('fullpath'));
outputDir = fullfile(base, 'figures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

means = [16 32 64];
pdIs = {'1.05', '1.1', '1.15', '1.2', '1.25', '1.3'};
colors = [0.0000 0.4470 0.7410; ...
          0.8500 0.3250 0.0980; ...
          0.9290 0.6940 0.1250; ...
          0.4940 0.1840 0.5560; ...
          0.4660 0.6740 0.1880; ...
          0.3010 0.7450 0.9330];
figures = gobjects(numel(means), 1);

for i = 1:numel(means)
    meanN = means(i);
    fig = figure('Color', 'w', 'Name', sprintf('PDMS Mn %d SZ series', meanN));
    ax = axes('Parent', fig);
    hold(ax, 'on');
    labels = cell(1, numel(pdIs));
    parameters = read_sz_parameters(fullfile(base, ...
        sprintf('table_N%d.csv', meanN)));
    upper = 2 * meanN;
    if meanN == 16
        upper = 34;
    end
    caseCounts = cell(1, numel(pdIs));
    caseCurves = cell(1, numel(pdIs));
    yMax = 0;
    for j = 1:numel(pdIs)
        caseName = sprintf('N%d_PDI%s', meanN, pdIs{j});
        counts = read_chain_counts(fullfile(base, caseName, 'model.conf'));
        upperForCase = 2 * meanN;
        if meanN == 16 && j == numel(pdIs)
            upperForCase = 34;
        end
        discreteN = 4:upperForCase;
        curveN = linspace(4, upperForCase, 241);
        shape = parameters(j, 1);
        rate = parameters(j, 2);
        weights = discreteN.^(shape - 1) .* exp(-rate .* discreteN);
        curve = sum(counts(:, 2)) .* ...
            (curveN.^(shape - 1) .* exp(-rate .* curveN)) ./ sum(weights);
        caseCounts{j} = counts;
        caseCurves{j} = [curveN(:), curve(:)];
        yMax = max([yMax; counts(:, 2); curve(:)]);
        labels{j} = sprintf('PDI %s', pdIs{j});
    end
    if yMax <= 0
        error('plot_pdms_series:InvalidCounts', 'No positive chain counts found');
    end
    rawStep = yMax / 5;
    tickBase = 10 ^ floor(log10(rawStep));
    choices = [1 2 5 10] * tickBase;
    yStep = choices(find(choices >= rawStep, 1, 'first'));
    yLimit = ceil(yMax / yStep) * yStep;
    plot(ax, [meanN meanN], [0 yLimit], 'k--', ...
         'LineWidth', 1.4, 'HandleVisibility', 'off');
    for j = 1:numel(pdIs)
        plot(ax, [parameters(j, 3) parameters(j, 3)], [0 yLimit], '--', ...
             'Color', colors(j, :), 'LineWidth', 1.2, ...
             'HandleVisibility', 'off');
    end
    legendHandles = gobjects(1, numel(pdIs));
    for j = 1:numel(pdIs)
        curve = caseCurves{j};
        counts = caseCounts{j};
        plot(ax, curve(:, 1), curve(:, 2), ...
             'Color', colors(j, :), 'LineWidth', 2.0, ...
             'HandleVisibility', 'off');
        scatter(ax, counts(:, 1), counts(:, 2), 16, colors(j, :), ...
                'filled', 'MarkerEdgeColor', 'w', ...
                'HandleVisibility', 'off');
        legendHandles(j) = plot(ax, NaN, NaN, '-o', ...
            'Color', colors(j, :), 'LineWidth', 2.0, ...
            'MarkerFaceColor', colors(j, :), 'MarkerSize', 5);
    end
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [4 upper]);
    ylim(ax, [0 yLimit]);
    set(ax, 'YTick', 0:yStep:yLimit);
    xlabel(ax, '$n$', 'Interpreter', 'latex');
    ylabel(ax, '$M(n)$', 'Interpreter', 'latex');
    legend(ax, legendHandles, labels, ...
           'Location', 'northeast', 'Box', 'off', 'Color', 'none');
    set(ax, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.0);
    set(fig, 'Position', [100 + 40*i, 100 + 30*i, 1100, 825]);
    savefig(fig, fullfile(outputDir, sprintf('SZ_N%d.fig', meanN)));
    print(fig, fullfile(outputDir, sprintf('SZ_N%d.tif', meanN)), ...
          '-dtiff', '-r300');
    figures(i) = fig;
end
end

function parameters = read_sz_parameters(path)
fid = fopen(path, 'r');
if fid < 0
    error('plot_pdms_series:MissingTable', 'Cannot open %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fgetl(fid); % CSV header
columns = textscan(fid, '%s%f%f%f%f%f%f%s%f%f%f', ...
                   'Delimiter', ',', 'ReturnOnError', false);
if numel(columns{5}) ~= 6 || numel(columns{9}) ~= 6 || numel(columns{11}) ~= 6
    error('plot_pdms_series:InvalidTable', 'Expected six SZ rows in %s', path);
end
parameters = [columns{9}, columns{11}, columns{5}]; % k, rate, realized Mw
end

function counts = read_chain_counts(path)
if ~isfile(path)
    error('plot_pdms_series:MissingConfig', 'Missing config: %s', path);
end
raw = fileread(path);
tokens = regexp(raw, ...
    '(?m)^[ \t]*chain_count[ \t]*=[ \t]*(\d+)[ \t]+(\d+)[ \t]*(?:#.*)?$', ...
    'tokens');
if isempty(tokens)
    error('plot_pdms_series:MissingCounts', 'No chain_count rows in %s', path);
end
counts = zeros(numel(tokens), 2);
for row = 1:numel(tokens)
    counts(row, :) = [str2double(tokens{row}{1}), str2double(tokens{row}{2})];
end
end
