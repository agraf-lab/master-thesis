%% Developmental duration: EV versus cap-1 RNAi
% Requires the Statistics and Machine Learning Toolbox.
% Saves one SVG file. No CSV or MAT files are created.

clear;

%% 1. SETTINGS YOU MAY WANT TO CHANGE

data_file = 'Y:\shared\agraf\conditions_struct_cap1_RNAi_Anna_fresh.mat';
save_folder = 'Y:\shared\agraf';
svg_name = 'developmental_duration_cap1_RNAi_Anna_adj.svg';

% Events: 0 = Hatch, 1 = M1, 2 = M2, 3 = M3, 4 = M4.
% [3 4] means M3-M4.
% [1 3] means M1-M3.
% To show both intervals, use {[1 3], [3 4]}.
stage_intervals = {[3 4]};

outlier_threshold = 3;  % Remove extreme outliers using log-duration IQR.
data_fraction = 0.66;   % Lower 66% for data; upper 34% for brackets.

%% 2. CONDITIONS AND COLORS

% condition_map maps plot positions 1-6 to the 0-based condition indices
% in the .mat file. Correct order confirmed from Python labels:
%   condition_0: EV raga-1          -> plot position 1
%   condition_1: cap-1 raga-1       -> plot position 2
%   condition_2: EV col-10          -> plot position 3
%   condition_5: cap-1 col-10       -> plot position 4
%   condition_3: EV myo-2           -> plot position 5
%   condition_4: cap-1 myo-2        -> plot position 6
condition_map = [0 1 2 5 3 4];

background_names = {'raga1', 'col10', 'myo2'};
background_number = [1 1 2 2 3 3];

treatment_names = {'EV', 'cap1'};
treatment_number = [1 2 1 2 1 2];

stage_names = {'Hatch', 'M1', 'M2', 'M3', 'M4'};
duration_field = 'durations';

pair_labels = {
    {'raga-1::AID'; '500 uM IAA'}
    {'raga-1::AID'; 'col-10::TIR'; '500 uM IAA'}
    {'raga-1::AID'; 'myo-2::TIR'; '100 uM IAA'}
};

% Dark colors = EV; light colors = cap-1 RNAi.
colors = [
    0.08 0.20 0.45
    0.45 0.72 0.92
    0.35 0.10 0.55
    0.85 0.70 0.95
    0.05 0.35 0.15
    0.55 0.85 0.50
];

% Make the horizontal point jitter reproducible.
rng(1, 'twister');

validateattributes(data_fraction, {'numeric'}, ...
    {'scalar', '>', 0, '<', 1});

validateattributes(outlier_threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});

data = load(data_file);

% Verify that each condition index in condition_map points to the expected
% experimental group. This script accesses conditions by numeric index, not
% by name, so if the .mat file is re-exported with conditions in a different
% order the wrong animals would be analysed without any error or warning.
% The check below catches that silently wrong situation at load time.
expected_descriptions = {
    'EV, raga-1::AID, 500uM IAA'
    'cap-1 RNAi, raga-1::AID, 500uM IAA'
    'EV, raga-1::AID, col-10::TIR, 500uM IAA'
    'cap-1 RNAi, raga-1::AID, col-10::TIR, 500uM IAA'
    'EV, raga-1::AID, myo-2::TIR, 100uM IAA'
    'cap-1 RNAi, raga-1::AID, myo-2::TIR, 100uM IAA'
};
for group = 1:6
    fname = sprintf('condition_%d', condition_map(group));
    if isfield(data.(fname), 'description')
        actual = strtrim(data.(fname).description);
        expected = expected_descriptions{group};
        assert(strcmp(actual, expected), ...
            ['Condition mismatch at plot position %d.\n' ...
             '  Expected: %s\n' ...
             '  Found:    %s\n' ...
             'Update condition_map or expected_descriptions to match ' ...
             'the current .mat file.'], ...
            group, expected, actual);
    end
end

n_intervals = numel(stage_intervals);

assert(n_intervals > 0, 'Choose at least one stage interval.');

%% 3. CREATE THE FIGURE

% Figure dimensions: convert millimetres to inches.
figure_width = 150 * n_intervals / 25.4;
figure_height = 155 / 25.4;

fig = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1 1 figure_width figure_height]);

set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperSize', [figure_width figure_height], ...
    'PaperPosition', [0 0 figure_width figure_height], ...
    'PaperPositionMode', 'manual');

for panel = 1:n_intervals

    %% 4. CALCULATE AND CLEAN DURATIONS

    interval = stage_intervals{panel};

    validateattributes(interval, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', '>=', 0, '<=', 4});

    assert(interval(2) > interval(1), ...
        'The end stage must follow the start.');

    % Preserve the original stage-to-column convention.
    columns = (interval(1) + 1):interval(2);

    % One cell holds the cleaned durations for each condition.
    durations = cell(1, 6);

    for group = 1:6
        field_name = sprintf('condition_%d', condition_map(group));
        raw = double(data.(field_name).(duration_field));
        selected = raw(:, columns);

        % Exclude animals with any missing or infinite selected stage.
        complete = all(~isnan(selected) & ~isinf(selected), 2);

        % Add the selected stage durations for each animal.
        y = sum(selected(complete, :), 2);

        % Remove nonpositive durations and extreme outliers.
        durations{group} = clean_durations(y, outlier_threshold);
    end

    %% 5. DRAW THE SIX Box Plots

    ax = axes( ...
        'Parent', fig, ...
        'Units', 'normalized', ...
        'Position', [(panel - 1 + 0.15)/n_intervals, ...
                     0.26, 0.80/n_intervals, 0.68], ...
        'PositionConstraint', 'innerposition');

    hold(ax, 'on');

    set(ax, ...
        'YScale', 'log', ...
        'Box', 'on', ...
        'FontName', 'Arial', ...
        'FontSize', 9, ...
        'XTick', 1:6, ...
        'XTickLabel', ...
            {'EV', 'cap-1 RNAi', 'EV', 'cap-1 RNAi', 'EV', 'cap-1 RNAi'}, ...
        'XTickLabelRotation', 35, ...
        'TickLabelInterpreter', 'none');

    xlim(ax, [0.3 6.7]);

    for group = 1:6
        draw_box(ax, group, durations{group}, colors(group, :));
    end

    % Reserve space inside the plot for the statistical brackets.
    % Calculating this in log units keeps the y-axis accurate.
    all_y = vertcat(durations{:});

    if isempty(all_y)
        data_limits = [1 10];
    else
        data_limits = [min(all_y)/1.08, max(all_y)*1.08];
    end

    log_low = log(data_limits(1));
    log_height = log(data_limits(2)) - log_low;

    ylim(ax, exp([log_low, log_low + log_height/data_fraction]));

    title(ax, sprintf('%s-%s', ...
        stage_names{interval(1)+1}, stage_names{interval(2)+1}), ...
        'FontWeight', 'normal', ...
        'FontSize', 10);

    ylabel(ax, 'Developmental duration (h)');
    set_log_ticks(ax);

    % Add genotype and IAA labels below each pair.
    for pair = 1:3
        pair_center = 2*pair - 0.5;

        text(ax, (pair_center - 0.3)/6.4, -0.19, pair_labels{pair}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontName', 'Arial', ...
            'FontSize', 8, ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
    end

    %% 6. COMPARE EV AND cap-1 WITHIN EACH BACKGROUND

    % Five bracket heights, expressed as fractions of the plot height.
    bracket_levels = data_fraction + ...
        (1-data_fraction)*[0.22 0.34 0.46 0.60 0.80];

    fprintf('\nPanel %s-%s\n', stage_names{interval(1)+1}, stage_names{interval(2)+1});
    fprintf('%-8s  %5s  %5s  %12s  %12s\n', 'Pair', 'n EV', 'n cap1', 'p raw', 'p Holm');

    % Collect raw p-values for Holm correction across pairs.
    p_raw = nan(1, 3);
    for pair = 1:3
        ev   = 2*pair - 1;
        cap1 = 2*pair;
        if ~isempty(durations{ev}) && ~isempty(durations{cap1})
            p_raw(pair) = ranksum(durations{ev}, durations{cap1}, 'tail', 'both');
        end
    end

    % Holm correction across the three pairs.
    p_holm = holm_correct(p_raw);

    for pair = 1:3
        ev   = 2*pair - 1;
        cap1 = 2*pair;
        fprintf('%-8s  %5d  %5d  %12.4g  %12.4g\n', ...
            background_names{pair}, numel(durations{ev}), numel(durations{cap1}), ...
            p_raw(pair), p_holm(pair));

        draw_bracket(ax, ev, cap1, ...
            bracket_levels(pair), significance(p_holm(pair)));
    end

    %% 7. TEST WHETHER THE cap-1 EFFECT DIFFERS BETWEEN BACKGROUNDS

    % Model: logDuration ~ Background * CAP1.
    % P-values are unadjusted; animals are treated as independent.
    % Background includes IAA dose, which differs for the myo-2 pair.

    group_sizes = cellfun(@numel, durations);

    if all(group_sizes > 0) && sum(group_sizes) > 6

        % Build one table containing all retained animals.
        logDuration = [];
        Background = cell(0, 1);
        CAP1 = cell(0, 1);

        for group = 1:6
            n = numel(durations{group});

            logDuration = [logDuration; log(durations{group})]; %#ok<AGROW>

            Background = [Background; repmat( ...
                background_names(background_number(group)), n, 1)]; %#ok<AGROW>

            CAP1 = [CAP1; repmat( ...
                treatment_names(treatment_number(group)), n, 1)]; %#ok<AGROW>
        end

        % Set raga1 and EV as the reference categories.
        Background = reordercats( ...
            categorical(Background), background_names);

        CAP1 = reordercats(categorical(CAP1), treatment_names);

        model_data = table(logDuration, Background, CAP1);
        model = fitlm(model_data, 'logDuration ~ Background * CAP1');

        % Find the two interaction coefficients by their names.
        names = model.Coefficients.Properties.RowNames;
        H = zeros(2, numel(names));

        for comparison = 1:2
            term = ['Background_' background_names{comparison+1}];

            index = find( ...
                strcmp(names, [term ':CAP1_cap1']) | ...
                strcmp(names, ['CAP1_cap1:' term]));

            assert(numel(index) == 1, ...
                'Interaction coefficient not found.');

            H(comparison, index) = 1;
        end

        % Overall interaction test. These variables are not saved to a file.
        [p_interaction, F_interaction] = coefTest(model, H);

        % Compare each background's cap-1/EV fold change with raga1.
        for comparison = 1:2
            contrast = H(comparison, :);

            fc_ratio = exp(contrast * model.Coefficients.Estimate);
            p = coefTest(model, contrast);

            label = sprintf('FC_{ratio}=%.2f %s', ...
                fc_ratio, significance(p));

            draw_bracket(ax, 1.5, 2*comparison + 1.5, ...
                bracket_levels(comparison+3), label);
        end

    else
        text(ax, 0.5, 0.90, 'Interaction model unavailable', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'FontName', 'Arial', ...
            'FontSize', 8);
    end

    hold(ax, 'off');
end

%% 8. SAVE THE SVG — THE ONLY FILE WRITTEN BY THIS SCRIPT

if ~isfolder(save_folder)
    mkdir(save_folder);
end

drawnow;
print(fig, fullfile(save_folder, svg_name), '-dsvg');

%% HELPER FUNCTIONS
% Keep these at the bottom. You normally do not need to edit them.

function y = clean_durations(y, threshold)
    % Remove invalid durations.
    y = y(:);
    y = y(~isnan(y) & ~isinf(y) & y > 0);

    % Do not estimate outliers from fewer than four observations.
    if numel(y) < 4
        return;
    end

    % Apply the original IQR filter on log duration.
    log_y = log(y);
    q = prctile(log_y, [25 75]);
    spread = q(2) - q(1);

    if spread > 0
        keep = log_y >= q(1) - threshold*spread & ...
               log_y <= q(2) + threshold*spread;

        y = y(keep);
    end
end

function p_holm = holm_correct(p_raw)
    % Holm-Bonferroni correction across a vector of p-values.
    % NaN entries (missing tests) are left as NaN.
    n = numel(p_raw);
    p_holm = p_raw;
    valid = ~isnan(p_raw);
    idx = find(valid);
    [sorted_p, sort_order] = sort(p_raw(idx));
    m = numel(sorted_p);
    adjusted = zeros(1, m);
    running_max = 0;
    for k = 1:m
        adj = sorted_p(k) * (m - k + 1);
        running_max = max(running_max, adj);
        adjusted(k) = min(running_max, 1);
    end
    p_holm(idx(sort_order)) = adjusted;
end

function draw_box(ax, x, y, color)
    if isempty(y)
        return;
    end

    box_half = 0.25;   % half-width of the box (violins used 0.35)
    cap_half = 0.10;   % half-width of the whisker caps

    if numel(y) >= 2
        q = prctile(y, [25 50 75]);
        iqr_val = q(3) - q(1);

        % Whiskers: most extreme points within 1.5 * IQR of the box.
        w_low  = min(y(y >= q(1) - 1.5*iqr_val));
        w_high = max(y(y <= q(3) + 1.5*iqr_val));

        % Box (Q1 to Q3)
        patch(ax, ...
            [x-box_half x+box_half x+box_half x-box_half], ...
            [q(1) q(1) q(3) q(3)], ...
            color, ...
            'FaceAlpha', 0.25, ...
            'EdgeColor', color, ...
            'LineWidth', 0.8);

        % Whiskers and caps
        plot(ax, [x x], [q(3) w_high], '-', 'Color', color, 'LineWidth', 0.8);
        plot(ax, [x x], [w_low q(1)],  '-', 'Color', color, 'LineWidth', 0.8);
        plot(ax, [x-cap_half x+cap_half], [w_high w_high], '-', ...
            'Color', color, 'LineWidth', 0.8);
        plot(ax, [x-cap_half x+cap_half], [w_low w_low], '-', ...
            'Color', color, 'LineWidth', 0.8);
    end

    % Individual animals, jittered within the box width
    point_x = x + 0.8*box_half*(2*rand(size(y)) - 1);
    scatter(ax, point_x(:), y(:), 5, color, ...
        'filled', 'MarkerFaceAlpha', 0.5);

    % Median line drawn last so it sits on top of the points
    if numel(y) >= 2
        plot(ax, [x-box_half x+box_half], [q(2) q(2)], '-k', 'LineWidth', 1.2);
    end
end

function draw_bracket(ax, left, right, height, label)
    % Convert a fraction of plot height to a value on the log y-axis.
    limits = log(ylim(ax));
    to_y = @(fraction) exp(limits(1) + fraction*diff(limits));

    plot(ax, ...
        [left left right right], ...
        to_y([height-0.015 height height height-0.015]), ...
        '-k', 'LineWidth', 0.6);

    text(ax, mean([left right]), to_y(height+0.008), label, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', ...
        'FontSize', 8, ...
        'Interpreter', 'tex', ...
        'Clipping', 'on');
end

function set_log_ticks(ax)
    % Show intermediate values instead of only powers of ten.
    limits = ylim(ax);
    decades = 10.^(floor(log10(limits(1))):ceil(log10(limits(2))));

    ticks = sort(reshape([1; 2; 3; 5].*decades, 1, []));
    ticks = ticks(ticks >= limits(1) & ticks <= limits(2));

    if numel(ticks) < 4
        ticks = exp(linspace(log(limits(1)), log(limits(2)), 5));
        ticks([1 end]) = limits;
    end

    yticks(ax, ticks);
    yticklabels(ax, arrayfun(@(v) sprintf('%.3g', v), ...
        ticks, 'UniformOutput', false));
end

function label = significance(p)
    if isnan(p) || isinf(p)
        label = 'NA';
    elseif p < 0.0001
        label = '****';
    elseif p < 0.001
        label = '***';
    elseif p < 0.01
        label = '**';
    elseif p < 0.05
        label = '*';
    else
        label = 'ns';
    end
end

%"I found that the conditions in the .mat file were not in the order the script assumed — specifically, 
% col-10 cap-1 and the two myo-2 conditions were in a different sequence than expected. 
% This meant the script was plotting and testing the wrong animals for those groups, which is why the myo-2 significance differed 
% from the pipeline. I fixed this by adding an explicit condition mapping and a verification check that confirms each condition index points 
% to the right experimental group at load time, so the error cannot go undetected if the file is re-exported in a different order in the future.

% noticed the p-value for myo-2 differed between the MATLAB script and the ALIGN pipeline, traced it to a difference in n, and cross-checked 
% the condition labels in Python to find the ordering mismatch.