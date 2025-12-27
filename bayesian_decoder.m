%% BAYESIAN DECODER - Position Decoding from Hippocampal Neural Ensembles
%
% Description:
%   Decodes animal spatial position from population neural calcium activity
%   using naive Bayes decoder with Poisson model. Includes random sampling
%   of neurons across multiple iterations and shuffle control for 
%   statistical significance testing.
%
% Reference:
%   Yoo et al. - "Selective roles of astrocytes in the formation and 
%   stabilization of new hippocampal place fields"
%   Figure 1i-j, Figure 4: Decoding error analysis
%
% =========================================================================

function [results] = bayesian_decoder(data_path, varargin)
    %
    % Usage:
    %   results = bayesian_decoder(data_path)
    %   results = bayesian_decoder(data_path, 'num_iterations', 50, 'num_shuffles', 1000)
    %
    % Input:
    %   data_path [string] Path to folder containing 'information_data.mat'
    %       The MAT file must include the following variables:
    %         - SmoothMat_Total : 1×N cell array
    %             Each cell contains a 2D spatial firing-rate map
    %             for one neuron (e.g., 50×50 bins; units = Hz).
    %         - spkMat_Total    : 1×N cell array
    %             Each cell contains spike/event times and positions
    %             for one neuron in the format [t, x, y], where t is
    %             time in seconds and x,y are spatial coordinates (cm).
    %         - posTable_array  : T×3 double array
    %             Behavioral trajectory of the animal in the format
    %             [t, x, y], where t is time in seconds and x,y are
    %             positions (cm) sampled at a constant rate.
    %         - calmea          : N×2 double array
    %             Per-neuron summary metrics:
    %               calmea(:,1) = mean firing rate (Hz)
    %               calmea(:,2) = spatial information score
    %                             (bits/event, Skaggs measure).
    %
    % Optional name-value pairs:
    %   'num_iterations' [int] Number of random sampling iterations (default: 50)
    %   'num_shuffles' [int] Number of shuffle controls per iteration (default: 1000)
    %   'target_sample_size' [int] Neurons per iteration (default: 100)
    %   'min_required_cells' [int] Minimum eligible cells needed (default: 100)
    %   'time_bin_size' [int] Temporal bin in seconds (default: 2)
    %   'min_firing_rate' [double] Minimum firing rate threshold Hz (default: 0.01)
    %   'verbose' [logical] Print progress (default: true)
    %
    % Output:
    %   results [struct] Contains:
    %       .baseline_error_mean      - Mean decoding error (cm)
    %       .baseline_error_sem       - SEM of decoding error
    %       .shuffled_error_mean      - Mean shuffle error (cm)
    %       .shuffled_error_sem       - SEM of shuffle error
    %       .p_value                  - t-test significance
    %       .cohens_d                 - Effect size
    %       .decoding_error_by_minute - Time-resolved error [minutes × 1]
    %       .num_neurons_used         - Number of eligible neurons
    %       .num_iterations           - Iterations completed
    %       .analysis_time_min        - Total computation time (minutes)
    %
    % =========================================================================
    
    %% Parse inputs
    p = inputParser;
    addParameter(p, 'num_iterations', 50, @isnumeric);
    addParameter(p, 'num_shuffles', 1000, @isnumeric);
    addParameter(p, 'target_sample_size', 100, @isnumeric);
    addParameter(p, 'min_required_cells', 100, @isnumeric);
    addParameter(p, 'time_bin_size', 2, @isnumeric);
    addParameter(p, 'min_firing_rate', 0.01, @isnumeric);
    addParameter(p, 'verbose', true, @islogical);
    parse(p, varargin{:});
    
    params = p.Results;
    
    %% Load data
    if params.verbose
        fprintf('\n========================================\n');
        fprintf('Loading data from: %s\n', data_path);
    end
    
    load([data_path 'information_data.mat']);
    
    SmoothMat_Total = SmoothMat_Total;
    cell_numbers = size(SmoothMat_Total, 1);
    spkMat_Total = spkMat_Total;
    positionMat = posTable_array;
    calmea = calmea;
    
    if params.verbose
        fprintf('Loaded %d neurons\n', cell_numbers);
    end
    
    %% Time parameters
    test_start_timestamp = min(positionMat(:, 1));
    test_end_timestamp = max(positionMat(:, 1));
    test_duration = test_end_timestamp - test_start_timestamp;
    num_time_bins = floor(test_duration / params.time_bin_size);
    
    minute_bin_size = 60;
    bins_per_minute = minute_bin_size / params.time_bin_size;
    num_minute_bins = ceil(num_time_bins / bins_per_minute);
    
    %% Select eligible cells (firing rate threshold)
    min_firing_rate = params.min_firing_rate;
    eligible_cells = find(calmea(:, 1) >= min_firing_rate)';
    
    if length(eligible_cells) < params.min_required_cells
        fprintf('\nWarning: Insufficient eligible cells (%d < %d)\n', ...
                length(eligible_cells), params.min_required_cells);
        results = struct('status', 'insufficient_cells', ...
                        'available_cells', length(eligible_cells), ...
                        'required_cells', params.min_required_cells);
        return;
    end
    
    actual_sample_size = min(params.target_sample_size, length(eligible_cells));
    
    if params.verbose
        fprintf('Eligible neurons: %d (firing rate >= %.3f Hz)\n', ...
                length(eligible_cells), min_firing_rate);
        fprintf('\n========== BAYESIAN DECODING ========\n');
        fprintf('Iterations: %d | Shuffles: %d\n', params.num_iterations, params.num_shuffles);
        fprintf('Sample size: %d neurons per iteration\n', actual_sample_size);
        fprintf('Time bins: %d (%.1f min total)\n', num_time_bins, test_duration/60);
        fprintf('=====================================\n\n');
    end
    
    %% Initialize result storage
    all_error_distances_baseline = zeros(params.num_iterations, 1);
    all_error_distances_shuffled = zeros(params.num_iterations, params.num_shuffles);
    all_minute_errors_baseline = zeros(params.num_iterations, num_minute_bins);
    all_minute_errors_shuffled = zeros(params.num_iterations, num_minute_bins, params.num_shuffles);
    used_cells_per_iteration = cell(params.num_iterations, 1);
    
    total_start_time = tic;
    
    %% Main iteration loop
    for iter = 1:params.num_iterations
        iter_start_time = tic;
        
        if params.verbose && mod(iter, 10) == 1
            fprintf('Iteration %d/%d (%.1f%%)...\n', iter, params.num_iterations, (iter/params.num_iterations)*100);
        end
        
        % Random sampling of neurons
        rand_indices = eligible_cells(randperm(length(eligible_cells), actual_sample_size));
        used_cells_per_iteration{iter} = rand_indices;
        
        % Build population rate map
        pop_ratemap = [];
        for cell_idx = 1:actual_sample_size
            file_iter = rand_indices(cell_idx);
            if ~isempty(SmoothMat_Total{file_iter})
                curr_map = reshape(SmoothMat_Total{file_iter}, 1, []);
                pop_ratemap(cell_idx, :) = curr_map;
            end
        end
        
        if size(pop_ratemap, 1) < 10
            all_error_distances_baseline(iter) = NaN;
            all_error_distances_shuffled(iter, :) = NaN;
            all_minute_errors_baseline(iter, :) = NaN;
            all_minute_errors_shuffled(iter, :, :) = NaN;
            continue;
        end
        
        % Remove NaN locations (outside arena)
        nan_index = find(isnan(pop_ratemap(1, :)))';
        pop_ratemap(:, nan_index) = [];
        
        % Baseline decoding
        [baseline_error, baseline_minute_errors, decoded_pos, real_pos] = ...
            perform_baseline_decoding(pop_ratemap, spkMat_Total, rand_indices, ...
                                     positionMat, test_start_timestamp, ...
                                     params.time_bin_size, num_time_bins, ...
                                     bins_per_minute, num_minute_bins);
        
        all_error_distances_baseline(iter) = baseline_error;
        all_minute_errors_baseline(iter, :) = baseline_minute_errors;
        
        % Shuffle control
        for shuffle_iter = 1:params.num_shuffles
            shift_amount = randi([1, num_time_bins - 1]);
            shuffled_decoded_pos = circshift(decoded_pos, shift_amount);
            
            % Calculate shuffled errors
            shuffled_errors = zeros(num_time_bins, 1);
            for t = 1:num_time_bins
                if ~any(isnan(shuffled_decoded_pos(t, :))) && ~any(isnan(real_pos(t, :)))
                    curr_error = sqrt((shuffled_decoded_pos(t, 1) - real_pos(t, 1))^2 + ...
                                     (shuffled_decoded_pos(t, 2) - real_pos(t, 2))^2);
                    shuffled_errors(t) = curr_error;
                end
            end
            
            valid_shuffled = shuffled_errors(~isnan(shuffled_errors));
            if ~isempty(valid_shuffled)
                all_error_distances_shuffled(iter, shuffle_iter) = mean(valid_shuffled);
            end
            
            % 1-minute resolution for shuffled
            for minute_idx = 1:num_minute_bins
                start_bin = (minute_idx - 1) * bins_per_minute + 1;
                end_bin = min(minute_idx * bins_per_minute, num_time_bins);
                minute_errors = shuffled_errors(start_bin:end_bin);
                valid_minute = minute_errors(~isnan(minute_errors));
                if ~isempty(valid_minute)
                    all_minute_errors_shuffled(iter, minute_idx, shuffle_iter) = mean(valid_minute);
                end
            end
        end
        
        iter_time = toc(iter_start_time);
        if params.verbose && mod(iter, 10) == 0
            fprintf('  Completed in %.1f sec | Baseline error: %.3f cm\n', iter_time, baseline_error);
        end
    end
    
    total_time = toc(total_start_time);
    
    %% Statistical analysis
    valid_baseline = all_error_distances_baseline(~isnan(all_error_distances_baseline));
    valid_shuffled = all_error_distances_shuffled(~isnan(all_error_distances_shuffled));
    
    baseline_mean = mean(valid_baseline);
    baseline_std = std(valid_baseline);
    baseline_sem = baseline_std / sqrt(length(valid_baseline));
    
    shuffled_mean = mean(valid_shuffled);
    shuffled_std = std(valid_shuffled);
    shuffled_sem = shuffled_std / sqrt(length(valid_shuffled));
    
    [~, p_value] = ttest2(valid_baseline, valid_shuffled);
    
    pooled_std = sqrt(((length(valid_baseline) - 1) * baseline_std^2 + ...
                       (length(valid_shuffled) - 1) * shuffled_std^2) / ...
                       (length(valid_baseline) + length(valid_shuffled) - 2));
    if pooled_std > 0
        cohens_d = abs(baseline_mean - shuffled_mean) / pooled_std;
    else
        cohens_d = 0;
    end
    
    % 1-minute statistics
    mean_minute_baseline = nanmean(all_minute_errors_baseline, 1);
    mean_minute_shuffled = nanmean(nanmean(all_minute_errors_shuffled, 3), 1);
    
    if params.verbose
        fprintf('\n========== RESULTS ==========\n');
        fprintf('Baseline:  %.3f ± %.3f cm (n=%d)\n', baseline_mean, baseline_sem, length(valid_baseline));
        fprintf('Shuffled:  %.3f ± %.3f cm (n=%d)\n', shuffled_mean, shuffled_sem, length(valid_shuffled));
        fprintf('p-value: %.6f\n', p_value);
        fprintf('Cohen''s d: %.3f\n', cohens_d);
        fprintf('Time: %.1f min\n', total_time/60);
        fprintf('=============================\n\n');
    end
    
    %% Output structure
    results = struct();
    results.baseline_error_mean = baseline_mean;
    results.baseline_error_sem = baseline_sem;
    results.baseline_error_std = baseline_std;
    results.shuffled_error_mean = shuffled_mean;
    results.shuffled_error_sem = shuffled_sem;
    results.shuffled_error_std = shuffled_std;
    results.p_value = p_value;
    results.cohens_d = cohens_d;
    results.decoding_error_by_timebin = all_error_distances_baseline;
    results.decoding_error_by_minute = mean_minute_baseline;
    results.decoding_error_by_minute_shuffled = mean_minute_shuffled;
    results.num_neurons_used = length(eligible_cells);
    results.num_iterations = params.num_iterations;
    results.num_shuffles = params.num_shuffles;
    results.analysis_time_min = total_time / 60;
    results.parameters = params;
    results.all_baseline_errors = valid_baseline;
    results.all_shuffled_errors = valid_shuffled;
end

% =========================================================================
% Helper function
% =========================================================================

function [baseline_error, minute_errors, decoded_pos, real_pos] = ...
    perform_baseline_decoding(pop_ratemap, spkMat_Total, cell_indices, ...
                             positionMat, test_start_timestamp, ...
                             time_bin_size, num_time_bins, ...
                             bins_per_minute, num_minute_bins)
    
    num_cells = size(pop_ratemap, 1);
    decoded_pos = zeros(num_time_bins, 2);
    real_pos = zeros(num_time_bins, 2);
    timepoint_errors = zeros(num_time_bins, 1);
    
    % Decoding loop
    for time_bin = 1:num_time_bins
        current_start_time = test_start_timestamp + (time_bin - 1) * time_bin_size;
        current_end_time = current_start_time + time_bin_size;
        
        % Count spikes
        spike_counts = zeros(num_cells, 1);
        for i = 1:num_cells
            cell_idx = cell_indices(i);
            if ~isempty(spkMat_Total{cell_idx})
                spkMat = spkMat_Total{cell_idx};
                spkMat = [spkMat(:, 3), spkMat(:, 1:2)];  % Reorder: [timestamp, ...]
                spike_counts(i) = sum(spkMat(:, 1) >= current_start_time & ...
                                     spkMat(:, 1) < current_end_time);
            end
        end
        
        % Naive Bayes: log posterior = log likelihood + log prior
        log_likelihood = zeros(1, size(pop_ratemap, 2));
        for bin_idx = 1:size(pop_ratemap, 2)
            lambda = pop_ratemap(:, bin_idx);
            % Poisson: n*log(λ) - λ
            log_likelihood(bin_idx) = sum(spike_counts .* log(lambda + 1e-10) - lambda);
        end
        
        % Find decoded position
        [~, max_bin_1d] = max(log_likelihood);
        [decoded_y, decoded_x] = ind2sub([50, 50], max_bin_1d);
        decoded_pos(time_bin, :) = [decoded_x, decoded_y];
        
        % Real position
        pos_mask = positionMat(:, 1) >= current_start_time & ...
                   positionMat(:, 1) < current_end_time;
        if any(pos_mask)
            real_pos(time_bin, 1) = mean(positionMat(pos_mask, 2));
            real_pos(time_bin, 2) = mean(positionMat(pos_mask, 3));
            
            % Error
            error = sqrt((decoded_x - real_pos(time_bin, 1))^2 + ...
                        (decoded_y - real_pos(time_bin, 2))^2);
            timepoint_errors(time_bin) = error;
        else
            real_pos(time_bin, :) = [NaN, NaN];
            timepoint_errors(time_bin) = NaN;
        end
    end
    
    % Overall error
    valid_errors = timepoint_errors(~isnan(timepoint_errors));
    if ~isempty(valid_errors)
        baseline_error = mean(valid_errors);
    else
        baseline_error = NaN;
    end
    
    % 1-minute resolution
    minute_errors = zeros(num_minute_bins, 1);
    for minute_idx = 1:num_minute_bins
        start_bin = (minute_idx - 1) * bins_per_minute + 1;
        end_bin = min(minute_idx * bins_per_minute, num_time_bins);
        minute_data = timepoint_errors(start_bin:end_bin);
        valid_minute = minute_data(~isnan(minute_data));
        if ~isempty(valid_minute)
            minute_errors(minute_idx) = mean(valid_minute);
        else
            minute_errors(minute_idx) = NaN;
        end
    end

end
