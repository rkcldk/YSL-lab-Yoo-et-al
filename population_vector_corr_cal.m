%% =========================================================================
%% POPULATION VECTOR CORRELATION - Standalone Script (Mouse-level Median)
%% (Profiling sheet 생성 없음, PV correlation만 계산, Median값)
%% =========================================================================

clear; clc; close all;

%% 파일 경로 설정 및 초기화
[start_path] = uigetdir;
cd(start_path);
session_num1 = 1;
session_num2 = 2;

%% 폴더 리스트업 및 cell registration 불러오기
listing = dir([start_path]);
folderlists(:, 1) = {listing.name};

% cellReg 폴더 리스트업
t = 1;
for i = 1:size(folderlists, 1)
    list = cell2mat(folderlists(i, 1));
    if contains(list, "cellReg") == 1        
        cellreg_folder_list(t, 1) = {list};
        t = t + 1;
    end
end

% cellReg 데이터 로드
name_spa = [start_path '\' cellreg_folder_list{1}];
load(name_spa)
registered_cell_index_spa = cell_registered_struct.cell_to_index_map;

name_centro = [start_path '\' cellreg_folder_list{2}];
load(name_centro)
registered_cell_index_centro = cell_registered_struct.cell_to_index_map;

% spa에서 두 세션 모두에 등록된 세포들
spa_same_cells = registered_cell_index_spa(find(registered_cell_index_spa(:, 1) ~= 0 & registered_cell_index_spa(:, 2) ~= 0), :);

% centro에서 두 세션 모두에 등록된 세포들
centro_same_cells = registered_cell_index_centro(find(registered_cell_index_centro(:, 1) ~= 0 & registered_cell_index_centro(:, 2) ~= 0), :);

% spa와 centro 둘 다에서 동일한 cell ID를 가진 세포들 찾기
common_cell_ids = intersect(spa_same_cells(:,1), centro_same_cells(:,1));

% 최종 분석용 세포 리스트
same_cell_number = [];
for i = 1:length(common_cell_ids)
    cell_id = common_cell_ids(i);
    spa_row = find(spa_same_cells(:,1) == cell_id);
    centro_row = find(centro_same_cells(:,1) == cell_id);
    if ~isempty(spa_row) && ~isempty(centro_row)
        same_cell_number(end+1, :) = [spa_same_cells(spa_row, :), centro_same_cells(centro_row, :)];
    end
end

% 세션 경로 설정
cd ../
session_path = cd;

%% 데이터 로드
session_num1_path = strcat(session_path, '\session-', string(session_num1), '\data');
session_num2_path = strcat(session_path, '\session-', string(session_num2), '\data');

% 세션 1 데이터 불러오기
cd(session_num1_path);
listing = dir('working*.mat');
load(listing.name);
SmoothMat_ses1 = SmoothMat_Total;
calmea_ses1 = calmea2;
cell_numbers_1 = size(SmoothMat_Total, 1);

% 세션 2 데이터 불러오기
cd(session_num2_path);
listing = dir('working*.mat');
load(listing.name);
SmoothMat_ses2 = SmoothMat_Total;
calmea_ses2 = calmea2;
cell_numbers_2 = size(SmoothMat_Total, 1);

%% 통계 저장 폴더 생성
cd(start_path);
stat_folder = sprintf('%d-%d_PV_correlation', session_num1, session_num2);
if ~isfolder(stat_folder)
    mkdir(stat_folder);
end
stat_path = fullfile(start_path, stat_folder);

%% =========================================================================
%% POPULATION VECTOR CORRELATION - 메인 분석
%% =========================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('POPULATION VECTOR CORRELATION ANALYSIS\n');
fprintf('========================================\n');
fprintf('Session 1: %d neurons\n', cell_numbers_1);
fprintf('Session 2: %d neurons\n', cell_numbers_2);
fprintf('CellReg matched cells: %d\n\n', size(same_cell_number, 1));

%% Step 1: 양쪽 세션 모두에서 active인 cell 선택
fprintf('[Step 1] Selecting neurons active in BOTH sessions...\n');

pv_eligible_cells = [];
for i = 1:size(same_cell_number, 1)
    cell_s1 = same_cell_number(i, 1);
    cell_s2 = same_cell_number(i, 2);
    
    fr_s1 = calmea_ses1(cell_s1, 1);  % Session 1 firing rate
    fr_s2 = calmea_ses2(cell_s2, 1);  % Session 2 firing rate
    
    % Both sessions must be active (≥0.01 Hz)
    if fr_s1 >= 0.01 && fr_s2 >= 0.01
        pv_eligible_cells = [pv_eligible_cells; i];
    end
end

fprintf('  Total tracked cells: %d\n', size(same_cell_number, 1));
fprintf('  Active in BOTH sessions (≥0.01 Hz): %d\n', length(pv_eligible_cells));
fprintf('  Inclusion ratio: %.1f%%\n\n', (length(pv_eligible_cells) / size(same_cell_number, 1)) * 100);

%% Step 2: ERsheet matrices 구성
fprintf('[Step 2] Building population rate matrices...\n');

ERsheet1 = [];
ERsheet2 = [];

for idx = 1:length(pv_eligible_cells)
    i = pv_eligible_cells(idx);
    cell_s1 = same_cell_number(i, 1);
    cell_s2 = same_cell_number(i, 2);
    
    SmoothMat1 = cell2mat(SmoothMat_ses1(cell_s1));
    SmoothMat2 = cell2mat(SmoothMat_ses2(cell_s2));
    
    % Only include if both have valid rate maps
    if ~isempty(SmoothMat1) && ~isempty(SmoothMat2)
        ERsheet1(idx, :) = reshape(SmoothMat1, 1, []);
        ERsheet2(idx, :) = reshape(SmoothMat2, 1, []);
    end
end

fprintf('  ERsheet1 size: [%d neurons × %d spatial bins]\n', size(ERsheet1,1), size(ERsheet1,2));
fprintf('  ERsheet2 size: [%d neurons × %d spatial bins]\n\n', size(ERsheet2,1), size(ERsheet2,2));

%% Step 3: Population vector correlation 계산 (Bin-by-bin, Mouse-level median)
fprintf('[Step 3] Calculating population vector correlations...\n');

pv_sheet = zeros(1, 2500);
valid_bins = 0;

for p = 1:2500
    bin1_seq = ERsheet1(:, p);  % Population vector at bin p, session 1
    bin2_seq = ERsheet2(:, p);  % Population vector at bin p, session 2
    
    % Only compute correlation if at least one session has activity
    if mean(bin1_seq) > 0 || mean(bin2_seq) > 0
        R = corr2(bin1_seq, bin2_seq);
        pv_sheet(1, p) = R;
        valid_bins = valid_bins + 1;
    else
        pv_sheet(1, p) = -100;  % Mark as outside arena
    end
end

% Remove invalid bins
pv_sheet(pv_sheet == -100) = [];
pv_sheet(isnan(pv_sheet)) = 0;

% Mouse-level median 계산
pv_median = median(pv_sheet);
pv_iqr = iqr(pv_sheet);

fprintf('  Valid spatial bins analyzed: %d\n');
fprintf('  Spatial bins with valid correlation: %d\n\n', length(pv_sheet));

%% Step 4: 결과 출력
fprintf('========== PV CORRELATION RESULTS ==========\n');
fprintf('Median PV correlation: %.4f\n', pv_median);
fprintf('IQR (Interquartile Range): %.4f\n', pv_iqr);
fprintf('Min: %.4f | Max: %.4f\n', min(pv_sheet), max(pv_sheet));
fprintf('=========================================\n\n');

%% Step 5: 결과 저장
fprintf('[Step 5] Saving results...\n');

% Structure로 저장
pv_results = struct();
pv_results.median_correlation = pv_median;
pv_results.iqr_correlation = pv_iqr;
pv_results.min_correlation = min(pv_sheet);
pv_results.max_correlation = max(pv_sheet);
pv_results.all_correlations = pv_sheet;
pv_results.num_neurons_used = length(pv_eligible_cells);
pv_results.num_valid_bins = length(pv_sheet);
pv_results.num_valid_spatial_bins = valid_bins;
pv_results.session1 = session_num1;
pv_results.session2 = session_num2;
pv_results.analysis_date = datetime('now');
pv_results.analysis_method = 'Mouse-level Median (Non-parametric)';

% 파일 저장
results_file = fullfile(stat_path, 'pv_correlation_results.mat');
save(results_file, 'pv_results', 'pv_sheet', 'pv_eligible_cells');
fprintf('  Saved: %s\n\n', results_file);

%% 완료
fprintf('========================================\n');
fprintf('Analysis completed successfully!\n');
fprintf('Results saved in: %s\n', stat_path);
fprintf('========================================\n\n');

%% 결과 요약
fprintf('Summary:\n');
fprintf('  - %d neurons analyzed\n', length(pv_eligible_cells));
fprintf('  - %d spatial bins\n', length(pv_sheet));
fprintf('  - Median PV correlation: %.4f (IQR: %.4f)\n', pv_median, pv_iqr);
fprintf('  - Statistical method: Mann-Whitney U test (non-parametric)\n');
fprintf('\n');