% BABYMACHINE_CCR_ANALYSIS
% Input: Table S4 beside this file or in an ancestor directory.
% TABLE_S4 may specify an explicit workbook path. Requires Optimization Toolbox.
% Usage: Babymachine_CCR_analysis('all') -- cell-number and RNA/DNA analyses.
%        Babymachine_CCR_analysis('omics') -- correction, fits and CCR only.
%        Babymachine_CCR_analysis('sync') -- cell-number fitting only.
% Output: results/baby_machine; run Babymachine_CCR_clustering_COG.R next.

function Babymachine_CCR_analysis(mode)
% Reproduce the published baby-machine analyses from Table S4.
% Usage: Babymachine_CCR_analysis('all')  % default; includes the full v4 grid search
%        Babymachine_CCR_analysis('omics')
%        Babymachine_CCR_analysis('sync')
% Run Babymachine_CCR_clustering_COG.R after the omics stage.
if nargin<1, mode='all'; end
if ~ismember(mode,{'all','omics','sync','plots','check'}), error('Unknown mode'); end
here=fileparts(mfilename('fullpath')); standardRoot=here;
inputFile=getenv('TABLE_S4');
if isempty(inputFile)
    searchDir=here; searched={};
    for depth=0:8
        searched{end+1}=searchDir; %#ok<AGROW>
        f=dir(fullfile(searchDir,'TableS4*.xlsx'));
        f=f(~[f.isdir]);
        if numel(f)>1
            error('Multiple Table S4 workbooks in %s. Set TABLE_S4 to the intended file. Candidates: %s', ...
                searchDir,strjoin({f.name},'; '));
        elseif numel(f)==1
            inputFile=fullfile(f(1).folder,f(1).name); break
        end
        parentDir=fileparts(searchDir);
        if strcmp(parentDir,searchDir), break; end
        searchDir=parentDir;
    end
    if isempty(inputFile)
        error('Table S4 not found. Set TABLE_S4 or place TableS4*.xlsx beside this script. Searched: %s',strjoin(searched,'; '));
    end
end
if ~isfile(inputFile), error('TABLE_S4 file does not exist: %s',inputFile); end
out=fullfile(standardRoot,'results','baby_machine');
if strcmp(mode,'check')
    fprintf('Input: %s\nOutput: %s\n',inputFile,out);
    return
end
work=fullfile(out,'work');
if ~isfolder(work), mkdir(work); end
fprintf('MATLAB %s\nInput: %s\n',version,inputFile);
if strcmp(mode,'plots')
    plot_dosage_diagnostics(fullfile(work,'TableS4_DNA_RNA_adapter.xlsx'),out, ...
        fullfile(out,'figures','copy_number_and_correlations'));
    return
end
if ismember(mode,{'all','omics'})
    D=readtable(inputFile,'Sheet','BabyM_DNA','VariableNamingRule','preserve');
    R=readtable(inputFile,'Sheet','BabyM_RNA','VariableNamingRule','preserve');
    % Original algorithms use exactly the first four annotations and 21 samples.
    times=0:11:220; sampleNames=cellstr(string(times)+" min");
    required=[{'ID','Gene name','Start','End'},sampleNames];
    if ~all(ismember(required,D.Properties.VariableNames)) || ...
            ~all(ismember(required,R.Properties.VariableNames))
        error('Unexpected TableS4 baby-machine column structure');
    end
    adapted=fullfile(work,'TableS4_DNA_RNA_adapter.xlsx');
    if isfile(adapted), delete(adapted); end
    % The original correction reads RNA relative position. Table S4 stores
    % this annotation in BabyM_DNA; join it by stable ID, never by row order.
    [positionPresent,positionIndex]=ismember(string(R.ID),string(D.ID));
    if any(R.Start(positionPresent)~=D.Start(positionIndex(positionPresent)) | ...
            R.End(positionPresent)~=D.End(positionIndex(positionPresent)))
        error('DNA and RNA reference coordinates disagree for matching IDs');
    end
    % Recover absent annotations from reference coordinates using the same
    % oriC and genome size as the supplied algorithms. Validate this
    % reconstruction against every supplied DNA position before using it.
    oriC=round((4579349+4579580)/2); genomeSize=4621599;
    dpos=((D.Start+D.End)/2-oriC)/(genomeSize/2);
    dpos(dpos < -1)=dpos(dpos < -1)+2;
    if max(abs(dpos-D.('Relative position')))>1e-10
        error('Published relative positions do not match the reference coordinate formula');
    end
    rpos=((R.Start+R.End)/2-oriC)/(genomeSize/2);
    rpos(rpos < -1)=rpos(rpos < -1)+2;
    rpos(positionPresent)=D.('Relative position')(positionIndex(positionPresent));
    R.('Relative position')=rpos;
    withPosition=[required,{'Relative position'}];
    writetable(D(:,withPosition),adapted,'Sheet','Gene_copy_num_M_Ala');
    writetable(R(:,withPosition),adapted,'Sheet','TPM_total_M_Ala');
    fit_gene_trajectories(adapted,fullfile(out,'BabyMachineDataAnalysis_output'));
    plot_dosage_diagnostics(adapted,out,fullfile(out,'figures','copy_number_and_correlations'));
end
if ismember(mode,{'all','sync'})
    P=readtable(inputFile,'Sheet','BabyM_Phen','VariableNamingRule','preserve');
    use=strcmp(string(P.Type),"Expt. Data");
    if sum(use)<2, error('No experimental cell-number data'); end
    % Fitting rows are a reference, never input observations.
    T=table(P.('Time (min)')(use),repmat("Replicate1",sum(use),1), ...
        P.('Relative cell number')(use), ...
        'VariableNames',{'Time (min)','Replicates','RelCellNO'});
    syncInput=fullfile(work,'Time_RelCellNO_from_TableS4.xlsx');
    writetable(T,syncInput,'Sheet','Sheet1');
    syncOut=fullfile(out,'synchrony_fit');
    if ~isfolder(syncOut), mkdir(syncOut); end
    % Checkpoint reuse must match both the observed data and the engine.
    signature=[sprintf('%.17g,',[T.('Time (min)'),T.RelCellNO]'), ...
        fileread([mfilename('fullpath') '.m'])];
    signatureFile=fullfile(syncOut,'checkpoint_input_and_engine_signature.txt');
    if isfile(signatureFile) && ~strcmp(fileread(signatureFile),signature)
        error('Existing synchrony results use different inputs/settings. Select a new output folder.');
    end
    fid=fopen(signatureFile,'w'); fprintf(fid,'%s',signature); fclose(fid);
    fit_cell_number(syncInput,syncOut);
    export_sync_tables(syncOut);
end
fprintf('Completed requested mode: %s\n',mode);
end

function s=baby_machine_config()
% Numerical defaults from the supplied paper scripts.
s.N0=1000; s.idt_grid=80:2:140;
s.cv_grid=0:.01:.30; s.asynchronous_grid=0:.01:.30;
s.substeps=20; s.plot_dt=.01; s.random_seed=20260915;
s.resume=true; s.minimum_free_memory_GB=1;
s.expression_max=10; s.expression_cycle_max=5;
s.division_boundary_min=118;
s.fc_threshold=1.5; s.damping_threshold=.75; s.ssr_ratio_threshold=.7;
% Original bounded sinusoidal model and peak/trough rules remain in
% modules/original_baby_analysis.m.
end
function fit_gene_trajectories(inputXlsx,outDir)
%% BabyMachineDataAnalysis
% Standalone analysis pipeline for synchronized baby-machine DNA-seq/RNA-seq data.
%
% Input:
%   total_genes_M_Ala_240627.xlsx
%   Required sheets:
%     1. Gene_copy_num_M_Ala
%     2. TPM_total_M_Ala
%
% Main outputs:
%   1. Gene_copy_num_Simulation
%   2. Low-expression filter QC and sensitivity summaries
%   3. TPM_Simulation for expression-eligible genes
%   4. 50-kb binned DNA copy-number and TPM matrices
%   5. DNA copy-number-corrected TPM matrix
%   6. corrected_TPM_Simulation
%   7. CCR gene tables
%
% The correction uses regional DNA copy-number trajectories averaged in 50-kb bins,
% rather than single-gene DNA trajectories, to reduce gene-level copy-number noise.



dir_mkdir(outDir);
outXlsx=fullfile(outDir,'BabyMachineDataAnalysis_results.xlsx');
if isfile(outXlsx), delete(outXlsx); end
settings=baby_machine_config();
timeVals = 0:11:220;
binSize = 50000;
oriC = round((4579349 + 4579580) / 2);
genomeSize = 4621599;

% Low-expression filter applied to raw TPM before copy-number correction
% and periodicity analysis. The complete time course must reach TPM 10,
% and each of the two division cycles must independently reach TPM 5.
% With the measured division time near 118 min, sampled times 0-110 min
% define cycle 1 and 121-220 min define cycle 2.
divisionBoundaryMin = settings.division_boundary_min;
expressionMaxThreshold = settings.expression_max;
expressionPerCycleThreshold = settings.expression_cycle_max;

% Threshold pairs [whole-time-course maximum, per-cycle maximum] used only
% for sensitivity reporting. The middle row is the primary analysis.
expressionSensitivityThresholds = [
     5,  2;
    10,  5;
    20, 10
];
expressionSensitivityLabels = ["lenient"; "primary"; "stringent"];

% CCR screening criteria
% fcThreshold: simulated peak-to-valley fold change in at least one
% synchronous cycle.
% dampingThreshold: maximum allowed absolute value of p6, which controls
% cycle-to-cycle amplitude damping or amplification in the sinusoidal model.
% sinLinearSsrRatioThreshold: maximum allowed SSR_sin/SSR_lin ratio; values
% below this threshold indicate that the sinusoidal model fits better than a
% linear trend.
fcThreshold = settings.fc_threshold;
dampingThreshold = settings.damping_threshold;
sinLinearSsrRatioThreshold = settings.ssr_ratio_threshold;

fprintf('BabyMachineDataAnalysis started.\n');
fprintf('Input workbook:\n%s\n', inputXlsx);
fprintf('Output folder:\n%s\n\n', outDir);

%% Read input sheets
fprintf('Reading input sheets...\n');
dfDPM = readtable(inputXlsx, 'Sheet', 'Gene_copy_num_M_Ala', ...
    'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
dfTPM = readtable(inputXlsx, 'Sheet', 'TPM_total_M_Ala', ...
    'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
dfTPM = remove_extra_annotation_columns(dfTPM);

% Sort RNA table by ID to make downstream outputs deterministic.
dfTPM = sortrows(dfTPM, {'ID'});

dpmRaw = table2array(dfDPM(:, 5:25));
tpmRawAll = table2array(dfTPM(:, 5:25));

%% Filter low-expression genes before RNA correction and periodicity fitting
fprintf('Filtering low-expression genes from raw TPM trajectories...\n');
cycle1Mask = timeVals < divisionBoundaryMin;
cycle2Mask = timeVals >= divisionBoundaryMin;
if ~any(cycle1Mask) || ~any(cycle2Mask)
    error('The division boundary does not split timeVals into two non-empty cycles.');
end

maxTpmAll = max(tpmRawAll, [], 2, 'omitnan');
maxTpmCycle1 = max(tpmRawAll(:, cycle1Mask), [], 2, 'omitnan');
maxTpmCycle2 = max(tpmRawAll(:, cycle2Mask), [], 2, 'omitnan');
validTpmCount = sum(isfinite(tpmRawAll), 2);

passWholeCourse = isfinite(maxTpmAll) & maxTpmAll >= expressionMaxThreshold;
passCycle1 = isfinite(maxTpmCycle1) & maxTpmCycle1 >= expressionPerCycleThreshold;
passCycle2 = isfinite(maxTpmCycle2) & maxTpmCycle2 >= expressionPerCycleThreshold;
keepExpression = passWholeCourse & passCycle1 & passCycle2;

expressionFilterStatus = [dfTPM(:, 1:4), table( ...
    maxTpmAll, maxTpmCycle1, maxTpmCycle2, validTpmCount, ...
    passWholeCourse, passCycle1, passCycle2, keepExpression, ...
    'VariableNames', {'max_TPM_all', 'max_TPM_cycle1', 'max_TPM_cycle2', ...
    'n_finite_timepoints', 'pass_max_all', 'pass_cycle1', 'pass_cycle2', ...
    'keep_for_CCR_analysis'})];

expressionFilterSummary = table( ...
    height(dfTPM), sum(keepExpression), sum(~keepExpression), ...
    100 * sum(keepExpression) / height(dfTPM), ...
    expressionMaxThreshold, expressionPerCycleThreshold, divisionBoundaryMin, ...
    min(timeVals(cycle1Mask)), max(timeVals(cycle1Mask)), ...
    min(timeVals(cycle2Mask)), max(timeVals(cycle2Mask)), ...
    'VariableNames', {'n_input_genes', 'n_retained_genes', 'n_removed_genes', ...
    'pct_retained', 'max_TPM_all_threshold', 'max_TPM_each_cycle_threshold', ...
    'division_boundary_min', 'cycle1_start_min', 'cycle1_end_min', ...
    'cycle2_start_min', 'cycle2_end_min'});

nSensitivity = size(expressionSensitivityThresholds, 1);
sensitivityRetained = zeros(nSensitivity, 1);
for i = 1:nSensitivity
    maxThreshold = expressionSensitivityThresholds(i, 1);
    cycleThreshold = expressionSensitivityThresholds(i, 2);
    sensitivityRetained(i) = sum( ...
        isfinite(maxTpmAll) & maxTpmAll >= maxThreshold & ...
        isfinite(maxTpmCycle1) & maxTpmCycle1 >= cycleThreshold & ...
        isfinite(maxTpmCycle2) & maxTpmCycle2 >= cycleThreshold);
end
expressionFilterSensitivity = table( ...
    expressionSensitivityLabels, expressionSensitivityThresholds(:, 1), ...
    expressionSensitivityThresholds(:, 2), sensitivityRetained, ...
    100 * sensitivityRetained / height(dfTPM), ...
    'VariableNames', {'threshold_set', 'max_TPM_all_threshold', ...
    'max_TPM_each_cycle_threshold', 'n_retained_genes', 'pct_retained'});

writetable(expressionFilterStatus, outXlsx, 'Sheet', 'expression_filter_status');
writetable(expressionFilterSummary, outXlsx, 'Sheet', 'expression_filter_summary');
writetable(expressionFilterSensitivity, outXlsx, 'Sheet', 'expression_filter_sensitivity');
writetable(expressionFilterStatus, fullfile(outDir, 'expression_filter_status.csv'));
writetable(expressionFilterSummary, fullfile(outDir, 'expression_filter_summary.csv'));
writetable(expressionFilterSensitivity, fullfile(outDir, 'expression_filter_sensitivity.csv'));
writetable(dfTPM(~keepExpression, :), fullfile(outDir, 'low_expression_genes.csv'));

fprintf(['  Retained %d / %d genes (%.1f%%): max TPM >= %g overall and ', ...
    'max TPM >= %g in each cycle.\n'], ...
    sum(keepExpression), height(dfTPM), ...
    100 * sum(keepExpression) / height(dfTPM), ...
    expressionMaxThreshold, expressionPerCycleThreshold);

dfTPM = dfTPM(keepExpression, :);
tpmRaw = table2array(dfTPM(:, 5:25));
writetable(dfTPM, fullfile(outDir, 'expression_eligible_TPM_total_M_Ala.csv'));

%% Generate sinusoidal simulation tables
fprintf('Fitting DNA copy-number trajectories...\n');
geneCopyNumSimulation = fit_sinusoidal_table(dfDPM(:, 1:4), dpmRaw, timeVals);
writetable(geneCopyNumSimulation, outXlsx, 'Sheet', 'Gene_copy_num_Simulation');
writetable(geneCopyNumSimulation, fullfile(outDir, 'Gene_copy_num_Simulation.csv'));

fprintf('Fitting uncorrected TPM trajectories...\n');
tpmSimulation = fit_sinusoidal_table(dfTPM(:, 1:4), tpmRaw, timeVals);
writetable(tpmSimulation, outXlsx, 'Sheet', 'TPM_Simulation');
writetable(tpmSimulation, fullfile(outDir, 'TPM_Simulation.csv'));

%% Generate 50-kb binned DNA copy-number and TPM matrices
fprintf('Generating 50-kb binned DNA copy-number and TPM matrices...\n');
binnedDNA = make_binned_matrix(dfDPM, timeVals, binSize, oriC, genomeSize);
binnedTPM = make_binned_matrix(dfTPM, timeVals, binSize, oriC, genomeSize);
log2RelBinnedDNA = make_log2_relative_matrix(binnedDNA, timeVals);
log2RelBinnedTPM = make_log2_relative_matrix(binnedTPM, timeVals);

writetable(binnedDNA, outXlsx, 'Sheet', 'Binned_DNA_copy_number');
writetable(binnedTPM, outXlsx, 'Sheet', 'Binned_TPM');
writetable(log2RelBinnedDNA, outXlsx, 'Sheet', 'log2_rel_binned_DNA');
writetable(log2RelBinnedTPM, outXlsx, 'Sheet', 'log2_rel_binned_TPM');
writetable(binnedDNA, fullfile(outDir, 'Binned_DNA_copy_number.csv'));
writetable(binnedTPM, fullfile(outDir, 'Binned_TPM.csv'));
writetable(log2RelBinnedDNA, fullfile(outDir, 'log2_relative_binned_DNA_matrix.csv'));
writetable(log2RelBinnedTPM, fullfile(outDir, 'log2_relative_binned_TPM_matrix.csv'));

%% Correct TPM using regional DNA copy-number trajectories
fprintf('Correcting TPM by regional DNA copy-number trajectories...\n');
[correctedTPM, tpmToDnaBinMatch] = correct_tpm_by_binned_dna(dfTPM, binnedDNA, timeVals);
writetable(correctedTPM, outXlsx, 'Sheet', 'corrected_TPM_total_M_Ala');
writetable(tpmToDnaBinMatch, outXlsx, 'Sheet', 'TPM_to_DNA_bin_match');
writetable(correctedTPM, fullfile(outDir, 'corrected_TPM_total_M_Ala.csv'));
writetable(tpmToDnaBinMatch, fullfile(outDir, 'TPM_to_DNA_bin_match.csv'));

fprintf('Fitting corrected TPM trajectories...\n');
correctedTpmSimulation = fit_sinusoidal_table(correctedTPM(:, 1:4), table2array(correctedTPM(:, 5:25)), timeVals);
writetable(correctedTpmSimulation, outXlsx, 'Sheet', 'corrected_TPM_Simulation');
writetable(correctedTpmSimulation, fullfile(outDir, 'corrected_TPM_Simulation.csv'));

%% Identify CCR genes from corrected TPM trajectories
fprintf('Screening CCR genes from corrected TPM trajectories...\n');
ccrMask = identify_ccr_genes(correctedTpmSimulation, fcThreshold, dampingThreshold, sinLinearSsrRatioThreshold);
ccrSummary = make_ccr_summary(correctedTPM, correctedTpmSimulation, ccrMask);
nonCcrSummary = make_ccr_summary(correctedTPM, correctedTpmSimulation, ~ccrMask);

writetable(ccrSummary, outXlsx, 'Sheet', 'CCR_genes');
writetable(nonCcrSummary, outXlsx, 'Sheet', 'non_CCR_genes');
writetable(ccrSummary, fullfile(outDir, 'CCR_genes.csv'));
writetable(nonCcrSummary, fullfile(outDir, 'non_CCR_genes.csv'));

fprintf('\nAnalysis completed.\n');
fprintf('CCR genes identified: %d / %d\n', sum(ccrMask), numel(ccrMask));
fprintf('Low-expression genes excluded before correction and fitting: %d\n', sum(~keepExpression));
fprintf('Main output workbook:\n%s\n', outXlsx);

end

function dir_mkdir(p)
if ~isfolder(p), mkdir(p); end
end

%% Local functions
function simTable = fit_sinusoidal_table(geneInfo, rawTimeData, timeVals)
    % Sinusoidal model:
    % y = p1 * (1 - p6)^(t / p2) * sin(2*pi*t/p2 + 2*pi*p3) + p4*t/p2 + p5
    %
    % Parameter meanings:
    % p1: initial oscillation amplitude.
    % p2: oscillation period in minutes.
    % p3: phase shift, expressed as a fraction of one oscillation period.
    % p4: linear baseline slope across the time series.
    % p5: baseline offset.
    % p6: damping parameter; positive values reduce oscillation amplitude
    %     between cycles, whereas negative values increase amplitude.
    % p7: reserved parameter retained for compatibility with the original
    %     output table; it is fixed at 0 and is not used in the fitted model.
    p0 = [0.1, 110, -0.5, 0, 0, 0.2, 0];
    pars = [0, 5; 100, 130; -1, 0; -0.5, 0.5; -2, 2; -1, 1; 0, 0];
    nGenes = size(rawTimeData, 1);
    simCurve = nan(nGenes, 39);
    X = timeVals;

    for n = 1:nGenes
        if mod(n, 500) == 0
            fprintf('  fitted %d / %d\n', n, nGenes);
        end
        y = rawTimeData(n, :);
        if any(~isfinite(y)) || any(y <= 0)
            continue;
        end
        Y = log2(y ./ mean(y));
        try
            evalc('[Sim, ~, R2, P, PV_t] = Sinusoidal_simulation_embedded(X, Y, p0, pars, 11);');
            simY = 2 .^ Sim .* mean(y);
            fc = nan(1, 2);
            fc(1) = value_at_time(simY, X, PV_t(1)) ./ value_at_time(simY, X, PV_t(2));
            fc(2) = value_at_time(simY, X, PV_t(3)) ./ value_at_time(simY, X, PV_t(4));
            simCurve(n, :) = [simY, R2(1:3), P, PV_t, fc, R2(2) ./ R2(5)];
        catch ME
            warning('Sinusoidal fitting failed for row %d: %s', n, ME.message);
        end
    end

    simTable = [geneInfo, array2table(simCurve)];
    simTable = splitvars(simTable);
    simTable.Properties.VariableNames = [geneInfo.Properties.VariableNames, ...
        string(timeVals), {'R^2_sin', 'SSR_sin', 'SST_sin'}, ...
        {'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7'}, ...
        {'peak_time1', 'valley_time1', 'peak_time2', 'valley_time2', 'replicate_time'}, ...
        {'FC1_Sim', 'FC2_Sim', 'SSR_sin/SSR_lin'}];
end

function v = value_at_time(y, X, t)
    [~, idx] = min(abs(X - t));
    v = y(idx);
end

function [SimCurve, LinCurve, R2, p_sin, T] = Sinusoidal_simulation_embedded(X, Y, p0, parameters, dt)
    fun = @(par, xdata) par(1) .* (1 - par(6)) .^ (xdata / par(2)) .* ...
        sin(2 * pi / par(2) .* xdata + 2 * pi * par(3)) + ...
        par(4) * xdata / par(2) + par(5);

    p_sin = lsqcurvefit(fun, p0, X, Y, parameters(:, 1), parameters(:, 2));
    SimCurve = fun(p_sin, X);

    SSR_sin = sum((SimCurve - Y) .^ 2);
    SST_sin = sum((Y - mean(Y)) .^ 2);
    R2_sin = 1 - SSR_sin / SST_sin;

    p_lin = polyfit(X, Y, 1);
    LinCurve = polyval(p_lin, X);

    SSR_lin = sum((LinCurve - Y) .^ 2);
    SST_lin = sum((Y - mean(Y)) .^ 2);
    R2_lin = 1 - SSR_lin / SST_lin;

    R2 = [R2_sin, SSR_sin, SST_sin, R2_lin, SSR_lin, SST_lin];

    t = min(X):dt:max(X);
    tlist = t / p_sin(2) + p_sin(3) + 1;
    while sum(tlist > 1) > 0
        tlist(tlist > 1) = tlist(tlist > 1) - 1;
    end

    peak = t(abs(tlist - 0.25) < dt / p_sin(2) / 2);
    if isempty(peak)
        [~, idx] = max(fun(p_sin, X));
        peak = X(idx);
    end
    if numel(peak) == 1
        if fun(p_sin, X(1)) >= fun(p_sin, X(end))
            peak = [X(1), peak];
        else
            peak = [peak, X(end)];
        end
    end

    valley = t(abs(tlist - 0.75) < dt / p_sin(2) / 2);
    if isempty(valley)
        [~, idx] = min(fun(p_sin, X));
        valley = X(idx);
    end
    if numel(valley) == 1
        if fun(p_sin, X(1)) <= fun(p_sin, X(end))
            valley = [X(1), valley];
        else
            valley = [valley, X(end)];
        end
    end

    if peak(1) < 0
        peak(1) = 0;
        peak(2) = peak(2) + dt;
    end
    if peak(2) > X(end)
        peak(2) = X(end);
        peak(1) = peak(1) - dt;
    end
    if valley(1) < 0
        valley(1) = 0;
        valley(2) = valley(2) + dt;
    end
    if valley(2) > X(end)
        valley(2) = X(end);
        valley(1) = valley(1) - dt;
    end

    replicate_t = (0.25 - 1 / 4 - p_sin(3)) * p_sin(2);
    while replicate_t < 0
        replicate_t = replicate_t + p_sin(2);
    end
    while replicate_t >= p_sin(2)
        replicate_t = replicate_t - p_sin(2);
    end

    T = [peak(1), valley(1), peak(2), valley(2), replicate_t];
end

function binned = make_binned_matrix(df, timeVals, binSize, oriC, genomeSize)
    binCenter = [sort(oriC - (1:round(genomeSize / 2 / binSize)) * binSize, 'ascend'), ...
        oriC, oriC + (1:round(genomeSize / 2 / binSize)) * binSize - 1];
    binCenter(binCenter > genomeSize) = binCenter(binCenter > genomeSize) - genomeSize;
    binCenter(binCenter <= 0) = binCenter(binCenter <= 0) + genomeSize;
    binLeft = binCenter - binSize / 2;
    binRight = binCenter + binSize / 2 - 1;
    binLeft(binLeft <= 0) = binLeft(binLeft <= 0) + genomeSize;
    binRight(binRight > genomeSize) = binRight(binRight > genomeSize) - genomeSize;
    relPos = (binCenter - oriC) / (genomeSize / 2);
    relPos(relPos < -1) = relPos(relPos < -1) + 2;

    genePos = round((table2array(df(:, {'Start'})) + table2array(df(:, {'End'}))) / 2);
    values = table2array(df(:, 5:25));
    out = nan(numel(binCenter), numel(timeVals));

    for i = 1:numel(binCenter)
        if binLeft(i) <= binRight(i)
            inBin = genePos >= binLeft(i) & genePos <= binRight(i);
        else
            inBin = genePos >= binLeft(i) | genePos <= binRight(i);
        end
        if any(inBin)
            out(i, :) = mean(values(inBin, :), 1, 'omitnan');
        end
    end

    binned = array2table([binCenter(:), binLeft(:), binRight(:), relPos(:), out], ...
        'VariableNames', [{'Bin_center', 'Bin_start', 'Bin_end', 'Relative_position'}, string(timeVals)]);
end

function log2Rel = make_log2_relative_matrix(binned, timeVals)
    values = table2array(binned(:, 5:25));
    rel = log2(values ./ mean(values, 2, 'omitnan'));
    log2Rel = [binned(:, 1:4), array2table(rel, 'VariableNames', string(timeVals))];
end

function [corrTPM, matchTable] = correct_tpm_by_binned_dna(dfTPM, binnedDNA, timeVals)
    tpmData = table2array(dfTPM(:, 5:25));
    dnaData = table2array(binnedDNA(:, 5:25));
    dnaRel = dnaData ./ mean(dnaData, 2, 'omitnan');
    tpmPos = table2array(dfTPM(:, {'Relative position'}));
    binPos = table2array(binnedDNA(:, {'Relative_position'}));

    corrected = nan(size(tpmData));
    corrCoef = nan(height(dfTPM), 1);
    matchedBin = nan(height(dfTPM), 1);
    matchedBinPos = nan(height(dfTPM), 1);

    for i = 1:height(dfTPM)
        [~, idx] = min(abs(binPos - tpmPos(i)));
        matchedBin(i) = idx;
        matchedBinPos(i) = binPos(idx);
        corrected(i, :) = tpmData(i, :) ./ dnaRel(idx, :);
        corrCoef(i) = corr(tpmData(i, :)', dnaRel(idx, :)', 'Rows', 'complete');
    end

    corrTPM = [dfTPM(:, 1:4), array2table(corrected, 'VariableNames', string(timeVals)), dfTPM(:, 26:end)];
    if any(strcmp(corrTPM.Properties.VariableNames, 'corr_DNA'))
        corrTPM(:, {'corr_DNA'}) = array2table(corrCoef);
    end
    matchTable = table(table2array(dfTPM(:, {'ID'})), table2array(dfTPM(:, {'Gene name'})), tpmPos, matchedBin, matchedBinPos, ...
        'VariableNames', {'ID', 'Gene_name', 'TPM_relative_position', 'Matched_DNA_bin_index', 'Matched_DNA_bin_relative_position'});
end

function ccrMask = identify_ccr_genes(simTable, fcThreshold, dampingThreshold, sinLinearSsrRatioThreshold)
    fc1 = table2array(simTable(:, {'FC1_Sim'}));
    fc2 = table2array(simTable(:, {'FC2_Sim'}));
    p6 = table2array(simTable(:, {'p6'}));
    ratio = table2array(simTable(:, {'SSR_sin/SSR_lin'}));
    ccrMask = ((fc1 >= fcThreshold) | (fc2 >= fcThreshold)) & ...
        abs(p6) <= dampingThreshold & ratio < sinLinearSsrRatioThreshold;
end

function summaryTable = make_ccr_summary(exprTable, simTable, mask)
    summaryTable = [exprTable(mask, :), ...
        array2table(table2array(simTable(mask, 5:25)), 'VariableNames', strcat('S', string(0:11:220))), ...
        simTable(mask, 26:end)];
end

function tbl = remove_extra_annotation_columns(tbl)
    names = tbl.Properties.VariableNames;
    drop = startsWith(names, 'Unnamed') | startsWith(names, 'Var') | strcmp(names, 'CRP_regulated');
    tbl(:, drop) = [];
end
function fit_cell_number(inputFile,outputDir)
%% BabyMachine_Fitting_FromRawData_v4_checkpoint
% Independent, checkpointed fitting for Time_RelCellNO.xlsx.
% v4 corrects the legacy 20-25 min IDT search to a range that includes
% the approximately 110 min doubling time supplied for this experiment.


runStart = tic;

settings=baby_machine_config();
%% User settings
inputFileName = inputFile;
sheetName = 'Sheet1';
outputFolderName = 'BabyMachine_Fitting_v4_checkpoint_outputs';

N0 = settings.N0;                       % simulated founder cells; not measured cell count
idtlist = settings.idt_grid;              % candidate mean inter-division time (min)
CVlist = settings.cv_grid;            % candidate division-time coefficient of variation
USlist = settings.asynchronous_grid;            % candidate initially unsynchronized fraction
substepsPerInterval = settings.substeps;        % preserve original simulation resolution
plotDt = settings.plot_dt;                   % display-grid spacing only (min)
randomSeed = settings.random_seed;           % reproducible Monte Carlo seed
resumeFromCheckpoint = settings.resume;     % continue after an interrupted MATLAB session
minimumFreeMemoryGB = settings.minimum_free_memory_GB;       % save checkpoint and stop before low-memory state

if ~isfolder(outputDir), mkdir(outputDir); end

%% Read and validate raw data
T = readtable(inputFile, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
requiredVariables = {'Time (min)', 'Replicates', 'RelCellNO'};
if ~all(ismember(requiredVariables, T.Properties.VariableNames))
    error('Workbook must contain exactly these required columns: %s', ...
        strjoin(requiredVariables, ', '));
end
timeAll = T.('Time (min)');
replicateAll = string(T.('Replicates'));
relCellAll = T.('RelCellNO');
if ~isnumeric(timeAll) || ~isnumeric(relCellAll) || ...
        any(~isfinite(timeAll)) || any(~isfinite(relCellAll)) || ...
        any(timeAll < 0) || any(relCellAll <= 0) || any(ismissing(replicateAll))
    error('Required columns contain invalid time, replicate, or RelCellNO values.');
end
replicateNames = unique(replicateAll, 'stable');
nReplicates = numel(replicateNames);
if nReplicates == 0
    error('No replicate data found.');
end

summaryReplicate = strings(nReplicates, 1);
summaryN = nan(nReplicates, 1);
summaryTEnd = nan(nReplicates, 1);
summarySeed = nan(nReplicates, 1);
summaryIDT = nan(nReplicates, 1);
summaryCV = nan(nReplicates, 1);
summaryUS = nan(nReplicates, 1);
summaryRD = nan(nReplicates, 1);

fprintf('Input: %s\n', inputFile);
fprintf('IDT search: %s min\n', mat2str(idtlist));
fprintf('CV-US combinations per IDT: %d\n', numel(CVlist) * numel(USlist));

%% Independent fitting of each replicate
for iRep = 1:nReplicates
    replicateName = replicateNames(iRep);
    safeName = matlab.lang.makeValidName(char(replicateName));
    workspaceFile = fullfile(outputDir, [safeName '_fit_workspace.mat']);
    checkpointFile = fullfile(outputDir, [safeName '_checkpoint.mat']);
    progressFile = fullfile(outputDir, [safeName '_progress.csv']);

    useRows = replicateAll == replicateName;
    exp_t = timeAll(useRows);
    exp_N = relCellAll(useRows);
    [exp_t, sortIndex] = sort(exp_t(:));
    exp_N = exp_N(sortIndex);
    if numel(exp_t) < 2 || any(diff(exp_t) <= 0)
        error('%s must contain at least two strictly increasing time points.', replicateName);
    end
    t_end = max(exp_t);
    replicateSeed = randomSeed + iRep - 1;

    % Skip an already complete replicate only when all model settings match.
    if resumeFromCheckpoint && isfile(workspaceFile)
        completed = load(workspaceFile, 'idtlist', 'CVlist', 'USlist', ...
            'bestIDT', 'bestCV', 'bestUS', 'bestRD', ...
            'N0','substepsPerInterval','replicateSeed','exp_t','exp_N');
        if isequal(completed.idtlist, idtlist) && isequal(completed.CVlist, CVlist) && ...
                isequal(completed.USlist, USlist) && ...
                isequal(completed.N0,N0) && isequal(completed.substepsPerInterval,substepsPerInterval) && ...
                isequal(completed.replicateSeed,replicateSeed) && ...
                isequal(completed.exp_t,exp_t) && isequal(completed.exp_N,exp_N)
            fprintf('\n%s already complete; retaining existing result.\n', replicateName);
            summaryReplicate(iRep) = replicateName;
            summaryN(iRep) = numel(exp_t);
            summaryTEnd(iRep) = t_end;
            summarySeed(iRep) = replicateSeed;
            summaryIDT(iRep) = completed.bestIDT;
            summaryCV(iRep) = completed.bestCV;
            summaryUS(iRep) = completed.bestUS;
            summaryRD(iRep) = completed.bestRD;
            continue
        end
    end

    rd = nan(numel(idtlist), numel(CVlist), numel(USlist));
    nextIDTIndex = 1;
    if resumeFromCheckpoint && isfile(checkpointFile)
        partial = load(checkpointFile, 'rd', 'nextIDTIndex', 'idtlist', 'CVlist', 'USlist');
        if isequal(partial.idtlist, idtlist) && isequal(partial.CVlist, CVlist) && ...
                isequal(partial.USlist, USlist) && isequal(size(partial.rd), size(rd))
            rd = partial.rd;
            nextIDTIndex = partial.nextIDTIndex;
            fprintf('\nResuming %s at IDT %d of %d.\n', ...
                replicateName, nextIDTIndex, numel(idtlist));
        else
            warning('Ignoring incompatible checkpoint for %s.', replicateName);
        end
    end

    fprintf('\nFitting %s: n=%d, t_end=%.6g min, seed=%d\n', ...
        replicateName, numel(exp_t), t_end, replicateSeed);
    replicateStart = tic;

    for iIDT = nextIDTIndex:numel(idtlist)
        idt = idtlist(iIDT);
        fprintf('  IDT %.6g min (%d/%d)\n', idt, iIDT, numel(idtlist));

        for iCV = 1:numel(CVlist)
            checkAvailableMemory(minimumFreeMemoryGB, checkpointFile, rd, iIDT, ...
                idtlist, CVlist, USlist, replicateName);
            fprintf('    CV %.2f (%d/%d); parameter sets %d-%d of %d\n', ...
                CVlist(iCV), iCV, numel(CVlist), ...
                (iCV - 1) * numel(USlist) + 1, iCV * numel(USlist), ...
                numel(CVlist) * numel(USlist));

            for iUS = 1:numel(USlist)
                % Common reproducible random stream for every candidate.
                rng(replicateSeed, 'twister');
                rd(iIDT, iCV, iUS) = simulateRelativeDeviation( ...
                    exp_t, exp_N, N0, idt, CVlist(iCV), USlist(iUS), ...
                    substepsPerInterval, false);
            end
        end

        nextIDTIndex = iIDT + 1;
        save(checkpointFile, 'rd', 'nextIDTIndex', 'idtlist', 'CVlist', 'USlist', ...
            'replicateName', 'replicateSeed', 'exp_t', 'exp_N', 't_end', '-v7');
        progressTable = table(replicateName, iIDT, numel(idtlist), idt, ...
            100 * iIDT / numel(idtlist), toc(replicateStart) / 60, ...
            'VariableNames', {'Replicate', 'CompletedIDTIndex', 'TotalIDTCount', ...
            'CompletedIDT_min', 'PercentComplete', 'ElapsedMinutes'});
        writetable(progressTable, progressFile);
        fprintf('  Checkpoint saved: %.1f%% of %s complete.\n', ...
            100 * iIDT / numel(idtlist), replicateName);
    end

    %% Derive best candidate and rerun it once for plotting
    [bestRD, linearIndex] = min(rd(:));
    [bestIDTIndex, bestCVIndex, bestUSIndex] = ind2sub(size(rd), linearIndex);
    bestIDT = idtlist(bestIDTIndex);
    bestCV = CVlist(bestCVIndex);
    bestUS = USlist(bestUSIndex);
    rng(replicateSeed, 'twister');
    [~, bestTDense, bestNDense, bestNAtExperimentalTimes] = simulateRelativeDeviation( ...
        exp_t, exp_N, N0, bestIDT, bestCV, bestUS, substepsPerInterval, true);

    %% Save numerical outputs
    [IDTgrid, CVgrid, USgrid] = ndgrid(idtlist, CVlist, USlist);
    gridResults = table(IDTgrid(:), CVgrid(:), USgrid(:), rd(:), ...
        'VariableNames', {'IDT_min', 'CV', 'UnsynchronizedFraction', 'RelativeDeviation'});
    writetable(gridResults, fullfile(outputDir, [safeName '_parameter_grid.csv']));
    save(workspaceFile, 'replicateName', 'exp_t', 'exp_N', 't_end', 'N0', ...
        'idtlist', 'CVlist', 'USlist', 'substepsPerInterval', 'plotDt', ...
        'replicateSeed', 'rd', 'bestRD', 'bestIDT', 'bestCV', 'bestUS', ...
        'bestTDense', 'bestNDense', 'bestNAtExperimentalTimes', '-v7');

    %% Save figures without keeping graphical windows open
    minRDByIDT = squeeze(min(min(rd, [], 3), [], 2));
    rdAtBestIDT = squeeze(rd(bestIDTIndex, :, :));
    h1 = figure('Color', 'w', 'Visible', 'off');
    plot(idtlist, minRDByIDT, 'o-', 'LineWidth', 2); grid on;
    xlabel('Inter-division time (min)'); ylabel('Best fitted relative deviation');
    title(sprintf('%s: minimum deviation versus IDT', replicateName), 'Interpreter', 'none');
    saveFigurePair(h1, outputDir, [safeName '_min_rd']); close(h1);

    h2 = figure('Color', 'w', 'Visible', 'off');
    imagesc(CVlist, USlist, rdAtBestIDT');
    set(gca, 'XDir', 'normal', 'YDir', 'normal');
    xlim([0.10, 0.20]);
    ylim([0, 0.20]);
    xticks(0.10:0.02:0.20);
    yticks(0:0.05:0.20);
    hold on;
    scatter(bestCV, bestUS, 45, 'rx', 'LineWidth', 1.5); colorbar;
    xlabel('CV'); ylabel('Unsynchronized fraction');
    title(sprintf('%s: IDT = %.1f min', replicateName, bestIDT), 'Interpreter', 'none');
    saveFigurePair(h2, outputDir, [safeName '_heatmap']); close(h2);

    expNorm = exp_N / exp_N(1);
    denseNorm = bestNDense / bestNAtExperimentalTimes(1);
    plotTime = unique([(bestTDense(1):plotDt:bestTDense(end))'; bestTDense(end)]);
    plotNorm = interp1(bestTDense, denseNorm, plotTime, 'previous');
    h3 = figure('Color', 'w', 'Visible', 'off');
    semilogy(exp_t, expNorm, 'o', 'LineWidth', 2, 'DisplayName', 'Experimental data'); hold on;
    semilogy(plotTime, plotNorm, '-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Fitted simulation; RD = %.4f', bestRD));
    xlabel('Time (min)'); ylabel('Normalized cell number'); grid on; legend('Location', 'best');
    title(sprintf('%s: IDT = %.1f, CV = %.2f, unsynchronized = %.2f', ...
        replicateName, bestIDT, bestCV, bestUS), 'Interpreter', 'none');
    saveFigurePair(h3, outputDir, [safeName '_best_fitting']); close(h3);

    if isfile(checkpointFile)
        delete(checkpointFile);
    end
    summaryReplicate(iRep) = replicateName;
    summaryN(iRep) = numel(exp_t);
    summaryTEnd(iRep) = t_end;
    summarySeed(iRep) = replicateSeed;
    summaryIDT(iRep) = bestIDT;
    summaryCV(iRep) = bestCV;
    summaryUS(iRep) = bestUS;
    summaryRD(iRep) = bestRD;
end

%% Final summary
summaryTable = table(summaryReplicate, summaryN, summaryTEnd, summarySeed, ...
    summaryIDT, summaryCV, summaryUS, summaryRD, ...
    'VariableNames', {'Replicate', 'ObservationCount', 't_end_min', 'RandomSeed', ...
    'BestIDT_min', 'BestCV', 'BestUnsynchronizedFraction', 'BestRelativeDeviation'});
writetable(summaryTable, fullfile(outputDir, 'independent_fit_summary.csv'));
matlabVersion = string(version);
settingsTable = table(string(inputFileName), string(sheetName), N0, ...
    string(mat2str(idtlist)), string(mat2str(CVlist)), string(mat2str(USlist)), ...
    substepsPerInterval, plotDt, randomSeed, minimumFreeMemoryGB, matlabVersion, ...
    'VariableNames', {'InputFile', 'Sheet', 'N0', 'IDTGrid', 'CVGrid', 'USGrid', ...
    'SubstepsPerInterval', 'PlotDt_min', 'BaseRandomSeed', 'MinimumFreeMemoryGB', ...
    'MATLABVersion'});
writetable(settingsTable, fullfile(outputDir, 'analysis_settings.csv'));
disp(summaryTable);
fprintf('\nCompleted. Outputs saved to: %s\n', outputDir);
toc(runStart);

end

%% Local functions
function varargout = simulateRelativeDeviation(exp_t, exp_N, N0, idt, CV, uscell, substepsPerInterval, returnTrajectory)
    anchorTimes = unique([0; exp_t(:)], 'stable');
    tDense = anchorTimes(1);
    for i = 2:numel(anchorTimes)
        intervalGrid = linspace(anchorTimes(i - 1), anchorTimes(i), substepsPerInterval + 1)';
        tDense = [tDense; intervalGrid(2:end)]; %#ok<AGROW>
    end
    CV_idt = idt * CV;
    NDense = nan(size(tDense)); NDense(1) = N0;
    divisionTime = idt + randn(N0, 1) * CV_idt;
    nUnsynchronized = ceil(N0 * uscell);
    if nUnsynchronized > 0
        divisionTime(1:nUnsynchronized) = idt * rand(nUnsynchronized, 1);
    end
    divisionTime(divisionTime < 0) = 0;
    timeToDivision = divisionTime;
    currentN = N0;
    for i = 2:numel(tDense)
        timeToDivision = timeToDivision - (tDense(i) - tDense(i - 1));
        dividing = timeToDivision < 0;
        nDividing = sum(dividing);
        if nDividing > 0
            timeToDivision(dividing) = idt + randn(nDividing, 1) * CV_idt;
            timeToDivision = [timeToDivision; idt + randn(nDividing, 1) * CV_idt]; %#ok<AGROW>
            currentN = currentN + nDividing;
        end
        NDense(i) = currentN;
    end
    NAtExp = interp1(tDense, NDense, exp_t, 'previous');
    simulatedNorm = NAtExp / NAtExp(1);
    experimentalNorm = exp_N / exp_N(1);
    rd = mean(abs(simulatedNorm(:) - experimentalNorm(:)) ./ experimentalNorm(:));
    if returnTrajectory
        varargout = {rd, tDense, NDense, NAtExp};
    else
        varargout = {rd};
    end
end

function checkAvailableMemory(minimumFreeMemoryGB, checkpointFile, rd, nextIDTIndex, idtlist, CVlist, USlist, replicateName)
    if ~ispc
        return
    end
    mem = memory;
    if mem.MemAvailableAllArrays < minimumFreeMemoryGB * 1024^3
        save(checkpointFile, 'rd', 'nextIDTIndex', 'idtlist', 'CVlist', 'USlist', 'replicateName', '-v7');
        error('Stopped safely: available MATLAB memory is below %.2f GB. Restart this script to resume.', minimumFreeMemoryGB);
    end
end

function saveFigurePair(h, outputDir, baseName)
    savefig(h, fullfile(outputDir, [baseName '.fig']));
    saveas(h, fullfile(outputDir, [baseName '.png']));
end
function plot_dosage_diagnostics(inputXlsx,dataDir,outDir)
expressionFilterCsv=fullfile(dataDir,'BabyMachineDataAnalysis_output','expression_filter_status.csv');
if ~isfolder(outDir), mkdir(outDir); end
timeVals = 0:11:220;
timeNames = string(timeVals);
oriC = round((4579349 + 4579580) / 2);
genomeSize = 4621599;
binSize = 50000;
axisSizeCm = 5.04;
axisLineWidthPt = 0.5;
fontName = 'Arial';
fontSizePt = 8;
relTicks = [-1, -0.5, 0, 0.5, 1];
corrTicks = [-1, 0, 1];
relSignalTicks = [-0.5, 0, 0.5];

fprintf('Reading input workbook...\n%s\n', inputXlsx);
dfDPM = readtable(inputXlsx, 'Sheet', 'Gene_copy_num_M_Ala', ...
    'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
dfTPM = readtable(inputXlsx, 'Sheet', 'TPM_total_M_Ala', ...
    'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
dfTPM = plot_remove_extra_annotation_columns(dfTPM);

% Apply exactly the same upstream low-expression filter used for the CCR
% analysis. Keeping the filter in BabyMachineDataAnalysis.m as the single
% source of truth prevents the plotting and screening pipelines from drifting.
if exist(expressionFilterCsv, 'file') ~= 2
    error(['Expression-filter status file not found:\n%s\n', ...
        'Run BabyMachineDataAnalysis.m before this plotting script.'], ...
        expressionFilterCsv);
end
expressionFilter = readtable(expressionFilterCsv, ...
    'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
requiredFilterColumns = {'ID', 'keep_for_CCR_analysis'};
if ~all(ismember(requiredFilterColumns, expressionFilter.Properties.VariableNames))
    error('Expression-filter status file is missing ID or keep_for_CCR_analysis.');
end

filterFlagRaw = expressionFilter.keep_for_CCR_analysis;
if islogical(filterFlagRaw)
    filterFlag = filterFlagRaw;
elseif isnumeric(filterFlagRaw)
    filterFlag = isfinite(filterFlagRaw) & filterFlagRaw ~= 0;
else
    filterFlag = ismember(lower(strtrim(string(filterFlagRaw))), ["true", "1", "yes"]);
end
eligibleIds = string(expressionFilter.ID(filterFlag));
inputIds = string(dfTPM.ID);
keepExpression = ismember(inputIds, eligibleIds);
missingEligibleIds = setdiff(eligibleIds, inputIds);
if ~isempty(missingEligibleIds)
    error('%d expression-eligible IDs were not found in TPM_total_M_Ala.', ...
        numel(missingEligibleIds));
end
fprintf('Applying upstream RNA expression filter: retained %d / %d genes.\n', ...
    sum(keepExpression), height(dfTPM));
dfTPM = dfTPM(keepExpression, :);

dpmRaw = table2array(dfDPM(:, 5:25));
tpmRaw = table2array(dfTPM(:, 5:25));
dpmRel = dpmRaw ./ mean(dpmRaw, 2, 'omitnan');
tpmRel = tpmRaw ./ mean(tpmRaw, 2, 'omitnan');

%% Sample correlations and interval summaries: sample-sample periodic reproducibility
fprintf('Computing Sample correlations and interval summaries...\n');
dnaSampleCorr = corr(dpmRel, dpmRel, 'Rows', 'pairwise');
rnaSampleCorr = corr(tpmRel, tpmRel, 'Rows', 'pairwise');
[dnaInterval, dnaMean, dnaSd] = plot_summarize_corr_by_interval(dnaSampleCorr, timeVals);
[rnaInterval, rnaMean, rnaSd] = plot_summarize_corr_by_interval(rnaSampleCorr, timeVals);

fig27 = plot_make_manual_figure();

axes('Parent', fig27, 'Units', 'centimeters', 'Position', plot_panel_position(1));
imagesc(timeVals, timeVals, dnaSampleCorr);
axis image; set(gca, 'YDir', 'normal');
title('DNA-seq'); xlabel('Time (min)'); ylabel('Time (min)');
colormap(gca, plot_redblue(256)); caxis([-1, 1]); cb = colorbar; cb.Ticks = corrTicks;
plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig27, 'Units', 'centimeters', 'Position', plot_panel_position(2));
errorbar(dnaInterval, dnaMean, dnaSd, '-o', 'Color', [0.2, 0.35, 0.85], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1);
xlabel('Sampling interval (min)'); ylabel('Pearson correlation');
title('DNA-seq'); xlim([0, 220]); ylim([-1, 1]); grid on;
plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig27, 'Units', 'centimeters', 'Position', plot_panel_position(3));
imagesc(timeVals, timeVals, rnaSampleCorr);
axis image; set(gca, 'YDir', 'normal');
title('RNA-seq'); xlabel('Time (min)'); ylabel('Time (min)');
colormap(gca, plot_redblue(256)); caxis([-1, 1]); cb = colorbar; cb.Ticks = corrTicks;
plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig27, 'Units', 'centimeters', 'Position', plot_panel_position(4));
errorbar(rnaInterval, rnaMean, rnaSd, '-o', 'Color', [0.2, 0.35, 0.85], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1);
xlabel('Sampling interval (min)'); ylabel('Pearson correlation');
title('RNA-seq'); xlim([0, 220]); ylim([-1, 1]); grid on;
plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

plot_save_figure(fig27, outDir, 'Figs_periodic_reproducibility');

%% Binned DNA and RNA patterns: binned copy-number and transcriptional patterns
fprintf('Computing Binned DNA and RNA patterns...\n');
binnedDPM = plot_make_binned_matrix(dfDPM, timeVals, binSize, oriC, genomeSize);
binnedTPM = plot_make_binned_matrix(dfTPM, timeVals, binSize, oriC, genomeSize);
binnedDpmVals = table2array(binnedDPM(:, 5:25));
binnedTpmVals = table2array(binnedTPM(:, 5:25));
binnedDpmLogRel = log2(binnedDpmVals ./ mean(binnedDpmVals, 2, 'omitnan'));
binnedTpmLogRel = log2(binnedTpmVals ./ mean(binnedTpmVals, 2, 'omitnan'));
binPos = table2array(binnedDPM(:, {'Relative_position'}));
[binPosSorted, binOrder] = sort(binPos, 'ascend');
binnedDpmLogRelSorted = binnedDpmLogRel(binOrder, :);
binnedTpmLogRelSorted = binnedTpmLogRel(binOrder, :);
dnaBinCorr = corr(binnedDpmLogRelSorted', binnedDpmLogRelSorted', 'Rows', 'pairwise');
rnaDnaBinCorr = corr(binnedTpmLogRelSorted', binnedDpmLogRelSorted', 'Rows', 'pairwise');

fig28 = plot_make_manual_figure();

axes('Parent', fig28, 'Units', 'centimeters', 'Position', plot_panel_position(1));
imagesc(timeVals, binPosSorted, binnedDpmLogRelSorted);
set(gca, 'YDir', 'normal'); xlabel('Time (min)'); ylabel('Relative genome position');
title('DNA-seq'); colormap(gca, plot_redblue(256)); caxis([-0.5, 0.5]); cb = colorbar; cb.Ticks = relSignalTicks;
yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig28, 'Units', 'centimeters', 'Position', plot_panel_position(2));
imagesc(binPosSorted, binPosSorted, dnaBinCorr);
axis image; set(gca, 'YDir', 'normal');
xlabel('Relative genome position of DNA bins'); ylabel('Relative genome position of DNA bins');
title('DNA-DNA correlation'); colormap(gca, plot_redblue(256)); caxis([-1, 1]); cb = colorbar; cb.Ticks = corrTicks;
xticks(relTicks); yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig28, 'Units', 'centimeters', 'Position', plot_panel_position(3));
imagesc(timeVals, binPosSorted, binnedTpmLogRelSorted);
set(gca, 'YDir', 'normal'); xlabel('Time (min)'); ylabel('Relative genome position');
title('RNA-seq'); colormap(gca, plot_redblue(256)); caxis([-0.5, 0.5]); cb = colorbar; cb.Ticks = relSignalTicks;
yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig28, 'Units', 'centimeters', 'Position', plot_panel_position(4));
imagesc(binPosSorted, binPosSorted, rnaDnaBinCorr);
axis image; set(gca, 'YDir', 'normal');
xlabel('Relative genome position of DNA bins'); ylabel('Relative genome position of RNA bins');
title('RNA-DNA correlation'); colormap(gca, plot_redblue(256)); caxis([-1, 1]); cb = colorbar; cb.Ticks = corrTicks;
xticks(relTicks); yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

plot_save_figure(fig28, outDir, 'Figs_replication_transcription_patterns');

%% DNA dosage correction: copy-number correction at gene level
fprintf('Computing DNA dosage correction...\n');
[corrTPM, matchTable] = plot_correct_tpm_by_binned_dpm(dfTPM, binnedDPM, timeVals);
corrTpmRaw = table2array(corrTPM(:, 5:25));
corrTpmRel = corrTpmRaw ./ mean(corrTpmRaw, 2, 'omitnan');
corrTpmLogRel = log2(corrTpmRel);
tpmLogRel = log2(tpmRel);
genePos = table2array(dfTPM(:, {'Relative position'}));
[genePosSorted, geneOrder] = sort(genePos, 'ascend');
tpmLogRelSorted = tpmLogRel(geneOrder, :);
corrTpmLogRelSorted = corrTpmLogRel(geneOrder, :);
uncorrectedGeneCorr = corr(tpmLogRelSorted', tpmLogRelSorted', 'Rows', 'pairwise');
correctedGeneCorr = corr(corrTpmLogRelSorted', corrTpmLogRelSorted', 'Rows', 'pairwise');

fig29 = plot_make_manual_figure();

axes('Parent', fig29, 'Units', 'centimeters', 'Position', plot_panel_position(1));
imagesc(timeVals, genePosSorted, tpmLogRelSorted);
set(gca, 'YDir', 'normal'); xlabel('Time (min)'); ylabel('Relative genome position');
title('Uncorrected RNA'); colormap(gca, plot_redblue(256)); caxis([-0.5, 0.5]); cb = colorbar; cb.Ticks = relSignalTicks;
yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig29, 'Units', 'centimeters', 'Position', plot_panel_position(2));
imagesc(genePosSorted, genePosSorted, uncorrectedGeneCorr);
axis image; set(gca, 'YDir', 'normal');
xlabel('Relative genome position'); ylabel('Relative genome position');
title('Uncorrected RNA correlation'); colormap(gca, plot_redblue(256)); caxis([-0.5, 0.5]); cb = colorbar; cb.Ticks = relSignalTicks;
xticks(relTicks); yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig29, 'Units', 'centimeters', 'Position', plot_panel_position(3));
imagesc(timeVals, genePosSorted, corrTpmLogRelSorted);
set(gca, 'YDir', 'normal'); xlabel('Time (min)'); ylabel('Relative genome position');
title('Copy-number-corrected RNA'); colormap(gca, plot_redblue(256)); caxis([-0.5, 0.5]); cb = colorbar; cb.Ticks = relSignalTicks;
yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

axes('Parent', fig29, 'Units', 'centimeters', 'Position', plot_panel_position(4));
imagesc(genePosSorted, genePosSorted, correctedGeneCorr);
axis image; set(gca, 'YDir', 'normal');
xlabel('Relative genome position'); ylabel('Relative genome position');
title('Corrected RNA correlation'); colormap(gca, plot_redblue(256)); caxis([-0.5, 0.5]); cb = colorbar; cb.Ticks = relSignalTicks;
xticks(relTicks); yticks(relTicks); plot_style_axes(gca, fontName, fontSizePt, axisLineWidthPt);

plot_save_figure(fig29, outDir, 'Figs_copy_number_correction');

%% Save intermediate matrices for inspection
fprintf('Writing intermediate matrices...\n');
writetable(array2table(dnaSampleCorr, 'VariableNames', matlab.lang.makeValidName("T" + timeNames)), ...
    fullfile(outDir, 'Figs_DNA_sample_correlation.csv'));
writetable(table(dnaInterval(:), dnaMean(:), dnaSd(:), ...
    'VariableNames', {'Sampling_interval_min', 'Mean_Pearson_r', 'SD'}), ...
    fullfile(outDir, 'Figs_DNA_interval_summary.csv'));
writetable(array2table(rnaSampleCorr, 'VariableNames', matlab.lang.makeValidName("T" + timeNames)), ...
    fullfile(outDir, 'Figs_RNA_sample_correlation.csv'));
writetable(table(rnaInterval(:), rnaMean(:), rnaSd(:), ...
    'VariableNames', {'Sampling_interval_min', 'Mean_Pearson_r', 'SD'}), ...
    fullfile(outDir, 'Figs_RNA_interval_summary.csv'));
writetable(table(dnaInterval(:), dnaMean(:), dnaSd(:), rnaMean(:), rnaSd(:), ...
    'VariableNames', {'Sampling_interval_min', 'DNA_mean_Pearson_r', 'DNA_SD', 'RNA_mean_Pearson_r', 'RNA_SD'}), ...
    fullfile(outDir, 'Figs_DNA_RNA_interval_summary.csv'));
writetable(binnedDPM, fullfile(outDir, 'Binned_DNA_copy_number.csv'));
writetable(binnedTPM, fullfile(outDir, 'Binned_TPM.csv'));
writetable(corrTPM, fullfile(outDir, 'Corrected_TPM_by_binned_DNA.csv'));
writetable(matchTable, fullfile(outDir, 'TPM_to_DNA_bin_match.csv'));

fprintf('Done. Outputs written to:\n%s\n', outDir);

end

%% Local functions
function [intervals, means, sds] = plot_summarize_corr_by_interval(corrMat, timeVals)
    intervals = unique(abs(timeVals(:) - timeVals(:)'));
    intervals = intervals(:);
    means = nan(size(intervals));
    sds = nan(size(intervals));
    for i = 1:numel(intervals)
        % A zero sampling interval represents each sample correlated with
        % itself. Its Pearson correlation is identically 1 with zero
        % dispersion. The previous implementation removed the diagonal,
        % leaving no values and therefore returning NaN.
        if intervals(i) == 0
            means(i) = 1;
            sds(i) = 0;
            continue;
        end
        mask = abs(timeVals(:) - timeVals(:)') == intervals(i);
        vals = corrMat(mask);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            means(i) = mean(vals);
            sds(i) = std(vals, 0);
        end
    end
end

function binned = plot_make_binned_matrix(df, timeVals, binSize, oriC, genomeSize)
    binCenter = [sort(oriC - (1:round(genomeSize / 2 / binSize)) * binSize, 'ascend'), ...
        oriC, oriC + (1:round(genomeSize / 2 / binSize)) * binSize - 1];
    binCenter(binCenter > genomeSize) = binCenter(binCenter > genomeSize) - genomeSize;
    binCenter(binCenter <= 0) = binCenter(binCenter <= 0) + genomeSize;
    binLeft = binCenter - binSize / 2;
    binRight = binCenter + binSize / 2 - 1;
    binLeft(binLeft <= 0) = binLeft(binLeft <= 0) + genomeSize;
    binRight(binRight > genomeSize) = binRight(binRight > genomeSize) - genomeSize;
    relPos = (binCenter - oriC) / (genomeSize / 2);
    relPos(relPos < -1) = relPos(relPos < -1) + 2;

    genePos = round((table2array(df(:, {'Start'})) + table2array(df(:, {'End'}))) / 2);
    values = table2array(df(:, 5:25));
    out = nan(numel(binCenter), numel(timeVals));
    for i = 1:numel(binCenter)
        if binLeft(i) <= binRight(i)
            inBin = genePos >= binLeft(i) & genePos <= binRight(i);
        else
            inBin = genePos >= binLeft(i) | genePos <= binRight(i);
        end
        if any(inBin)
            out(i, :) = mean(values(inBin, :), 1, 'omitnan');
        end
    end
    binned = array2table([binCenter(:), binLeft(:), binRight(:), relPos(:), out], ...
        'VariableNames', [{'Bin_center', 'Bin_start', 'Bin_end', 'Relative_position'}, string(timeVals)]);
end

function [corrTPM, matchTable] = plot_correct_tpm_by_binned_dpm(dfTPM, binnedDPM, timeVals)
    tpmData = table2array(dfTPM(:, 5:25));
    dnaData = table2array(binnedDPM(:, 5:25));
    dnaRel = dnaData ./ mean(dnaData, 2, 'omitnan');
    tpmPos = table2array(dfTPM(:, {'Relative position'}));
    binPos = table2array(binnedDPM(:, {'Relative_position'}));
    corrected = nan(size(tpmData));
    corrCoef = nan(height(dfTPM), 1);
    matchedBin = nan(height(dfTPM), 1);
    matchedBinPos = nan(height(dfTPM), 1);
    for i = 1:height(dfTPM)
        [~, idx] = min(abs(binPos - tpmPos(i)));
        matchedBin(i) = idx;
        matchedBinPos(i) = binPos(idx);
        corrected(i, :) = tpmData(i, :) ./ dnaRel(idx, :);
        corrCoef(i) = corr(tpmData(i, :)', dnaRel(idx, :)', 'Rows', 'complete');
    end
    corrTPM = [dfTPM(:, 1:4), array2table(corrected, 'VariableNames', string(timeVals)), dfTPM(:, 26:end)];
    if any(strcmp(corrTPM.Properties.VariableNames, 'corr_DNA'))
        corrTPM(:, {'corr_DNA'}) = array2table(corrCoef);
    end
    matchTable = table(table2array(dfTPM(:, {'ID'})), table2array(dfTPM(:, {'Gene name'})), tpmPos, matchedBin, matchedBinPos, ...
        'VariableNames', {'ID', 'Gene_name', 'TPM_relative_position', 'Matched_DNA_bin_index', 'Matched_DNA_bin_relative_position'});
end

function tbl = plot_remove_extra_annotation_columns(tbl)
    names = tbl.Properties.VariableNames;
    drop = startsWith(names, 'Unnamed') | startsWith(names, 'Var') | strcmp(names, 'CRP_regulated');
    tbl(:, drop) = [];
end

function cmap = plot_redblue(n)
    if nargin < 1
        n = 256;
    end
    bottom = [0.0, 0.2, 0.8];
    middle = [1, 1, 1];
    top = [0.85, 0.0, 0.0];
    n1 = floor(n / 2);
    n2 = n - n1;
    cmap = [interp1([1, n1], [bottom; middle], 1:n1); ...
        interp1([1, n2], [middle; top], 1:n2)];
end

function plot_style_axes(ax, fontName, fontSizePt, lineWidthPt)
    set(ax, 'FontName', fontName, 'FontSize', fontSizePt, 'LineWidth', lineWidthPt, 'Box', 'on', ...
        'TickDir', 'in', 'Layer', 'top');
    pbaspect(ax, [1, 1, 1]);
end

function fig = plot_make_manual_figure()
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 18.4, 15.4], 'PaperUnits', 'centimeters', ...
        'PaperPosition', [0, 0, 18.4, 15.4], 'PaperSize', [18.4, 15.4]);
end

function pos = plot_panel_position(panelIndex)
    axSize = 5.04;
    leftX = 1.8;
    rightX = 10.4;
    topY = 8.5;
    bottomY = 1.2;
    positions = [
        leftX, topY, axSize, axSize;
        rightX, topY, axSize, axSize;
        leftX, bottomY, axSize, axSize;
        rightX, bottomY, axSize, axSize
    ];
    pos = positions(panelIndex, :);
end

function plot_set_square_axes(fig, axisSizeCm)
    axesList = findall(fig, 'Type', 'axes');
    for i = 1:numel(axesList)
        axesList(i).Units = 'centimeters';
        pos = axesList(i).Position;
        axesList(i).Position = [pos(1), pos(2), axisSizeCm, axisSizeCm];
        pbaspect(axesList(i), [1, 1, 1]);
    end
end

function plot_save_figure(fig, outDir, baseName)
    pngPath = fullfile(outDir, [baseName, '.png']);
    pdfPath = fullfile(outDir, [baseName, '.pdf']);
    try
        exportgraphics(fig, pngPath, 'Resolution', 300);
        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    catch
        print(fig, pngPath, '-dpng', '-r300');
        print(fig, pdfPath, '-dpdf', '-painters');
    end
end

function export_sync_tables(outDir)



exportDir = outDir;
if ~isfolder(exportDir)
    mkdir(exportDir);
end

summary = readtable(fullfile(outDir, 'independent_fit_summary.csv'));
writetable(summary, fullfile(exportDir, 'best_parameters.csv'));

repNames = string(summary.Replicate);
idtGrid = [];
idtWide = table;

for iRep = 1:numel(repNames)
    rep = repNames(iRep);

    S = load(fullfile(outDir, rep + "_fit_workspace.mat"));
    G = readtable(fullfile(outDir, rep + "_parameter_grid.csv"));
    profile = groupsummary(G, 'IDT_min', 'min', 'RelativeDeviation');

    if iRep == 1
        idtGrid = profile.IDT_min;
        idtWide = table(idtGrid, 'VariableNames', {'IDT_min'});
    end
    idtWide.(char(rep)) = profile.min_RelativeDeviation;
    expT = S.exp_t(:);
    expY = S.exp_N(:) / S.exp_N(1);
    fitT = S.bestTDense(:);
    fitY = S.bestNDense(:) / S.bestNAtExperimentalTimes(1);

    nRow = max(numel(expT), numel(fitT));
    curveTable = table(nan(nRow,1), nan(nRow,1), nan(nRow,1), nan(nRow,1), ...
        'VariableNames', {'Experimental_Time_min', 'Experimental_NormalizedCellNo', ...
        'Fitted_Time_min', 'Fitted_NormalizedCellNo'});

    curveTable.Experimental_Time_min(1:numel(expT)) = expT;
    curveTable.Experimental_NormalizedCellNo(1:numel(expY)) = expY;
    curveTable.Fitted_Time_min(1:numel(fitT)) = fitT;
    curveTable.Fitted_NormalizedCellNo(1:numel(fitY)) = fitY;

    writetable(curveTable, fullfile(exportDir, rep + "_bestfit_curve.csv"));
    Gbest = G(G.IDT_min == S.bestIDT, :);
    cvGrid = unique(Gbest.CV);
    usGrid = unique(Gbest.UnsynchronizedFraction);
    heatmapMatrix = nan(numel(usGrid), numel(cvGrid));

    for iUS = 1:numel(usGrid)
        for iCV = 1:numel(cvGrid)
            idx = Gbest.UnsynchronizedFraction == usGrid(iUS) & ...
                  Gbest.CV == cvGrid(iCV);
            heatmapMatrix(iUS, iCV) = Gbest.RelativeDeviation(idx);
        end
    end

    heatmapTable = array2table(heatmapMatrix, ...
        'VariableNames', cellstr("CV_" + compose("%.2f", cvGrid)));
    heatmapTable = addvars(heatmapTable, usGrid, 'Before', 1, ...
        'NewVariableNames', 'UnsynchronizedFraction');

    writetable(heatmapTable, fullfile(exportDir, rep + "_heatmap.csv"));
end
writetable(idtWide, fullfile(exportDir, 'IDT_vs_minRD_all_replicates.csv'));

disp("Synchrony tables exported: " + exportDir);
end
