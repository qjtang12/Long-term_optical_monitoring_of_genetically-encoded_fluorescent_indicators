function output = generateFiberPhotometryTraces(preset, doPlot)
% generateFiberPhotometryTraces - Generates and optionally plots synthetic fiber photometry traces.
%
% Usage:
%   % Generate data AND create a plot (default behavior)
%   output = generateFiberPhotometryTraces('big_sine');
%
%   % Generate data only, without plotting
%   output = generateFiberPhotometryTraces('small_sine', false);
%
% Inputs:
%   preset - (Optional) 'big_sine' or 'small_sine'. Defaults to 'big_sine'.
%   doPlot - (Optional) true or false. If true, a plot is generated. Defaults to true.
%
% Outputs:
%   output - A structure containing the generated signals and metadata.

    %% --------------------------------------------------------------------
    %  Parameter Configuration
    %% --------------------------------------------------------------------
    if nargin < 1 || isempty(preset)
        preset = 'big_sine';
    end
    if nargin < 2
        doPlot = true; % Plot by default
    end
    
    params = configureParameters(preset);

    % Add the chosen preset to the params struct for use in the plot title
    params.preset = preset;

    if ~isempty(params.randomSeed)
        rng(params.randomSeed);
    else
        rng('shuffle');
    end

    %% --------------------------------------------------------------------
    %  Time Base and Slow Components
    %% --------------------------------------------------------------------
    t = (1:params.nSamples)';
    sineWave = params.sineAmplitude * sin(2 * pi * t / params.sinePeriodSamples);
    expDecay = exp(-t / params.decayTauSamples);

    %% --------------------------------------------------------------------
    %  Segmented Baselines for Each Channel
    %% --------------------------------------------------------------------
    segmentLength = params.nSamples / params.nSegments;
    segBaseCh1 = zeros(params.nSamples, 1);
    segBaseCh2 = zeros(params.nSamples, 1);

    for seg = 1:params.nSegments
        startIdx = round((seg - 1) * segmentLength) + 1;
        endIdx   = min(round(seg * segmentLength), params.nSamples);
        midPoint = round(mean([startIdx, endIdx]));
        segBaseCh1(startIdx:endIdx) = sineWave(midPoint);
        segBaseCh2(startIdx:endIdx) = params.ch2Baseline + params.ch2Drift * randn();
    end

    segBaseCh1 = segBaseCh1 .* expDecay;
    segBaseCh2 = segBaseCh2 .* expDecay;

    %% --------------------------------------------------------------------
    %  Initialize Signals and Add Noise
    %% --------------------------------------------------------------------
    ch1 = segBaseCh1 + params.noiseStd * randn(params.nSamples, 1);
    ch2 = segBaseCh2 + params.noiseStd * randn(params.nSamples, 1);

    %% --------------------------------------------------------------------
    %  Add Exponential Peaks to Channel 1
    %% --------------------------------------------------------------------
    for p = 1:params.nPeaks
        peakPlaced = false;
        while ~peakPlaced
            startIdx = randi([1, params.nSamples - max(params.peakWidthRange)]);
            localSine   = sineWave(startIdx);
            phaseWeight = (localSine / params.sineAmplitude + 1) / 2;
            probability = params.peakBaseProb * phaseWeight;

            if rand < probability
                peakPlaced = true;
                width   = randi(params.peakWidthRange);
                peakAmp = params.peakAmplitudeMu + params.peakAmplitudeSd * randn();
                tau     = -width / log(params.peakDecayLevel);

                endIdx = min(startIdx + width, params.nSamples);
                idxRange = startIdx:endIdx;
                timeKernel = (0:length(idxRange)-1)';
                decayKernel = peakAmp * exp(-timeKernel / tau);
                ch1(idxRange) = ch1(idxRange) + decayKernel;
            end
        end
    end

    %% --------------------------------------------------------------------
    %  Add Shared Blips (Artifacts) to Both Channels
    %% --------------------------------------------------------------------
    for b = 1:params.nBlips
        startIdx = randi([1, params.nSamples - max(params.blipWidthRange)]);
        width    = randi(params.blipWidthRange);
        endIdx   = min(startIdx + width - 1, params.nSamples);
        amp = params.blipAmpRange(1) + rand * diff(params.blipAmpRange);
        idxRange = startIdx:endIdx;
        shape = amp * linspace(1, 0, length(idxRange))';
        ch1(idxRange) = ch1(idxRange) + shape;
        ch2(idxRange) = ch2(idxRange) + shape;
    end

    %% --------------------------------------------------------------------
    %  Prepare Output and Optionally Plot
    %% --------------------------------------------------------------------
    output.signal1 = ch1;
    output.signal2 = ch2;
    output.timeVector = t;
    output.baselines.signal1 = segBaseCh1;
    output.baselines.signal2 = segBaseCh2;
    output.params = params;
    
    % Call the local plotting function if requested
    if doPlot
        plotTraces(output);
    end
end


function params = configureParameters(preset)
% (Local Function) Configures parameters for the simulation.
    params = struct();

    % --- Base parameters common to all presets ---
    params.nSamples = 100000;           % Total number of data points in the trace.
    params.daysToSimulate = 3;          % The total duration to simulate, used for x-axis scaling.
    params.randomSeed = 42;             % A fixed seed for reproducible random numbers. Set to [] for a different trace each time.
    params.nSegments = 96;              % Number of flat steps used to approximate the sine wave baseline. (96 segments / 3 days = 32 segments/day).
    params.sinePeriodSamples = params.nSamples / params.daysToSimulate; % Sets the sine wave period to exactly 1 day.
    params.decayTauSamples = params.nSamples * 2; % Time constant for photobleaching. A large value creates a gentle decay over the full trace.
    params.noiseStd = 0.01;             % Standard deviation of the additive Gaussian noise.
    params.ch2Baseline = 0.15;          % The mean baseline value for the isosbestic (control) channel.
    params.ch2Drift = 0.01;             % The amount of random vertical drift between segments in the control channel.
    params.nBlips = 10;                 % The total number of motion artifact "blips" to add to the trace.
    params.peakDecayLevel = 0.01;       % Defines how much a peak decays by its end (to 1% of its initial amplitude).

    % --- Preset-specific parameters ---
    switch lower(preset)
        case 'big_sine'
            % --- Settings for a trace with a prominent circadian rhythm ---

            % Sine Wave Modulation
            params.sineAmplitude = 0.5;         % Controls the max/min height of the sine wave. A large value creates a strong, visible circadian modulation.

            % Calcium-like Peaks (Signal 1)
            params.nPeaks = 500;                % Total number of calcium-like peaks to attempt to place.
            params.peakWidthRange = [100, 500]; % Defines the min/max duration (in samples) of each peak. Wider peaks for a more pronounced look.
            params.peakAmplitudeMu = 0.3;       % The mean height of the peaks.
            params.peakAmplitudeSd = 0.05;      % The standard deviation of peak height, adding variability.
            params.peakBaseProb = 0.05;         % A base probability factor for placing a peak. This is multiplied by the sine phase, so peaks are more likely at the top of the sine wave.

            % Motion Artifacts / Blips (Both Signals)
            params.blipWidthRange = [300, 800]; % The min/max duration of the motion artifacts.
            params.blipAmpRange = [-0.2, -0.1]; % The min/max depth of the negative blips. Larger negative values create more dramatic dips.

        case 'small_sine'
            % --- Settings for a trace with a subtle circadian rhythm ---

            % Sine Wave Modulation
            params.sineAmplitude = 0.05;        % A small value creates a subtle, almost flat baseline modulation.

            % Calcium-like Peaks (Signal 1)
            params.nPeaks = 800;                % More peaks are added to give the trace a busier, more continuously active look.
            params.peakWidthRange = [80, 200];  % Peaks are narrower, making them look sharper and more "spiky".
            params.peakAmplitudeMu = 0.15;      % The mean height of the peaks is smaller.
            params.peakAmplitudeSd = 0.05;      % Variability in peak height is kept the same.
            params.peakBaseProb = 0.05;         % The same base probability, but since the sine amplitude is small, the influence of the sine phase on peak placement is much weaker.

            % Motion Artifacts / Blips (Both Signals)
            params.blipWidthRange = [100, 200]; % Artifacts are shorter in duration.
            params.blipAmpRange = [-0.15, -0.05];% The negative dips are shallower and less pronounced.

        otherwise
            error('Unknown preset "%s". Use ''big_sine'' or ''small_sine''.', preset);
    end
end


function plotTraces(data)
% (Local Function) Plots the generated fiber photometry traces.
    figure('Color', 'w', 'Position', [100, 100, 900, 700]);
    time_days = data.timeVector / data.params.sinePeriodSamples;

    subplot(2, 1, 1);
    plot(time_days, data.signal1, 'Color', [0.2, 0.6, 1], 'LineWidth', 1);
    hold on;
    plot(time_days, data.signal2, 'Color', [0.8, 0, 0.8], 'LineWidth', 1.5);
    hold off;
    title(sprintf('Simulated Photometry Signals (Preset: %s)', data.params.preset), 'Interpreter', 'none');
    xlabel('Time (days)');
    ylabel('Fluorescence (A.U.)');
    legend({'Calcium-dependent', 'Isosbestic'}, 'Location', 'northeast');
    grid on; axis tight; set(gca, 'FontSize', 12);

    subplot(2, 1, 2);
    plot(time_days, data.baselines.signal1, 'k-', 'LineWidth', 1.5);
    hold on;
    plot(time_days, data.baselines.signal2, 'r-', 'LineWidth', 1.5);
    hold off;
    title('Underlying Baselines with Decay');
    xlabel('Time (days)');
    ylabel('Baseline (A.U.)');
    legend({'Channel 1 Baseline', 'Channel 2 Baseline'}, 'Location', 'northeast');
    grid on; axis tight; set(gca, 'FontSize', 12);
end