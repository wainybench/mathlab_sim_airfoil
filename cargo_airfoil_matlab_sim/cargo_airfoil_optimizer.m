
%% cargo_airfoil_optimizer.m
% Cargo aircraft + airfoil screening + 3S 1000 mAh battery range/power plot
%
% PURPOSE
%   1) Search wing geometry, speed, payload, and candidate airfoil polars.
%   2) Maximize INTERNAL payload while enforcing aerodynamic/propulsion limits.
%   3) Plot ELECTRICAL WATTS (Y) versus BATTERY RANGE (X, feet).
%
% IMPORTANT
%   - This model is a first-pass sizing/selection tool, not certification-level
%     analysis.
%   - Airfoil selection is only as good as the supplied polar data.
%   - Supply real XFOIL/wind-tunnel polar files in ./polars.
%   - Do NOT treat the placeholder geometric bounds below as final design
%     requirements; replace them with your actual project limits.
%
% POLAR CSV FORMAT (one file per airfoil)
%   Required columns:
%       Re, alpha_deg, CL, CD
%   Optional:
%       CM
%   Example filename:
%       polars/sd7037.csv
%
%   Data may contain multiple Reynolds numbers in the same file.
%   The script interpolates first in alpha and then in Re by choosing the
%   nearest available Re block. For final work, use sufficiently dense Re data.
%
% SOURCES FOR AIRFOIL COORDINATES/DATA:
%   UIUC Airfoil Data Site:
%   https://m-selig.ae.illinois.edu/ads/coord_database.html
%   UIUC low-speed airfoil test data:
%   https://m-selig.ae.illinois.edu/ads_faq.html
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% ========================= USER / MISSION INPUTS =========================
g = 9.80665;                % m/s^2
rho = 1.225;                % kg/m^3, sea-level standard density
mu = 1.81e-5;               % Pa*s, dynamic viscosity

% ---- Mass model ----
emptyMass_kg           = 1.20;
dragArticleMass_kg     = 1.20;
batteryMass_kg         = 0.150;     % 150 g
emptyMassIncludesBattery = false;   % change to true if the 1.20 kg already
                                    % includes the battery
internalPayloadRange_kg = 0:0.05:3.0;

% ---- Battery: 3S 1000 mAh ----
batteryCells            = 3;
batteryCapacity_Ah      = 1.000;
batteryNominalVoltage_V = 3.7*batteryCells;  % 11.1 V nominal
batteryEnergy_Wh       = batteryNominalVoltage_V*batteryCapacity_Ah;
batteryUsableFraction   = 0.80;       % reserve/de-rating assumption
usableEnergy_Wh         = batteryEnergy_Wh*batteryUsableFraction;

% ---- Propulsion/system assumptions ----
targetTW                = 0.50;       % required thrust-to-weight target
powertrainEfficiency    = 0.75;       % battery -> shaft/thrust useful efficiency
                                        % Replace with measured system efficiency.
useConstantTMax = true;                 % first-pass: T_max = targetTW * W
                                        % Better later: use measured T(V) curve.

% ---- Mission speed range ----
Vmin_ms       = 10;
Vcruise_ms    = 15;
Vmax_ms       = 25;
Vgrid_ms      = 10:0.5:25;

% ---- Wing design search (PLACEHOLDERS; replace with project limits) ----
Sgrid_m2      = 0.25:0.025:0.80;    % wing area
ARgrid        = 6:0.5:14;            % aspect ratio
maxSpan_m     = 2.00;                % maximum span
minChord_m    = 0.10;
maxChord_m    = 0.35;

% ---- Low-speed requirements ----
targetStall_ms = 8.0;                 % design constraint / placeholder
CL_margin      = 0.80;               % require CLreq <= 0.80*CLmax
Oswald_e       = 0.80;                % first-pass span efficiency

% ---- Structural placeholder ----
n_max = 3.0;                          % positive load factor
maxWingLoading_Npm2 = inf;            % set finite if required

% ---- Airfoil candidates ----
% Put matching CSV files inside ./polars.
candidateNames = { ...
    'SD7003', 'SD7037', 'S3003', 'S3021', ...
    'AG03', 'AG04', 'ClarkY', 'NACA2412'};

candidateFiles = { ...
    'polars/sd7003.csv', 'polars/sd7037.csv', 'polars/s3003.csv', ...
    'polars/s3021.csv', 'polars/ag03.csv', ...
    'polars/ag04.csv', 'polars/clarky.csv', 'polars/naca2412.csv'};

% Set true if you want to require a polar to exist before that airfoil is
% considered. This prevents accidental use of fake/default airfoil data.
requireRealPolarFiles = true;

%% ========================= VALIDATE INPUTS ================================
assert(~isempty(internalPayloadRange_kg), 'Payload range is empty.');
assert(batteryMass_kg > 0 && batteryCapacity_Ah > 0, 'Battery inputs invalid.');
assert(targetTW > 0, 'T/W target must be positive.');
assert(powertrainEfficiency > 0 && powertrainEfficiency <= 1, ...
    'Powertrain efficiency must be in (0,1].');

if emptyMassIncludesBattery
    fixedMass_kg = emptyMass_kg + dragArticleMass_kg;
else
    fixedMass_kg = emptyMass_kg + dragArticleMass_kg + batteryMass_kg;
end

%% ========================= LOAD AIRFOIL POLARS ============================
airfoils = struct('name', {}, 'file', {}, 'data', {});
for i = 1:numel(candidateNames)
    f = candidateFiles{i};
    if isfile(f)
        T = readtable(f);
        needed = {'Re','alpha_deg','CL','CD'};
        hasCols = all(ismember(needed, T.Properties.VariableNames));
        if ~hasCols
            warning('%s skipped: required columns Re, alpha_deg, CL, CD missing.', f);
            continue;
        end
        T = rmmissing(T, 'DataVariables', needed);
        airfoils(end+1).name = candidateNames{i}; %#ok<SAGROW>
        airfoils(end).file = f;
        airfoils(end).data = T;
    elseif requireRealPolarFiles
        fprintf('Polar file missing -> skipping %s: %s\n', candidateNames{i}, f);
    end
end

if isempty(airfoils)
    error(['No usable airfoil polar files were found. ' newline ...
           'Add CSV polar files to ./polars using columns: Re, alpha_deg, CL, CD' newline ...
           'and rerun.']);
end

fprintf('\nLoaded %d airfoil polar file(s).\n', numel(airfoils));

%% ========================= OPTIMIZATION SEARCH ===========================
%
% Objective:
%   maximize internal payload.
%
% Tie-breakers among equal payload:
%   lower cruise power, then higher L/D.
%
% Constraints:
%   L >= W
%   CLreq <= CLmax*CL_margin
%   Vs <= targetStall_ms
%   D <= T_available
%   b <= maxSpan
%   chord in bounds
%   Re inside available polar data
%   wing loading <= specified maximum, if finite
%   nmax used for a structural feasibility screen

results = table();

for ia = 1:numel(airfoils)
    af = airfoils(ia);
    data = af.data;
    Re_available = unique(data.Re);

    bestForAirfoil = [];

    for mip = internalPayloadRange_kg
        totalMass = fixedMass_kg + mip;
        W = totalMass*g;

        for S = Sgrid_m2
            if W/S > maxWingLoading_Npm2
                continue;
            end

            for AR = ARgrid
                b = sqrt(AR*S);
                if b > maxSpan_m
                    continue;
                end

                c = S/b;
                if c < minChord_m || c > maxChord_m
                    continue;
                end

                for V = Vgrid_ms
                    q = 0.5*rho*V^2;
                    CLreq = W/(q*S);
                    Re = rho*V*c/mu;

                    if Re < min(Re_available) || Re > max(Re_available)
                        continue;
                    end

                    % Get 2-D airfoil polar point near required CL at this Re.
                    [CD2D, alphaDeg, CL_used, ok, CLmax_local] = ...
                        interpPolarAtCL(data, Re, CLreq);

                    if ~ok
                        continue;
                    end

                    if CLreq > CL_margin*CLmax_local
                        continue;
                    end

                    % Wing-induced drag
                    CDi = CLreq^2/(pi*Oswald_e*AR);
                    CDtotal = CD2D + CDi;
                    D = q*S*CDtotal;

                    % Available thrust based on target T/W.
                    Tavailable = targetTW*W;

                    if D > Tavailable
                        continue;
                    end

                    % Stall estimate using airfoil CLmax
                    Vs = sqrt(2*W/(rho*S*CLmax_local));
                    if Vs > targetStall_ms
                        continue;
                    end

                    % Cruise L/D
                    LD = CLreq/CDtotal;

                    % Electrical power needed at cruise.
                    % P_shaft = D*V; P_batt = P_shaft/eta.
                    P_batt_W = D*V/powertrainEfficiency;

                    % Battery-only range at this condition.
                    usableEnergy_J = usableEnergy_Wh*3600;
                    endurance_s = usableEnergy_J/P_batt_W;
                    range_m = V*endurance_s;
                    range_ft = range_m*3.280839895;

                    candidate = struct( ...
                        'Airfoil', string(af.name), ...
                        'InternalPayload_kg', mip, ...
                        'TotalMass_kg', totalMass, ...
                        'WingArea_m2', S, ...
                        'AR', AR, ...
                        'Span_m', b, ...
                        'Chord_m', c, ...
                        'Speed_ms', V, ...
                        'Re', Re, ...
                        'Alpha_deg', alphaDeg, ...
                        'CL', CL_used, ...
                        'CLmax', CLmax_local, ...
                        'CD2D', CD2D, ...
                        'CDi', CDi, ...
                        'CDtotal', CDtotal, ...
                        'L_N', W, ...
                        'D_N', D, ...
                        'Tavail_N', Tavailable, ...
                        'Tmargin_N', Tavailable-D, ...
                        'Vs_ms', Vs, ...
                        'LD', LD, ...
                        'Power_W', P_batt_W, ...
                        'Range_ft', range_ft);

                    % Select the best feasible case for this airfoil.
                    if isempty(bestForAirfoil)
                        bestForAirfoil = candidate;
                    else
                        if isBetterCandidate(candidate, bestForAirfoil)
                            bestForAirfoil = candidate;
                        end
                    end

                    results = [results; struct2table(candidate)]; %#ok<AGROW>
                end
            end
        end
    end

    if ~isempty(bestForAirfoil)
        fprintf(['Best feasible point for %-8s: payload = %.2f kg, ' ...
                 'S = %.3f m^2, AR = %.1f, V = %.1f m/s, ' ...
                 'P = %.1f W, range = %.0f ft, L/D = %.1f\n'], ...
            bestForAirfoil.Airfoil, bestForAirfoil.InternalPayload_kg, ...
            bestForAirfoil.WingArea_m2, bestForAirfoil.AR, ...
            bestForAirfoil.Speed_ms, bestForAirfoil.Power_W, ...
            bestForAirfoil.Range_ft, bestForAirfoil.LD);
    else
        fprintf('No feasible design found for %s.\n', af.name);
    end
end

if isempty(results)
    error('No feasible wing/airfoil combinations found. Expand the design bounds or relax constraints.');
end

%% ========================= RANK AIRFOILS =================================
% For each airfoil, keep its maximum-payload feasible point.
bestRows = table();
names = unique(results.Airfoil);

for i = 1:numel(names)
    mask = results.Airfoil == names(i);
    R = results(mask,:);
    maxPayload = max(R.InternalPayload_kg);
    R = R(R.InternalPayload_kg == maxPayload,:);

    % Tie-breaker: minimum power, then maximum L/D.
    [~, idx] = min(R.Power_W + 0.001./max(R.LD, eps));
    bestRows = [bestRows; R(idx,:)]; %#ok<AGROW>
end

bestRows = sortrows(bestRows, {'InternalPayload_kg','Power_W'}, {'descend','ascend'});

disp(' ');
disp('================ AIRFOIL / DESIGN RANKING ================');
disp(bestRows(:, {'Airfoil','InternalPayload_kg','TotalMass_kg', ...
    'WingArea_m2','AR','Span_m','Chord_m','Speed_ms','Re', ...
    'CL','CLmax','CDtotal','LD','Power_W','Range_ft'}));

writetable(bestRows,'airfoil_ranking.csv');

%% ========================= GRAPH 1: AIRFOIL SCREENING =====================
figure('Name','Airfoil screening','Color','w');
scatter(bestRows.InternalPayload_kg, bestRows.Power_W, 70, 'filled');
grid on;
xlabel('Maximum feasible internal payload (kg)');
ylabel('Electrical cruise power (W)');
title('Airfoil screening: payload capability vs cruise power');

for i = 1:height(bestRows)
    text(bestRows.InternalPayload_kg(i), bestRows.Power_W(i), ...
        "  " + bestRows.Airfoil(i), 'Interpreter','none');
end

%% ========================= GRAPH 2: POWER VS BATTERY RANGE ================
% Use the maximum payload that is feasible for ALL airfoils being compared.
commonPayload = min(bestRows.InternalPayload_kg);

fprintf('\nCommon payload used for the power-vs-range comparison: %.2f kg\n', commonPayload);

figure('Name','Battery power vs range','Color','w');
hold on; grid on;

legendEntries = strings(0);

for i = 1:height(bestRows)
    afName = bestRows.Airfoil(i);

    % Select the geometry associated with this airfoil at common payload.
    candidates = results(results.Airfoil == afName & ...
                         results.InternalPayload_kg == commonPayload,:);

    if isempty(candidates)
        continue;
    end

    % Use the design with lowest cruise power at the nominal cruise speed.
    cruiseCandidates = candidates(abs(candidates.Speed_ms - Vcruise_ms) < 1e-9,:);
    if isempty(cruiseCandidates)
        cruiseCandidates = candidates;
    end
    [~, k] = min(cruiseCandidates.Power_W);
    base = cruiseCandidates(k,:);

    Vcurve = linspace(Vmin_ms,Vmax_ms,80);
    Pcurve = nan(size(Vcurve));
    Rcurve = nan(size(Vcurve));

    for j = 1:numel(Vcurve)
        V = Vcurve(j);
        q = 0.5*rho*V^2;
        CLreq = base.L_N/(q*base.WingArea_m2);
        Re = rho*V*base.Chord_m/mu;

        [CD2D, ~, ~, ok, CLmax_local] = interpPolarAtCL( ...
            airfoils(strcmp({airfoils.name},char(afName))).data, Re, CLreq);

        if ~ok || CLreq > CL_margin*CLmax_local
            continue;
        end

        CDi = CLreq^2/(pi*Oswald_e*base.AR);
        CDtotal = CD2D + CDi;
        D = q*base.WingArea_m2*CDtotal;

        if D > targetTW*base.L_N
            continue;
        end

        P = D*V/powertrainEfficiency;
        E_J = usableEnergy_Wh*3600;
        Rm = V*(E_J/P);

        Pcurve(j) = P;
        Rcurve(j) = Rm*3.280839895;
    end

    valid = isfinite(Pcurve) & isfinite(Rcurve);
    plot(Rcurve(valid), Pcurve(valid), 'LineWidth', 1.5);
    legendEntries(end+1) = afName; %#ok<AGROW>

    % Mark nominal cruise point
    idxV = find(abs(Vcurve - Vcruise_ms) == min(abs(Vcurve - Vcruise_ms)),1);
    if valid(idxV)
        plot(Rcurve(idxV),Pcurve(idxV),'o','HandleVisibility','off');
    end
end

xlabel('Battery range per charge (ft)');
ylabel('Electrical power required (W)');
title(sprintf(['3S 1000 mAh battery: Power vs range at %.2f kg common internal payload\n' ...
               'Battery mass = %.0f g, usable energy = %.2f Wh'], ...
               commonPayload, batteryMass_kg*1000, usableEnergy_Wh));
legend(legendEntries,'Location','best','Interpreter','none');

%% ========================= GRAPH 3: PAYLOAD VS POWER =======================
figure('Name','Payload vs power','Color','w');
hold on; grid on;
for i = 1:height(bestRows)
    afName = bestRows.Airfoil(i);
    mask = results.Airfoil == afName;
    R = results(mask,:);
    % Keep best power at each payload.
    payloads = unique(R.InternalPayload_kg);
    pbest = nan(size(payloads));
    for j = 1:numel(payloads)
        rr = R(R.InternalPayload_kg == payloads(j),:);
        pbest(j) = min(rr.Power_W);
    end
    plot(payloads,pbest,'LineWidth',1.5);
end
xlabel('Internal payload (kg)');
ylabel('Minimum feasible electrical cruise power (W)');
title('Payload vs minimum cruise power');
legend(bestRows.Airfoil,'Location','best','Interpreter','none');

%% ========================= SUMMARY ========================================
fprintf('\n================ BATTERY SUMMARY ================\n');
fprintf('Nominal voltage: %.2f V\n', batteryNominalVoltage_V);
fprintf('Capacity: %.2f Ah\n', batteryCapacity_Ah);
fprintf('Nominal energy: %.2f Wh\n', batteryEnergy_Wh);
fprintf('Usable energy assumption: %.2f Wh (%.0f%%)\n', ...
    usableEnergy_Wh,100*batteryUsableFraction);
fprintf('Battery mass: %.0f g\n', batteryMass_kg*1000);

fprintf('\n================ MODEL NOTES ================\n');
fprintf(['1) Replace placeholder S/AR/speed/stall/span bounds with your project values.\n' ...
         '2) Replace powertrainEfficiency with measured motor+ESC+propulsive efficiency.\n' ...
         '3) For final propulsion feasibility, replace constant T/W with measured T(V).\n' ...
         '4) Add cargo volume and CG constraints once cargo dimensions are known.\n' ...
         '5) Add structural/FEA results before treating the maximum payload as a design value.\n']);

%% ========================= LOCAL FUNCTIONS ================================
function [CD, alphaDeg, CLout, ok, CLmax] = interpPolarAtCL(T, ReTarget, CLTarget)
% Interpolate within the nearest available Reynolds-number polar block.

    ReVals = unique(T.Re);
    [~, idxRe] = min(abs(ReVals-ReTarget));
    ReUse = ReVals(idxRe);

    R = T(abs(T.Re-ReUse) < max(1e-12,0.001*ReUse), :);
    R = sortrows(R, 'CL');

    if height(R) < 4
        CD = NaN; alphaDeg = NaN; CLout = NaN; ok = false; CLmax = NaN;
        return;
    end

    CLmax = max(R.CL);

    % Do not interpolate beyond measured/predicted data.
    if CLTarget < min(R.CL) || CLTarget > max(R.CL)
        CD = NaN; alphaDeg = NaN; CLout = NaN; ok = false;
        return;
    end

    % CL can repeat after stall, so use only the pre-stall branch up to CLmax.
    % For a robust first pass, choose the point with the smallest CD among
    % samples close to the target CL.
    [~, order] = sort(R.CL);
    R = R(order,:);

    CD = interp1(R.CL, R.CD, CLTarget, 'linear', NaN);
    if ismember('alpha_deg', R.Properties.VariableNames)
        alphaDeg = interp1(R.CL, R.alpha_deg, CLTarget, 'linear', NaN);
    else
        alphaDeg = NaN;
    end

    CLout = CLTarget;
    ok = isfinite(CD) && isfinite(CLout);
end

function tf = isBetterCandidate(a,b)
% Primary: maximum internal payload.
% Tie-breaker 1: lower power.
% Tie-breaker 2: higher L/D.

    tol = 1e-9;

    if a.InternalPayload_kg > b.InternalPayload_kg + tol
        tf = true;
        return;
    end

    if abs(a.InternalPayload_kg-b.InternalPayload_kg) <= tol
        if a.Power_W < b.Power_W - 1e-6
            tf = true;
            return;
        elseif abs(a.Power_W-b.Power_W) <= 1e-6
            tf = a.LD > b.LD;
            return;
        end
    end

    tf = false;
end
