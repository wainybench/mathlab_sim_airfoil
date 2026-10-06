%% cargo_flight_sim.m
% Time-domain mission simulation for the cargo aircraft sized by
% cargo_airfoil_optimizer.m
%
% WHAT THIS ADDS OVER THE STEADY-STATE OPTIMIZER
%   - 3-DOF point-mass flight dynamics: x, y, h, V, gamma (flight path), psi (heading)
%   - Takeoff ground roll: rolling friction, ground effect, rotation at Vrot
%   - Climb to altitude, then racetrack laps with banked turns
%     (turn load factor raises CL and drag, so turns cost energy)
%   - Propulsion, two modes:
%       'model' : motor (Kv, Rm, I0) + propeller CT(J), CP(J), solved analytically
%       'table' : measured thrust-stand data T(V, throttle), I(V, throttle)
%   - LiPo battery: OCV(SOC) curve + internal resistance, so voltage sags under load
%   - Uses only the pre-stall branch of the polar
%   - Logs Re at every step and reports how much of the flight lies inside polar data
%   - Payload sweep: largest payload that passes the full mission
%
% REQUIRES: MATLAB R2019b+ (tiledlayout, local functions in scripts)
%
% THRUST TABLE CSV (for cfg.prop.mode = 'table'):
%   columns V_ms, throttle (0..1), thrust_N, current_A, measured at cfg.prop.tableVtest
% -------------------------------------------------------------------------

clear; clc; close all;

%% =============================== CONFIG ==================================
cfg = struct();

% ---- Environment ----
cfg.env.g    = 9.80665;
cfg.env.rho  = 1.225;
cfg.env.mu   = 1.81e-5;
cfg.env.wind = [0 0];            % air-mass velocity [wx wy] m/s; runway along +x

% ---- Pull design from optimizer output? ----
cfg.useOptimizerResult = false;  % true -> top row of airfoil_ranking.csv
cfg.optimizerFile      = 'airfoil_ranking.csv';

% ---- Aircraft / aerodynamics ----
cfg.ac.airfoil        = 's1210';
cfg.ac.polarFile      = 'polars/s1210_Re100k_200K.csv';
cfg.ac.span_m  = 1.49;   % FIXED wingspan (m)  [= 77.5 in]
cfg.ac.chord_m = 0.2;   % DESIGN VARIABLE: constant (rectangular) wing chord (m)
% S and AR are filled in by deriveWing():  S = span*chord,  AR = span/chord
cfg.ac.e              = 0.80;
cfg.ac.CD0_airframe   = 0.030;   % fuselage+tail+gear, referenced to S  (PLACEHOLDER)
cfg.ac.CDA_extra      = 0.000;   % m^2 extra drag area (external stores etc.)
cfg.ac.CLmax3D_factor = 0.90;    % wing CLmax ~ 0.9 * airfoil CLmax
cfg.ac.CL_clampFrac   = 0.95;    % autopilot never commands above this * CLmax3D
cfg.ac.CL_margin      = 0.80;    % optimizer's margin (flagged, not enforced)
cfg.ac.CL_groundRoll  = 0.40;    % CL in ground-roll attitude
cfg.ac.mu_roll        = 0.05;    % rolling friction (pavement ~0.04, grass ~0.08)
cfg.ac.wingHeight_m   = 0.15;    % wing height above ground on gear
cfg.ac.nMax           = 3.0;     % structural load-factor limit
cfg.ac.ReCDexp        = -0.5; % CD2D ~ Re^exp outside polar data (0 = off)
cfg.ac.CL_floor = 0.0; % control limit, independent of polar coverage
cfg.mass.wingKgPerM2 = 0.0;    % kg of structure per m^2 of wing area (MEASURE from your build)
cfg.mass.S_ref       = 0.323;  % wing area that empty_kg was estimated for
% ---- Mass ----l
cfg.mass.empty_kg        = 1.20;
cfg.mass.dragArticles_kg = 1.473;
cfg.ac.dragArt.V = [7 10];
cfg.ac.dragArt.F = [0.4412 0.8939];   % N, Pig + Banana + Chicken
cfg.mass.battery_kg      = 0.150;
cfg.mass.payload_kg      = 1.2;

% ---- Battery: 3S 1000 mAh LiPo ----
cfg.batt.cells        = 3;
cfg.batt.capacity_Ah  = 1.000;
cfg.batt.R_cell       = 0.015;   % ohm per cell (MEASURE: DC IR from charger)
cfg.batt.R_wiring     = 0.005;   % ohm, leads + connectors
cfg.batt.SOC0         = 1.00;
cfg.batt.reserveSOC   = 0.20;    % mission ends here (= your 80% usable)
cfg.batt.cutoffCell_V = 3.30;    % loaded cell-voltage floor
cfg.batt.maxC         = 80;      % continuous C rating
cfg.batt.ocvSOC  = [0    0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00];
cfg.batt.ocvCell = [3.27 3.61 3.69 3.73 3.77 3.79 3.82 3.87 3.92 3.98 4.08 4.20];

% ---- Propulsion ----
cfg.prop.mode       = 'ba1130';   % 'model' or 'table'
cfg.prop.Kv = 1130;
cfg.prop.Rm = 0.014;
cfg.prop.I0 = 2.30;      % datasheet value at 10 V
cfg.prop.etaESC     = 0.95;
cfg.prop.D_m = 12*0.0254;
cfg.prop.CT  = [0.10206 -0.03144 -0.17420];
cfg.prop.CP  = [0.0325  0.0637 -0.1676];
cfg.prop.CtFile = 'prop/apce_12x6.csv';
cfg.prop.tableFile  = 'thrust/thrust_table.csv';
cfg.prop.tableVtest = 11.1;      % pack voltage during thrust-stand test

% ---- buildModels(), add: ----
if isfield(cfg.prop,'ctFile') && ~isempty(cfg.prop.ctFile)
    Pt = readtable(cfg.prop.ctFile);
    M.fCT = griddedInterpolant(Pt.J, Pt.CT, 'pchip', 'linear');
else
    M.fCT = @(J) cfg.prop.CT(1) + cfg.prop.CT(2)*J + cfg.prop.CT(3)*J.^2;
end

% ---- propulsion(): pass M through ----
[T, Im, rpm] = motorProp(thr*Vt, V, cfg, M);

% ---- motorProp(): new signature, replace the thrust line ----
function [T, Im, rpm] = motorProp(Vm, V, ~, M)
    ...                                   % rpm solve unchanged (uses CP quadratic)
    J   = V/(n*D);
    T   = rho*n^2*D^4*M.fCT(J);           % identical to old formula when fCT is the quadratic
    Im  = max((Vm - Kt*2*pi*n)/p.Rm, 0);
    rpm = 60*n;
end

% ---- Mission ----
cfg.mis.dt              = 0.02;    % s
cfg.mis.tMax            = 1800;    % s
cfg.mis.targetTW        = 0.50;
cfg.mis.runway_m        = 15;   % 100 ft                    (PLACEHOLDER: use rules)
cfg.mis.VrotFactor      = 1.15;    % Vrot = factor * Vstall
cfg.mis.VminFactor = 1.3;     % cruise speed is never below VminFactor * Vstal
cfg.mis.cruiseAlt_m     = 30;
cfg.mis.climbGamMax_deg = 10;
cfg.mis.Vcruise_ms      = 11;
cfg.mis.legLength_m     = 152.4;   % 500 ft straights
cfg.mis.bank_deg        = 30;
cfg.mis.nLaps           = 4;     % inf -> fly until reserve/cutoff
cfg.mis.requiredLaps    = 1;       % pass/fail criterion
cfg.prop.thrMax = 1;

% ---- Autopilot gains ----
cfg.ap.kV_p   = 0.15;   % throttle per m/s speed error
cfg.ap.kV_i   = 0.05;   % throttle per (m/s * s)
cfg.ap.kh     = 0.04;   % rad gamma command per m altitude error
cfg.ap.kgam   = 2.0;    % 1/s flight-path tracking
cfg.ap.kpsi   = 1.5;    % rad bank per rad heading error
cfg.ap.kVprot = 0.05;   % rad gamma reduction per m/s underspeed in climb
cfg.ap.phiMaxStraight_deg = 15;

% ---- Payload sweep ----
cfg.sweep.enable   = true;
cfg.sweep.payloads = 0:0.05:1.50;
cfg.sweep.dt       = 0.05; % coarser step for speed
cfg.sweep.chordEnable = true;
cfg.sweep.chordStart  = 0.12;   % m  first chord in the sweep
cfg.sweep.chordEnd    = 0.28;   % m  last chord in the sweep
cfg.sweep.chordStep   = 0.01;   % m  step between chords
cfg.sweep.chords      = cfg.sweep.chordStart:cfg.sweep.chordStep:cfg.sweep.chordEnd;
% (same idea as  cfg.sweep.payloads = 0:0.05:1.50 ; for a fixed NUMBER of points use
%  cfg.sweep.chords = linspace(cfg.sweep.chordStart, cfg.sweep.chordEnd, 15); )
cfg.sweep.chordNLaps  = inf;    % inf = fly to battery reserve (range)

%% ================================ RUN ====================================
if cfg.useOptimizerResult
    cfg = applyOptimizerResult(cfg);
end

cfg = deriveWing(cfg);
M = buildModels(cfg);
fprintf('Polar: %s | Re block(s): %s\n', cfg.ac.polarFile, mat2str(M.polar.Re));

R = runMission(cfg, M);
printSummary(R, cfg);
selfCheck(R, cfg);
%%    and, once, to test the time step:        
convergenceCheck(cfg, M);
plotMission(R, cfg, M);

if cfg.sweep.enable
    SW = payloadSweep(cfg, M);
    if cfg.sweep.chordEnable
        CS = chordSweep(cfg, M);
        writetable(CS, 'chord_sweep.csv');
    end
    writetable(SW, 'payload_sweep.csv');
end

%% =========================== MISSION LOOP ================================
function R = runMission(cfg, M)
    g = cfg.env.g;  rho = cfg.env.rho;  dt = cfg.mis.dt;
    N = floor(cfg.mis.tMax/dt) + 1;

    mTot  = cfg.mass.empty_kg + cfg.mass.wingKgPerM2*(cfg.ac.S - cfg.mass.S_ref) + ...
            cfg.mass.dragArticles_kg + cfg.mass.battery_kg + cfg.mass.payload_kg;
    W     = mTot*g;
    S     = cfg.ac.S;
    span  = sqrt(cfg.ac.AR*S);
    chord = S/span;

    Vs   = sqrt(2*W/(rho*S*cfg.ac.CLmax3D_factor*M.CLmax2Dref));
    Vrot = cfg.mis.VrotFactor*Vs;
    Vcmd = max(cfg.mis.Vcruise_ms, cfg.mis.VminFactor*Vs);
    T0   = propulsion(1, 0, cfg.batt.SOC0, cfg, M);

    s = zeros(9,1);              % [x y h V gam psi Ah_used Wh_used dist]
    phase = 1;                   % 1 takeoff roll, 2 climb, 3 laps
    onGround = true;
    iV = 1;                      % speed-loop integrator (throttle units)
    legMode = 1;                 % 1 straight, 2 turn
    psiRef = 0; legStart = 0; turned = 0; psiPrev = 0;
    turns = 0; laps = 0;
    groundRoll = NaN; tLiftoff = NaN;
    status = 'END: tMax reached';

    hCmd    = cfg.mis.cruiseAlt_m;
    gamMax  = deg2rad(cfg.mis.climbGamMax_deg);
    bank    = deg2rad(cfg.mis.bank_deg);
    phiMaxS = deg2rad(cfg.ap.phiMaxStraight_deg);
    ap      = cfg.ap;
    u = struct('thr', 1, 'CL', cfg.ac.CL_groundRoll, 'phi', 0);

    names = {'t','x','y','h','V','gam','psi','thr','CL','CLmax3D','CD','LD', ...
             'phi','n','T','D','Ib','Vt','cellV','SOC','P','Re','inRange', ...
             'rpm','dist','phase'};
    for i = 1:numel(names), Lg.(names{i}) = nan(N,1); end

    kEnd = N;
    for k = 1:N
        t = (k-1)*dt;
        h = s(3); V = s(4); gam = s(5); psi = s(6);
        Va = max(V, 0.5);
        q  = 0.5*rho*Va^2;
        Re = rho*Va*chord/cfg.env.mu;
        [~, CLmax2D] = polarLookup(M.polar, 0, Re, 0);
        CLmax3D = cfg.ac.CLmax3D_factor*CLmax2D;

        % ---------------- phase transitions ----------------
        if phase == 1
            if V >= Vrot
                phase = 2; onGround = false;
                groundRoll = s(9); tLiftoff = t;
            elseif s(9) > cfg.mis.runway_m
                status = 'FAIL: runway exceeded before Vrot';
                kEnd = k - 1; break;
            end
        end
        if phase == 2 && h >= hCmd - 1
            phase = 3; legMode = 1; psiRef = 0; legStart = s(9); iV = 1;
        end

        % ---------------- guidance & control ----------------
        u.phi = 0;
        if phase == 1
            u.thr = cfg.prop.thrMax;
            u.CL = cfg.ac.CL_groundRoll;
        else
            if phase == 2
                u.thr = cfg.prop.thrMax;
                gamCmd = min(ap.kh*(hCmd - h), gamMax) - ap.kVprot*max(Vcmd - V, 0);
                gamCmd = max(gamCmd, 0);
                u.phi  = clamp(ap.kpsi*wrapPi(psiRef - psi), -phiMaxS, phiMaxS);
            else
                e      = Vcmd - V;
                iV     = clamp(iV + ap.kV_i*e*dt, 0, 1);
                u.thr = min(clamp(ap.kV_p*e + iV, 0, 1), cfg.prop.thrMax);
                gamCmd = clamp(ap.kh*(hCmd - h), -gamMax, gamMax);

                if legMode == 1
                    u.phi = clamp(ap.kpsi*wrapPi(psiRef - psi), -phiMaxS, phiMaxS);
                    if s(9) - legStart >= cfg.mis.legLength_m
                        legMode = 2; turned = 0; psiPrev = psi;
                    end
                else
                    u.phi  = bank;
                    turned = turned + wrapPi(psi - psiPrev);
                    psiPrev = psi;
                    if turned >= pi
                        legMode = 1; psiRef = wrapPi(psiRef + pi); legStart = s(9);
                        turns = turns + 1;
                        if mod(turns, 2) == 0, laps = laps + 1; end
                    end
                end
            end
            % CL needed to track gamma command in a coordinated bank
            Lreq = (W*cos(gam) + mTot*Va*ap.kgam*(gamCmd - gam))/cos(u.phi);
            CLhi = min(cfg.ac.CL_clampFrac*CLmax3D, cfg.ac.nMax*W/(q*S));
            u.CL = clamp(Lreq/(q*S), cfg.ac.CL_floor, CLhi);
        end

        % ---------------- RK4 integration (ZOH controls) ----------------
        [k1, o] = plant(s,             u, cfg, M, mTot, onGround);
        k2      = plant(s + 0.5*dt*k1, u, cfg, M, mTot, onGround);
        k3      = plant(s + 0.5*dt*k2, u, cfg, M, mTot, onGround);
        k4      = plant(s + dt*k3,     u, cfg, M, mTot, onGround);

        Lg.t(k)=t; Lg.x(k)=s(1); Lg.y(k)=s(2); Lg.h(k)=h; Lg.V(k)=V;
        Lg.gam(k)=gam; Lg.psi(k)=psi; Lg.thr(k)=u.thr; Lg.CL(k)=u.CL;
        Lg.CLmax3D(k)=CLmax3D; Lg.CD(k)=o.CD; Lg.LD(k)=o.L/max(o.D,eps);
        Lg.phi(k)=u.phi; Lg.n(k)=o.L/W; Lg.T(k)=o.T; Lg.D(k)=o.D;
        Lg.Ib(k)=o.Ib; Lg.Vt(k)=o.Vt; Lg.cellV(k)=o.Vt/cfg.batt.cells;
        Lg.SOC(k)=o.SOC; Lg.P(k)=o.Ib*o.Vt; Lg.Re(k)=o.Re;
        Lg.inRange(k)=o.inRange; Lg.rpm(k)=o.rpm; Lg.dist(k)=s(9);
        Lg.phase(k)=phase;

        s = s + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        if onGround, s(3) = 0; s(5) = 0; s(4) = max(s(4), 0); end

        % ---------------- termination ----------------
        if ~onGround && s(3) < -0.2
            status = 'FAIL: ground impact'; kEnd = k; break;
        end
        if o.SOC <= cfg.batt.reserveSOC
            status = 'END: battery reserve SOC reached'; kEnd = k; break;
        end
        if o.Vt/cfg.batt.cells < cfg.batt.cutoffCell_V
            status = 'END: loaded cell voltage below cutoff'; kEnd = k; break;
        end
        if laps >= cfg.mis.nLaps
            status = 'COMPLETE: commanded laps flown'; kEnd = k; break;
        end
    end

    kEnd = max(kEnd, 1);
    for i = 1:numel(names), Lg.(names{i}) = Lg.(names{i})(1:kEnd); end

    air = Lg.phase >= 2;
    R.log = Lg;  R.status = status;  R.mTot = mTot;  R.W = W;
    R.Vs = Vs;   R.Vrot = Vrot;
    R.staticThrust_N = T0;  R.staticTW = T0/W;
    R.groundRoll_m = groundRoll;  R.tLiftoff = tLiftoff;
    R.laps = laps;  R.turns = turns;
    R.distance_m   = s(9);
    R.flightTime_s = Lg.t(end) - tLiftoff;
    R.energy_Wh    = s(8);  R.charge_Ah = s(7);
    R.maxIb = max(Lg.Ib);   R.maxC = R.maxIb/cfg.batt.capacity_Ah;
    R.minCellV = min(Lg.cellV);
    R.ReInRangeFrac    = mean(Lg.inRange(air));
    R.CLmarginExceed_s = sum(Lg.CL(air) > cfg.ac.CL_margin*Lg.CLmax3D(air))*dt;
    R.CLclamp_s        = sum(Lg.CL(air) >= 0.999*cfg.ac.CL_clampFrac*Lg.CLmax3D(air))*dt;
    R.pass = ~startsWith(status, 'FAIL') && laps >= cfg.mis.requiredLaps && ...
             R.maxC <= cfg.batt.maxC;
    R.Vcmd = Vcmd;
end

function cfg = deriveWing(cfg)
    cfg.ac.S  = cfg.ac.span_m*cfg.ac.chord_m;   % m^2
    cfg.ac.AR = cfg.ac.span_m/cfg.ac.chord_m;   % rectangular wing
end
 
function CS = chordSweep(cfg, M)
    ft = 3.280839895;
    cfg.mis.dt    = cfg.sweep.dt;
    cfg.mis.nLaps = cfg.sweep.chordNLaps;
    C = cfg.sweep.chords(:);  n = numel(C);
    S = nan(n,1); AR = nan(n,1); Vs = nan(n,1); Vrot = nan(n,1);
    laps = zeros(n,1); dist_ft = nan(n,1); t_s = nan(n,1);
    roll_ft = nan(n,1); peakC = nan(n,1); mAh = nan(n,1);
    clFrac = nan(n,1); ok = false(n,1); status = strings(n,1);
 
    fprintf('\n============ CHORD SWEEP (span %.3f m, payload %.2f kg) ============\n', ...
        cfg.ac.span_m, cfg.mass.payload_kg);
    for i = 1:n
        c = cfg; 
        c.ac.chord_m = C(i); 
        c = deriveWing(c);
        r  = runMission(c, M);
        ph = r.log.phase == 3;
        S(i) = c.ac.S;  AR(i) = c.ac.AR;  Vs(i) = r.Vs;  Vrot(i) = r.Vrot;
        laps(i) = r.laps;  
        dist_ft(i) = r.distance_m*ft; 
        if startsWith(r.status, 'FAIL'), dist_ft(i) = NaN; end
        t_s(i) = r.flightTime_s;
        roll_ft(i) = r.groundRoll_m*ft;  peakC(i) = r.maxC;  mAh(i) = r.charge_Ah*1000;
        clFrac(i) = max([0; r.log.CL(ph)./r.log.CLmax3D(ph)]);   % worst CL/CLmax in the laps
        % acceptable = mission passes AND stays inside your own CL margin AND Vrot below cruise speed
        ok(i) = r.pass && clFrac(i) <= cfg.ac.CL_margin && r.Vrot <= r.Vcmd;
        status(i) = string(r.status);
        fprintf(['%5.1f cm | S %.3f | AR %5.1f | Vs %5.2f Vrot %5.2f | laps %3d | %7.0f ft | ' ...
                 'roll %5.1f ft | %4.1fC | CL/CLmax %.2f | %s | %s\n'], ...
            100*C(i), S(i), AR(i), Vs(i), Vrot(i), laps(i), dist_ft(i), ...
            roll_ft(i), peakC(i), clFrac(i), r.status, ...
            string(ifelse(ok(i), 'OK', 'check')));
    end
    CS = table(C, S, AR, Vs, Vrot, laps, dist_ft, t_s, roll_ft, peakC, mAh, clFrac, ok, status, ...
        'VariableNames', {'Chord_m','S_m2','AR','Vstall','Vrot','Laps','Distance_ft', ...
                          'FlightTime_s','GroundRoll_ft','PeakC','Charge_mAh','CLoverCLmax','Acceptable','Status'});
 
    figure('Name','Chord sweep','Color','w');
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile; plot(100*C, dist_ft, '-o'); hold on; grid on;
    plot(100*C(ok), dist_ft(ok), 'o', 'MarkerFaceColor', [0.2 0.7 0.2]);
    xlabel('Chord (cm)'); ylabel('Distance per charge (ft)');
    title('Range vs chord (filled = acceptable)');
    nexttile; plot(100*C, clFrac, '-o'); hold on; grid on;
    yline(cfg.ac.CL_margin, '--r', 'CL margin');
    xlabel('Chord (cm)'); ylabel('max C_L / C_{Lmax} during laps');
    title('Stall margin vs chord');
end
 
function s = ifelse(cond, a, b)
    if cond, s = a; else, s = b; end
end

%% ============================== PLANT ====================================
function [ds, o] = plant(s, u, cfg, M, mTot, onGround)
    g = cfg.env.g; rho = cfg.env.rho;
    h = s(3); V = s(4); gam = s(5); psi = s(6);
    Va  = max(V, 0.5);
    SOC = cfg.batt.SOC0 - s(7)/cfg.batt.capacity_Ah;
    W   = mTot*g;
    S   = cfg.ac.S; AR = cfg.ac.AR; b = sqrt(AR*S); c = S/b;
    q   = 0.5*rho*Va^2;
    Re  = rho*Va*c/cfg.env.mu;

    [CD2D, ~, ~, inRange] = polarLookup(M.polar, u.CL, Re, cfg.ac.ReCDexp);
    r   = (16*(max(h,0) + cfg.ac.wingHeight_m)/b)^2;
    kGE = r/(1 + r);                                   % ground-effect factor
    CD  = CD2D + cfg.ac.CD0_airframe + kGE*u.CL^2/(pi*cfg.ac.e*AR);
    L   = q*S*u.CL;
    D = q*(S*CD + M.fCdA(Va));

    [T, Ib, Vt, rpm] = propulsion(u.thr, Va, SOC, cfg, M);

    if onGround
        Nf   = max(W - L, 0);
        dV   = (T - D - cfg.ac.mu_roll*Nf)/mTot;
        dgam = 0; dpsi = 0; gE = 0;
    else
        dV   = (T - D - W*sin(gam))/mTot;
        dgam = (L*cos(u.phi) - W*cos(gam))/(mTot*Va);
        dpsi = L*sin(u.phi)/(mTot*Va*cos(gam));
        gE   = gam;
    end
    vx = V*cos(gE)*cos(psi) + cfg.env.wind(1);
    vy = V*cos(gE)*sin(psi) + cfg.env.wind(2);

    ds = [vx; vy; V*sin(gE); dV; dgam; dpsi; Ib/3600; Ib*Vt/3600; hypot(vx, vy)];

    if nargout > 1
        o = struct('L',L,'D',D,'T',T,'Ib',Ib,'Vt',Vt,'CD',CD,'Re',Re, ...
                   'inRange',inRange,'rpm',rpm,'SOC',SOC);
    end
end

%% ===================== PROPULSION + BATTERY ==============================
function [T, Ib, Vt, rpm] = propulsion(thr, V, SOC, cfg, M)
    ocv = M.ocvPack(min(max(SOC, 0), 1));
    Rp  = cfg.batt.cells*cfg.batt.R_cell + cfg.batt.R_wiring;
    Ib  = 0; T = 0; rpm = 0;
    for it = 1:4                          % battery <-> motor coupling
        Vt = max(ocv - Ib*Rp, 0);
        if strcmp(M.propMode, 'model')
            [T, Im, rpm] = motorProp(thr*Vt, V, cfg);
            Ib = thr*Im/cfg.prop.etaESC;
        else
            sc  = Vt/cfg.prop.tableVtest;  % first-order voltage correction
            T   = max(M.FT(V, thr), 0)*sc^2;
            Ib  = max(M.FI(V, thr), 0)*sc^2;
            rpm = NaN;
        end
    end
    Vt = max(ocv - Ib*Rp, 0);
end

function [T, Im, rpm] = motorProp(Vm, V, cfg)
% Torque balance  Qmotor(n) = Qprop(n), solved in closed form.
%   Qprop = rho*n^2*D^5*CP(J)/(2*pi),  J = V/(n*D)  -> quadratic in n
%   Qmotor = Kt*((Vm - Kt*2*pi*n)/Rm - I0)          -> linear in n
    p = cfg.prop; rho = cfg.env.rho; D = p.D_m;
    Kt = 60/(2*pi*p.Kv);
    A0 = Kt*(Vm/p.Rm - p.I0);
    A1 = 2*pi*Kt^2/p.Rm;
    kk = rho/(2*pi);
    a  = kk*p.CP(1)*D^5;
    b  = kk*p.CP(2)*D^4*V + A1;
    c  = kk*p.CP(3)*D^3*V^2 - A0;
    disc = b^2 - 4*a*c;
    if A0 <= 0 || disc < 0
        T = 0; Im = 0; rpm = 0; return;   % motor off / freewheel
    end
    n = (-b + sqrt(disc))/(2*a);          % rev/s
    if n <= 0
        T = 0; Im = 0; rpm = 0; return;
    end
    T   = rho*(p.CT(1)*n^2*D^4 + p.CT(2)*V*n*D^3 + p.CT(3)*V^2*D^2);
    Im  = max((Vm - Kt*2*pi*n)/p.Rm, 0);
    rpm = 60*n;
end

%% ============================ AIRFOIL POLAR ==============================
function P = loadPolar(file)
    if ~isfile(file), error('Polar file not found: %s', file); end
    T = readtable(file);
    need = {'Re','alpha_deg','CL','CD'};
    assert(all(ismember(need, T.Properties.VariableNames)), ...
        'Polar %s needs columns Re, alpha_deg, CL, CD.', file);
    T = rmmissing(T, 'DataVariables', need);
    ReKey  = round(T.Re, -2);              % absorb digitizing noise in Re
    ReVals = unique(ReKey);

    P.blocks = struct('Re',{},'fCD',{},'fA',{},'CLmax',{},'CLmin',{});
    for k = 1:numel(ReVals)
        Rk = sortrows(T(ReKey == ReVals(k), :), 'alpha_deg');
        [~, iMax] = max(Rk.CL);
        Rk = Rk(1:iMax, :);                 % pre-stall branch only
        cm   = cummax(Rk.CL);               % drop digitizer wiggles
        keep = [true; Rk.CL(2:end) > cm(1:end-1)];
        Rk   = Rk(keep, :);
        if height(Rk) < 4
            warning('Re = %g block in %s has < 4 usable points; skipped.', ReVals(k), file);
            continue;
        end
        P.blocks(end+1) = struct( ...
            'Re',    ReVals(k), ...
            'fCD',   griddedInterpolant(Rk.CL, Rk.CD, 'linear', 'nearest'), ...
            'fA',    griddedInterpolant(Rk.CL, Rk.alpha_deg, 'linear', 'nearest'), ...
            'CLmax', Rk.CL(end), ...
            'CLmin', Rk.CL(1));
    end
    assert(~isempty(P.blocks), 'No usable Re blocks in %s.', file);
    P.Re   = [P.blocks.Re];
    P.CLmin = max([P.blocks.CLmin]);
    P.file = file;
end

function [CD2D, CLmax2D, alpha, inRange] = polarLookup(P, CL, Re, ReExp)
    ReB = P.Re; nB = numel(ReB);
    if nB == 1 || Re <= ReB(1) || Re >= ReB(end)
        [~, k] = min(abs(log(ReB/Re)));
        bk = P.blocks(k);
        c  = min(max(CL, bk.CLmin), bk.CLmax);
        CD2D    = bk.fCD(c)*(Re/ReB(k))^ReExp;
        CLmax2D = bk.CLmax;
        alpha   = bk.fA(c);
        inRange = abs(log(Re/ReB(k))) < log(1.10);
    else
        k2 = find(ReB >= Re, 1); k1 = k2 - 1;
        w  = log(Re/ReB(k1))/log(ReB(k2)/ReB(k1));
        b1 = P.blocks(k1); b2 = P.blocks(k2);
        c1 = min(max(CL, b1.CLmin), b1.CLmax);
        c2 = min(max(CL, b2.CLmin), b2.CLmax);
        CD2D    = (1-w)*b1.fCD(c1) + w*b2.fCD(c2);
        CLmax2D = (1-w)*b1.CLmax   + w*b2.CLmax;
        alpha   = (1-w)*b1.fA(c1)  + w*b2.fA(c2);
        inRange = true;
    end
end

%% ============================== SETUP ====================================
function M = buildModels(cfg)
    M.polar      = loadPolar(cfg.ac.polarFile);

    CdA = cfg.ac.dragArt.F ./ (0.5*cfg.env.rho*cfg.ac.dragArt.V.^2);
    M.fCdA = griddedInterpolant(cfg.ac.dragArt.V, CdA, 'linear', 'nearest');

    M.CLmax2Dref = max([M.polar.blocks.CLmax]);
    M.ocvPack    = griddedInterpolant(cfg.batt.ocvSOC, ...
                     cfg.batt.cells*cfg.batt.ocvCell, 'linear', 'nearest');
    M.propMode   = lower(cfg.prop.mode);
    switch M.propMode
        case 'model'
        case 'table'
            T = readtable(cfg.prop.tableFile);
            need = {'V_ms','throttle','thrust_N','current_A'};
            assert(all(ismember(need, T.Properties.VariableNames)), ...
                'Thrust table needs columns V_ms, throttle, thrust_N, current_A.');
            M.FT = scatteredInterpolant(T.V_ms, T.throttle, T.thrust_N,  'linear', 'nearest');
            M.FI = scatteredInterpolant(T.V_ms, T.throttle, T.current_A, 'linear', 'nearest');
        otherwise
            error('cfg.prop.mode must be ''model'' or ''table''.');
    end
end

function cfg = applyOptimizerResult(cfg)
    names = {'af300','aquila','e387','fx63_187','s1210','s1211','s1223','sg6043','sh3055'};
    files = {'polars/af300_Re100000_digitized.csv', 'polars/aquila_Re100000_digitized.csv', ...
             'polars/e387_Re100000_digitized.csv',  'polars/fx63_137_Re100000_digitized.csv', ...
             'polars/s1210_Re100000_digitized.csv', 'polars/s1221_Re100000_digitized.csv', ...
             'polars/s1223_Re100000_digitized.csv', 'polars/sg6043_Re100000_digitized.csv', ...
             'polars/sh3055_Re100000_digitized.csv'};
    map = containers.Map(names, files);
    T = readtable(cfg.optimizerFile, 'TextType', 'string');
    r = T(1,:);
    cfg.ac.airfoil      = char(r.Airfoil);
    cfg.ac.polarFile    = map(cfg.ac.airfoil);
    cfg.ac.S            = r.WingArea_m2;
    cfg.ac.AR           = r.AR;
    cfg.mass.payload_kg = r.InternalPayload_kg;
    fprintf('Loaded optimizer design: %s, S = %.3f, AR = %.1f, payload = %.2f kg (opt. V = %.1f m/s)\n', ...
        cfg.ac.airfoil, cfg.ac.S, cfg.ac.AR, cfg.mass.payload_kg, r.Speed_ms);
end

%% ============================== OUTPUT ===================================
function printSummary(R, cfg)
    ft = 3.280839895;
    yn = {'NO','YES'};
    fprintf('\n================ MISSION SUMMARY ================\n');
    fprintf('Status:               %s\n', R.status);
    fprintf('Airfoil / S / AR:     %s / %.3f m^2 / %.1f\n', cfg.ac.airfoil, cfg.ac.S, cfg.ac.AR);
    fprintf('Total mass:           %.3f kg (payload %.3f kg)\n', R.mTot, cfg.mass.payload_kg);
    fprintf('Vstall(1g) / Vrot:    %.2f / %.2f m/s\n', R.Vs, R.Vrot);
    fprintf('Static thrust, T/W:   %.2f N, %.2f (target %.2f)\n', ...
        R.staticThrust_N, R.staticTW, cfg.mis.targetTW);
    fprintf('Ground roll:          %.1f ft (runway %.1f ft)\n', ...
        R.groundRoll_m*ft, cfg.mis.runway_m*ft);
    fprintf('Laps completed:       %d\n', R.laps);
    fprintf('Distance flown:       %.0f ft\n', R.distance_m*ft);
    fprintf('Flight time:          %.1f s\n', R.flightTime_s);
    fprintf('Energy used:          %.2f Wh (%.0f mAh)\n', R.energy_Wh, R.charge_Ah*1000);
    fprintf('Peak current:         %.1f A (%.1fC, rating %dC)\n', R.maxIb, R.maxC, cfg.batt.maxC);
    fprintf('Min loaded cell V:    %.2f V\n', R.minCellV);
    fprintf('Time CL > %.0f%% CLmax: %.1f s | time at CL clamp: %.1f s\n', ...
        100*cfg.ac.CL_margin, R.CLmarginExceed_s, R.CLclamp_s);
    fprintf('Airborne time with Re within +/-10%% of polar data: %.0f%%\n', 100*R.ReInRangeFrac);
    fprintf('Mission PASS (>= %d laps, no FAIL, current OK): %s\n', ...
        cfg.mis.requiredLaps, yn{R.pass + 1});
end

function plotMission(R, cfg, M)
    L = R.log; ft = 3.280839895; air = L.phase >= 2;

    figure('Name','Flight','Color','w');
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile; plot(L.x*ft, L.y*ft, 'LineWidth', 1.2); axis equal; grid on;
    xlabel('x (ft)'); ylabel('y (ft)'); title('Ground track');
    nexttile; plot(L.t, L.h*ft, 'LineWidth', 1.2); grid on;
    xlabel('t (s)'); ylabel('Altitude (ft)'); title('Altitude');
    nexttile; plot(L.t, L.V, 'LineWidth', 1.2); hold on; grid on;
    yline(R.Vs, '--', 'V_{stall}'); 
    yline(R.Vcmd, ':', 'V_{cmd}');
    xlabel('t (s)'); ylabel('Airspeed (m/s)'); title('Airspeed');
    nexttile; plot(L.t, L.CL, 'LineWidth', 1.2); hold on; grid on;
    plot(L.t, cfg.ac.CL_margin*L.CLmax3D, '--');
    plot(L.t, L.CLmax3D, ':');
    legend('CL', 'margin x CLmax', 'CLmax (3D)', 'Location', 'best');
    xlabel('t (s)'); ylabel('C_L'); title('Lift coefficient (turns show as steps)');

    figure('Name','Battery & propulsion','Color','w');
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile; plot(L.t, L.Ib); grid on;
    yline(cfg.batt.maxC*cfg.batt.capacity_Ah, '--r', 'C rating');
    xlabel('t (s)'); ylabel('Battery current (A)');
    nexttile; plot(L.t, L.cellV); grid on;
    yline(cfg.batt.cutoffCell_V, '--r', 'cutoff');
    xlabel('t (s)'); ylabel('Loaded cell voltage (V)');
    nexttile; plot(L.t, 100*L.SOC); grid on;
    yline(100*cfg.batt.reserveSOC, '--r', 'reserve');
    xlabel('t (s)'); ylabel('SOC (%)');
    nexttile; yyaxis left; plot(L.t, 100*L.thr); ylabel('Throttle (%)');
    yyaxis right; plot(L.t, L.T); ylabel('Thrust (N)'); grid on; xlabel('t (s)');

    figure('Name','Power vs distance','Color','w'); hold on; grid on;
    phNames = {'Takeoff roll','Climb','Laps'};
    for p = 1:3
        m = L.phase == p;
        if any(m)
            plot(L.dist(m)*ft, L.P(m), '.', 'MarkerSize', 4, 'DisplayName', phNames{p});
        end
    end
    xlabel('Distance flown (ft)'); ylabel('Electrical power (W)');
    title(sprintf('3S %.0f mAh: power vs distance | payload %.2f kg | %s', ...
        cfg.batt.capacity_Ah*1000, cfg.mass.payload_kg, cfg.ac.airfoil), 'Interpreter', 'none');
    legend('show', 'Location', 'best');

    figure('Name','Reynolds number','Color','w');
    plot(L.t(air), L.Re(air), 'LineWidth', 1.2); hold on; grid on;
    for r = M.polar.Re, yline(r, '--k', sprintf('polar Re %g', r)); end
    xlabel('t (s)'); ylabel('Re');
    title(sprintf('Re in flight: %.0f%% of airborne time within +/-10%% of polar data', ...
        100*R.ReInRangeFrac));

    figure('Name','Power vs Time','Color','w');
    plot(L.t, L.P, 'LineWidth', 1.2); 
    hold on; grid on;
    xlabel('Time (s)'); 
    ylabel('Electrical Power (W)');
    title(sprintf('Wattage vs Time | Total Consumed: %.0f mAh', R.charge_Ah*1000));
end

function SW = payloadSweep(cfg, M)
    cfg.sweep.chordEnable = true;
    cfg.sweep.chords      = 0.12:0.01:0.28;   % m
    cfg.sweep.chordNLaps  = 2;              % inf = fly to battery reserve (range)
    ft = 3.280839895;
    cfg.mis.dt = cfg.sweep.dt;
    P = cfg.sweep.payloads(:); n = numel(P);
    pass = false(n,1); laps = zeros(n,1); 
    dist_ft = nan(n,1);
    roll_ft = nan(n,1); maxC = nan(n,1); status = strings(n,1);

    passTxt = {'fail','PASS'};
    fprintf('\n================ PAYLOAD SWEEP ================\n');
    for i = 1:n
        c = cfg; c.mass.payload_kg = P(i);
        r = runMission(c, M);
        pass(i) = r.pass; laps(i) = r.laps; 
        dist_ft(i) = r.distance_m*ft;
        if startsWith(r.status, 'FAIL'), dist_ft(i) = NaN; end
        roll_ft(i) = r.groundRoll_m*ft; maxC(i) = r.maxC; status(i) = string(r.status);
        fprintf('%5.2f kg | %-40s | laps %3d | %7.0f ft | roll %5.1f ft | %4.1fC | %s\n', ...
            P(i), r.status, r.laps, dist_ft(i), roll_ft(i), r.maxC, ...
            passTxt{r.pass + 1});
    end
    SW = table(P, pass, laps, dist_ft, roll_ft, maxC, status, 'VariableNames', ...
        {'Payload_kg','Pass','Laps','Distance_ft','GroundRoll_ft','PeakC','Status'});

    if any(pass)
        fprintf('\nMax payload passing the full mission: %.2f kg\n', max(P(pass)));
    else
        fprintf('\nNo payload in the sweep passed the mission.\n');
    end

    figure('Name','Payload sweep','Color','w');
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile; plot(P, dist_ft, '-o'); hold on; grid on;
    plot(P(pass), dist_ft(pass), 'o', 'MarkerFaceColor', [0.2 0.7 0.2]);
    xlabel('Internal payload (kg)'); ylabel('Distance per charge (ft)');
    title('Range vs payload (filled = pass)');
    nexttile; plot(P, roll_ft, '-o'); hold on; grid on;
    yline(cfg.mis.runway_m*ft, '--r', 'Runway');
    xlabel('Internal payload (kg)'); ylabel('Ground roll (ft)');
    title('Takeoff ground roll vs payload');
end

function selfCheck(R, cfg)
    L = R.log;
    fprintf('\n================ SELF-CHECK ================\n');
 
    % 1) energy bookkeeping: integrated log vs the integrated state variables
    Wh_int = trapz(L.t, L.P)/3600;   Ah_int = trapz(L.t, L.Ib)/3600;
    fprintf('Energy: log %.3f Wh vs state %.3f Wh | charge: log %.4f Ah vs state %.4f Ah\n', ...
        Wh_int, R.energy_Wh, Ah_int, R.charge_Ah);
 
    % 2) steady straight cruise: thrust = drag, load factor = 1
    m = L.phase == 3 & abs(L.phi) < deg2rad(3) & abs(L.V - R.Vcmd) < 0.2;
    if any(m)
        fprintf('Cruise: median T %.2f N vs D %.2f N | median n %.3f (expect 1.000)\n', ...
            median(L.T(m)), median(L.D(m)), median(L.n(m)));
    else
        fprintf('Cruise: no steady straight segment found (check the speed plot)\n');
    end
 
    % 3) steady banked turns: n should be 1/cos(bank)
    tn = L.phase == 3 & abs(L.phi - deg2rad(cfg.mis.bank_deg)) < 1e-6;
    if any(tn)
        fprintf('Turn:   median n %.3f vs 1/cos(bank) = %.3f\n', ...
            median(L.n(tn)), 1/cos(deg2rad(cfg.mis.bank_deg)));
    end
 
    % 4) altitude hold during the laps
    a3 = L.phase == 3;
    if any(a3)
        err = max(abs(L.h(a3) - cfg.mis.cruiseAlt_m));
        fprintf('Max altitude error in laps: %.1f m %s\n', err, ...
            string(ifelse(err > 5, '  <-- WARNING: not holding altitude', '')));
    end
 
    % 5) lift coefficient pinned at the clamp = flying at the edge of stall
    if R.CLclamp_s > 0.5
        fprintf('WARNING: CL pinned at the clamp for %.1f s - results are not trustworthy\n', R.CLclamp_s);
    end
    if startsWith(R.status, 'FAIL')
        fprintf('WARNING: mission status is %s - distance/energy are for a crashed run\n', R.status);
    end
end
 
function convergenceCheck(cfg, M)
    dts = [0.04 0.02 0.01];
    fprintf('\n============ TIME-STEP CONVERGENCE ============\n');
    for i = 1:numel(dts)
        c = cfg;  c.mis.dt = dts(i);
        r = runMission(c, M);
        fprintf('dt %.3f s: %7.0f ft | %6.1f mAh | roll %5.1f ft | %s\n', ...
            dts(i), r.distance_m*3.280839895, r.charge_Ah*1000, r.groundRoll_m*3.280839895, r.status);
    end
    fprintf('Differences below about 1%% between 0.02 and 0.01 mean dt is fine.\n');
end

%% ============================== HELPERS ==================================
function y = clamp(x, lo, hi), y = min(max(x, lo), hi); end
function a = wrapPi(a), a = mod(a + pi, 2*pi) - pi; end