
CARGO AIRCRAFT / AIRFOIL MATLAB MODEL
======================================

Files
-----
cargo_airfoil_optimizer.m
polars/        <- put airfoil polar CSV files here

Goal
----
1. Maximize internal payload for a cargo aircraft.
2. Screen airfoils using actual CL-CD polar data at the simulated Reynolds number.
3. Include:
   - empty mass = 1.20 kg
   - fixed drag-article mass = 1.20 kg
   - battery mass = 0.150 kg
   - battery = 3S 1000 mAh = 11.1 Wh nominal
   - usable battery energy = 80% of nominal by default
   - target T/W = 0.50
4. Produce:
   - airfoil/design ranking
   - payload vs power plot
   - power (W) vs battery range (ft) plot

Polar CSV format
----------------
One CSV per airfoil.

Required columns:
Re, alpha_deg, CL, CD

Optional:
CM

Example:
Re,alpha_deg,CL,CD,CM
150000,-4,-0.20,0.030,-0.02
150000,-3,-0.10,0.024,-0.02
...

Use real XFOIL or wind-tunnel polar data. Do not create fake polar data for design decisions.

Useful source
-------------
UIUC Airfoil Data Site:
https://m-selig.ae.illinois.edu/ads/coord_database.html

UIUC FAQ / low-speed test-data information:
https://m-selig.ae.illinois.edu/ads_faq.html

Important assumptions to replace
--------------------------------
- wing-area search range
- aspect-ratio range
- max span
- chord bounds
- stall-speed target
- air density / altitude
- powertrain efficiency
- T/W target if different from 0.50
- battery usable-energy fraction
- structural load factor
- cargo volume/CG constraints
- measured motor/prop thrust curve

The current T/W model is a first-pass constraint:
    T_available = 0.5 * W

For final propulsion sizing, replace this with measured or manufacturer
thrust-vs-speed data.

Power/range relation
--------------------
Electrical cruise power:
    P_batt = D*V / eta_powertrain

Battery range:
    R = V*(E_usable_J/P_batt)

The graph therefore shows the tradeoff between cruise electrical power and
distance available from the fixed 3S 1000 mAh battery at each speed.
