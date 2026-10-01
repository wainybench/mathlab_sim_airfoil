AIRFOIL POLARS DIGITIZED FROM THE UPLOADED PAPER

Source:
  Singh et al., "Design of a low Reynolds number airfoil for small
  horizontal axis wind turbines", Renewable Energy 42 (2012), 66-76.

The paper compares AF300 with s1223, s1210, s1221, SH3055, FX 63-137,
E387, SG6043, and Aquila at Re = 100,000 using XFOIL/data. The paper
presents the CL-alpha, L/D-alpha, and drag-polar information as figures
plus Table 3 summary values. It does not provide these nine complete
polar datasets as downloadable CSV tables in the supplied PDF.

Therefore:
  - *_Re100000_digitized.csv are APPROXIMATE digitizations of Figures 3
    and 4 at Re = 100,000, sampled every 2 degrees from alpha=0 to 20.
  - CD was reconstructed approximately using CD = CL/(L/D).
  - CLmax values at/near the reported stall angles were pinned to the
    exact Table 3 values where available.
  - paper_table3_Re100000_summary.csv contains the exact summary values
    transcribed from Table 3.

DO NOT treat the digitized CSVs as raw experimental/XFOIL data. Use them
for screening/initial MATLAB work. For final airfoil selection, replace
them with original XFOIL/wind-tunnel polar files at the actual Reynolds
numbers of the aircraft.

Candidate airfoils:
  AF300, s1210, s1223, s1221, SH3055, FX63-137, E387, SG6043, Aquila

Note on SH3055:
  The paper states that SH3055 coordinates could not be obtained, so
  available lift/drag data were used for that airfoil.
