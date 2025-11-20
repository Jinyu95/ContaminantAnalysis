'''
Run contaminant calculation, ETACHA4 charge state distribution,
and beam tracking through post-stripper lattice.
!!! Before running this script, ensure that ETACHA4 and pyJuTrack are properly installed.
!!! Read README.md for setup instructions.
'''

from contaminant_etacha_integration import calculate_and_run_contaminants
import os
from beam_tracking import track_contaminants

## Calculate contaminants and run ETACHA4 to get charge state distributions
lithium_target = {
      'At': 6.94,       # Lithium
      'Zt': 3,          # Z = 3 (Lithium)
      'Thick': 0.5,     # mg/cm²
      'd': 0.51         # density g/cm³
}

prefix = "Zn64_analysis"

# This creates: `prefix_list.csv`,  `prefix_charge_distributions.csv`.
# The first file contains the calculated contaminants based on A/q difference.
# The second file contains the charge state distributions calculated by ETACHA4.
contaminants, results = calculate_and_run_contaminants(
    target_A=64,         # Main beam mass number at dipole
    target_q=19,         # Main beam charge state at dipole
    beam_energy=20.0,    # Energy in MeV/u at stripper
    target_material=lithium_target,
    output_prefix=prefix
)


## Beam tracking of main beam and contaminants through post-stripper lattice
# main beam at stripper
main_A = 64  
main_q = 28.5
main_energy = 19.85  # MeV/u


# import the lattice and aperture.
import stripper2CSS_lattice
lattice = stripper2CSS_lattice.create_stripper2CSS_lattice()

print(f"  Lattice length: {lattice.total_length():.2f} m")

# Beam parameters
emit_n_x = 0.15 * 1e-6  # m·rad (normalized emittance)
emit_n_y = 0.15 * 1e-6  # m·rad (normalized emittance)
beta_twiss_x = 0.2  # m
beta_twiss_y = 0.2  # m
alpha_twiss_x = 0.0  
alpha_twiss_y = 0.0  


csv_file = f"{prefix}_charge_distributions.csv" 

# make dir
output_dir = f"{prefix}_tracking_results"
os.makedirs(output_dir, exist_ok=True)
output_prefix = os.path.join(output_dir, "tracking")

track_contaminants(csv_file, 
                  lattice=lattice,
                  aperture=None,  # Not used. apertures are defined in lattice elements
                  main_A=main_A,
                  main_q=main_q,
                  emit_x=emit_n_x,
                  emit_y=emit_n_y,
                  beta_twiss_x=beta_twiss_x,
                  beta_twiss_y=beta_twiss_y,
                  alpha_twiss_x=alpha_twiss_x,
                  alpha_twiss_y=alpha_twiss_y,
                  output_prefix=output_prefix,
                  x_offset=0.0, 
                  y_offset=0.0,
                  s_plot_range_m=16.9,
                  plot_ylim_mm=80.0, 
                  main_energy_MeV_u=main_energy)  
