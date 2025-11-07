# Contaminant Calculator + ETACHA4 Integration
This code calculate possible ion contaminants after dipole based on the A/q difference. 
The charge state distribution of all contaminants after stripping is then calculated with ETACH4.

## Setup (Required First!)

1. Open `etacha_batch_run.py` and update these paths:
   - `ETACHA_EXE` - Path to your ETACHA4.exe
   - `RESULTS_DIR` - Where ETACHA saves results. Usually in Documents/LISEcute/results
2. Copy `template.etacha` to this directory
3. Contaminant species range is defined as DEFAULT_SPECIES in contaminant_calculator.py.
Change the definition of DEFAULT_SPECIES to search a customized species range.

## Usage

```python
from contaminant_etacha_integration import calculate_and_run_contaminants

# carbon stripper
carbon_target = {'At': 12, 'Zt': 6, 'Thick': 1.0, 'd': 2.26}

contaminants, results = calculate_and_run_contaminants(
    target_A=208,        # Main beam mass number at dipole
    target_q=64,         # Main beam charge state at dipole
    beam_energy=30.0,    # Energy in MeV/u at stripper
    target_material=carbon_target,
    output_prefix="my_analysis"
)
```

This creates: `my_analysis_list.csv`,  `my_analysis_charge_distributions.csv`.
The first one is contaminant list and the second one is charge state distribution of the contanminants after stripping.

To plot the results, run 
```python
from visualize_charge_distribution import plot_charge_distribution
contaminants_file = "my_analysis_list.csv"
charge_dist_file = "my_analysis_charge_distributions.csv"
plot_charge_distribution(contaminants_file, charge_dist_file, "charge_states.png")
```