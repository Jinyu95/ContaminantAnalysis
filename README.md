# Contaminant Analysis for heavy ion accelerators
This code calculates possible ion contaminants based on the A/q difference. 
The charge state distribution of all contaminants after stripping is calculated with ETACHA4.
Beam tracking and beam loss analysis are then implemented with JuTrack (pyJuTrack).

## Setup (Required First!)

1. Open `etacha_batch_run.py` and update these paths:
   - `ETACHA_EXE` - Path to your ETACHA4.exe
   - `RESULTS_DIR` - Where ETACHA saves results. Usually in Documents/LISEcute/results
2. Copy `template.etacha` to this directory
3. Contaminant species range is defined as DEFAULT_SPECIES in contaminant_calculator.py.
   Change the definition of DEFAULT_SPECIES to search a customized species range.
4. For beam tracking, install pyJuTrack (Python version 3.13 or later is required!):
   ```bash
   git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl.git
   cd JuTrack.jl/python_integration
   pip install -e .
   ```

## Usage

See `run_analysis.py`. The first run of JuTrack will install all dependencies that takes a few minutes.
Don't forget to clean up the results folder `RESULTS_DIR` after running the code.

## Visualize Charge Distributions before and after stripper

```python
from visualize_charge_distribution import plot_charge_distribution
contaminants_file = "my_analysis_list.csv"
charge_dist_file = "my_analysis_charge_distributions.csv"
plot_charge_distribution(contaminants_file, charge_dist_file, "charge_states.png")
```
