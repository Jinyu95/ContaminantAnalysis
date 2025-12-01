"""
Integration between Contaminant Calculator and ETACHA Batch Runner
Automatically creates ETACHA input files for calculated contaminants and processes results
Created by Jinyu Wan [wan@frib.msu.edu]
"""

import csv
import subprocess
from pathlib import Path
from contaminant_calculator import calculate_contaminants, save_results_to_csv
from etacha_batch_run import (
    create_case_files, run_case, compute_moments, 
    parse_forExcel, parse_eta, ETACHA_EXE, RESULTS_DIR
)


def create_betacha_from_contaminants(results, beam_energy, target_material, 
                                     output_file="contaminants.betacha"):
    """
    Create a .betacha input file from contaminant calculation results
    
    Parameters:
    -----------
    results : dict
        Results dictionary from calculate_contaminants()
    beam_energy : float or list
        Beam energy in MeV/u. Can be a single value or list of energies
    target_material : dict
        Target material parameters: {'At': int, 'Zt': int, 'Thick': float, 'd': float}
        At = target mass number
        Zt = target atomic number  
        Thick = thickness (mg/cm²)
        d = density (g/cm³)
    output_file : str
        Output .betacha filename
    
    Returns:
    --------
    Path : Path to created .betacha file
    """
    contaminants = results['contaminants']
    
    if isinstance(beam_energy, (int, float)):
        energies = [beam_energy]
    else:
        energies = beam_energy
    
    betacha_path = Path(output_file)
    
    with open(betacha_path, 'w') as f:
        f.write("! Ap, Zp, q,  E, At, Zt, Thick, d\n")
        f.write("! Auto-generated from contaminant calculator\n")
        f.write(f"! Target beam: A={results['target_A']}, q={results['target_q']}, "
                f"A/q={results['target_Aq']:.4f}\n")
        f.write(f"! Tolerance: ±{results['tolerance_percent']}%\n")
        f.write(f"! Found {results['num_contaminants']} contaminants\n")
        f.write("!\n")
        
        target_A = results['target_A']
        target_q = results['target_q']
        
        # Try to find the target species from the contaminants list
        target_species = None
        for cont in contaminants:
            if cont['mass_number'] == target_A and cont['charge'] == target_q:
                target_species = cont['species']
                break
        
        # If not found in contaminants, try to find a matching isotope in species_database
        if not target_species:
            from contaminant_calculator import species_database
            for element, data in species_database.items():
                for mass_amu, abundance in data['isotopes'].items():
                    if int(round(mass_amu)) == target_A:
                        target_species = element
                        break
                if target_species:
                    break
        
        # Get atomic number
        target_Z = get_atomic_number(target_species) if target_species else 0
        if target_Z == 0:
            # Fallback: estimate Z from A (this may be inaccurate)
            target_Z = int(target_A / 2.5)
            print(f"Warning: Could not determine Z for A={target_A}. Using estimate: Z={target_Z}")

        
        f.write(f"! Main beam:\n")
        for E in energies:
            f.write(f"{target_A} {target_Z} {target_q}  {E} {target_material['At']} "
                   f"{target_material['Zt']} {target_material['Thick']} "
                   f"{target_material['d']}\n")
        
        f.write(f"! Contaminants:\n")
        for cont in contaminants:
            Ap = cont['mass_number']  
            Zp = get_atomic_number(cont['species'])  
            q = cont['charge']  
            
            # Skip if this is the target beam (already written)
            if Ap == target_A and q == target_q:
                continue
            
            for E in energies:
                # Format: Ap Zp q E At Zt Thick d
                f.write(f"{Ap} {Zp} {q}  {E} {target_material['At']} "
                       f"{target_material['Zt']} {target_material['Thick']} "
                       f"{target_material['d']}\n")
    
    total_calcs = len(energies) + (len(contaminants) - 1) * len(energies)  # +1 for main beam, -1 if target in contaminants
    print(f"✓ Created ETACHA input file: {betacha_path}")
    print(f"  - Main beam + {len(contaminants)} contaminants × {len(energies)} energies = "
          f"{total_calcs} calculations")
    
    return betacha_path


def get_atomic_number(species_symbol):
    """Get atomic number (Z) from element symbol"""
    # Periodic table mapping
    Z_map = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,
        "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24,
        "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32,
        "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
        "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
        "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56,
        "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
        "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72,
        "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
        "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88,
        "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96,
        "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
        "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
        "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
    }
    return Z_map.get(species_symbol, 0)


def run_etacha_for_contaminants(betacha_file, contaminant_results, 
                                output_csv="contaminant_charge_distributions.csv",
                                min_ratio_percent=0.1):
    """
    Run ETACHA batch processing for contaminants and save results with abundance filtering
    
    Parameters:
    -----------
    betacha_file : str or Path
        Path to .betacha input file
    contaminant_results : dict
        Results from calculate_contaminants() containing abundance data
    output_csv : str
        Output CSV filename for results
    min_ratio_percent : float
        Minimum ratio (abundance × charge_fraction) in percent to include (default 0.1%)
    
    Returns:
    --------
    Path : Path to results CSV file with filtered charge states
    """
    if not Path(ETACHA_EXE).exists():
        print(f"ETACHA not found at: {ETACHA_EXE}")
        print("   Please update ETACHA_EXE path in etacha_batch_run.py")
        return None
    
    if not Path("template.etacha").exists():
        print("Missing template.etacha file")
        return None
    
    betacha_path = Path(betacha_file)
    
    # Create abundance lookup from contaminant results
    # Map by mass number only (isotope abundance doesn't depend on charge state)
    abundance_map = {}  # {A: abundance_%}
    for cont in contaminant_results['contaminants']:
        A = cont['mass_number']
        # Store the abundance for this mass number (same for all charge states)
        if A not in abundance_map:
            abundance_map[A] = cont.get('abundance_%', 0)
    
    # Read betacha file
    from etacha_batch_run import read_betacha
    rows_in = read_betacha(betacha_path)
    
    # Create case files
    print(f"\nCreating {len(rows_in)} ETACHA case files...")
    case_files = create_case_files(rows_in, betacha_path)
    
    # Run each case
    print(f"\nRunning ETACHA calculations...")
    results = []
    
    for i, (prefix, case_path, params) in enumerate(case_files, 1):
        print(f"  [{i}/{len(case_files)}] {prefix}")
        final_path = run_case(prefix, case_path)
        
        energy_res = None
        qmap = {}
        
        if final_path:
            if final_path.name.endswith("_forExcel.txt"):
                energy_res, qmap = parse_forExcel(final_path)
            else:
                energy_res, qmap = parse_eta(final_path)
        
        QM, QF, SIG, IQMAX = compute_moments(qmap)
        
        # Get natural abundance for this isotope (by mass number, not charge)
        A = int(params["Ap"])  # Convert from string to int for lookup
        q_initial = int(params["q"])
        abundance = abundance_map.get(A, 0)
        
        # Store the result with abundance and calculate combined ratios
        result_dict = {
            "Ab": A,
            "Zb": params["Zp"],
            "Qb": q_initial,
            "Eb": params["E"],
            "At": params["At"],
            "Zt": params["Zt"],
            "thick": params["Thick"],
            "density": params["d"],
            "abundance_%": abundance,  # Store abundance
            "QM": QM,
            "QF": QF,
            "sig(QF)": SIG,
            "iQmax": IQMAX,
            "EnergyResidual": energy_res,
            "EqThick-CS": "",
            "EqThick-Slope": "",
            "finalFile": str(final_path) if final_path else ""
        }
        
        # Store charge fractions and calculate combined ratios
        for q, charge_fraction in qmap.items():
            # Store charge fraction as percentage (qmap values are 0-1 decimals from ETACHA)
            result_dict[f"q{q}"] = charge_fraction * 100.0  # Convert to percentage for CSV
            
            # Calculate combined ratio = abundance% × charge_fraction(decimal)
            # Result is in percentage (e.g., 19.06 means 19.06%)
            combined_ratio = abundance * charge_fraction  # abundance is %, charge_fraction is 0-1
            result_dict[f"q{q}_ratio"] = combined_ratio
        
        results.append(result_dict)
    
    # Save results
    all_charge_states = set()
    for r in results:
        for key in r.keys():
            if key.startswith('q') and key[1:].lstrip('-').isdigit() and not key.endswith('_ratio'):
                q_num = key[1:].lstrip('-')
                if q_num.isdigit():
                    all_charge_states.add(int(q_num))
    
    charge_columns = [f"q{q}" for q in sorted(all_charge_states, reverse=True)]
    ratio_columns = [f"q{q}_ratio" for q in sorted(all_charge_states, reverse=True)]
    
    output_path = Path(output_csv)
    OUT_COLUMNS = [
        "Ab", "Zb", "Qb", "Eb",
        "At", "Zt", "thick", "density",
        "abundance_%",  # Add abundance column
        "QM", "QF", "sig(QF)", "iQmax",
        "EnergyResidual",
        "EqThick-CS", "EqThick-Slope",
        "finalFile",
    ] + charge_columns + ratio_columns  # Include both fractions and ratios
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        w.writeheader()
        for r in results:
            w.writerow(r)
    
    # Collect all charge state entries
    all_entries = []
    for r in results:
        A = r['Ab']
        Z = int(r['Zb'])  # Ensure Z is integer
        q_in = r['Qb']
        abundance = r.get('abundance_%', 0)
        species = get_element_symbol(Z)
        
        for q_state in sorted(all_charge_states, reverse=True):
            ratio_key = f"q{q_state}_ratio"
            charge_key = f"q{q_state}"
            # Filter for display: only show if ratio >= threshold
            if ratio_key in r:
                ratio = r.get(ratio_key, 0)
                if ratio >= min_ratio_percent:  # Apply threshold filter for display
                    all_entries.append({
                        'isotope': f"{species}-{A}",
                        'q_in': q_in,
                        'q_out': q_state,
                        'abundance': abundance,
                        'charge_fraction': r.get(charge_key, 0),
                        'ratio': ratio
                    })

    
    # Sort and display top 15
    all_entries.sort(key=lambda x: x['ratio'], reverse=True)
    for entry in all_entries[:15]:
        print(f"{entry['isotope']:<12} {entry['q_in']:<6} {entry['q_out']:<6} "
              f"{entry['abundance']:<10.3f} {entry['charge_fraction']:<10.3f} {entry['ratio']:<10.4f}")
    
    return output_path



def print_charge_distribution_summary(csv_file):
    """
    Read and print charge state distribution results from CSV
    
    Parameters:
    -----------
    csv_file : str or Path
        Path to results CSV file
    """
    csv_path = Path(csv_file)
    
    if not csv_path.exists():
        print(f"File not found: {csv_path}")
        return
    
    print("CHARGE STATE DISTRIBUTION RESULTS")
    
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        print("No results found in CSV file.")
        return
    
    print(f"\nTotal calculations: {len(rows)}")
    print("\n" + "-"*100)
    print(f"{'Isotope':<10} {'Zb':<4} {'Qb':<4} {'Energy':<8} {'Target':<8} "
          f"{'QM':<8} {'QF':<6} {'σ(QF)':<8} {'E_res':<10}")
    print("-"*100)
    
    for row in rows:
        Ab = row['Ab']
        Zb = row['Zb']
        Qb = row['Qb']
        Eb = row['Eb']
        At = row['At']
        Zt = row['Zt']
        
        # Get element symbol from Z
        elem = get_element_symbol(int(Zb)) if Zb else "?"
        isotope = f"{elem}-{Ab}"
        target = f"{get_element_symbol(int(Zt)) if Zt else '?'}-{At}"
        
        QM = float(row['QM']) if row['QM'] else None
        QF = int(float(row['QF'])) if row['QF'] else None
        sig = float(row['sig(QF)']) if row['sig(QF)'] else None
        E_res = float(row['EnergyResidual']) if row['EnergyResidual'] else None
        
        QM_str = f"{QM:.2f}" if QM is not None else "N/A"
        QF_str = f"{QF}" if QF is not None else "N/A"
        sig_str = f"{sig:.3f}" if sig is not None else "N/A"
        E_res_str = f"{E_res:.2f}" if E_res is not None else "N/A"
        
        print(f"{isotope:<10} {Zb:<4} {Qb:<4} {Eb:<8} {target:<8} "
              f"{QM_str:<8} {QF_str:<6} {sig_str:<8} {E_res_str:<10}")
    
    print("\nColumn definitions:")
    print("  QM = Mean charge state")
    print("  QF = Most frequent charge state")
    print("  sigma(QF) = Standard deviation of charge distribution")
    print("  E_res = Residual energy (MeV/u)")
    print()


def get_element_symbol(Z):
    """Get element symbol from atomic number"""
    symbols = [
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]
    return symbols[Z] if 0 <= Z < len(symbols) else "?"


def calculate_and_run_contaminants(target_A, target_q, beam_energy, target_material,
                                   tolerance_percent=1.0, species_list=None,
                                   output_prefix="contaminant_analysis",
                                   use_exact_mass=False):
    """
    Complete workflow: Calculate contaminants, create ETACHA input, run ETACHA simulations, 
    and display results
    
    Parameters:
    -----------
    target_A : int
        Target mass number
    target_q : int
        Target charge state
    beam_energy : float or list
        Beam energy in MeV/u (single value or list)
    target_material : dict
        {'At': int, 'Zt': int, 'Thick': float, 'd': float}
    tolerance_percent : float
        Tolerance for contaminant search (default 1.0%)
    species_list : list
        Species to check (default: Some pre-selected species)
    output_prefix : str
        Prefix for output files
    use_exact_mass : bool
        If True, use exact atomic mass. If False, use rounded mass number. (default: False)
    
    Returns:
    --------
    tuple : (contaminant_results, etacha_results_csv_path)
    """

    contaminant_results = calculate_contaminants(
        target_A=target_A,
        target_q=target_q,
        tolerance_percent=tolerance_percent,
        species_list=species_list,
        verbose=True,
        use_exact_mass=use_exact_mass
    )
    
    csv_file = f"{output_prefix}_list.csv"
    save_results_to_csv(contaminant_results, csv_file)
    
    if contaminant_results['num_contaminants'] == 0:
        print("\nNo contaminants found. Nothing to simulate.")
        return contaminant_results, None
    
    # Create ETACHA input file
    betacha_file = f"{output_prefix}.betacha"
    create_betacha_from_contaminants(
        contaminant_results,
        beam_energy,
        target_material,
        betacha_file
    )
    
    results_csv = f"{output_prefix}_charge_distributions.csv"
    etacha_results = run_etacha_for_contaminants(
        betacha_file, 
        contaminant_results,  # Pass contaminant results for abundance filtering
        results_csv,
        min_ratio_percent=0.1
    )
    
    if etacha_results:
        print_charge_distribution_summary(etacha_results)
    
    print(f"  • Contaminant list: {csv_file}")
    print(f"  • ETACHA input: {betacha_file}")
    if etacha_results:
        print(f"  • Charge distributions (with abundance filtering): {results_csv}")
    
    return contaminant_results, etacha_results


if __name__ == "__main__":
    # Example: U-238 at q=34 through Li target
    lithium_target = {
        'At': 6.94,      # Lithium
        'Zt': 3,         # Z = 3 (Lithium)
        'Thick': 0.905,  # mg/cm²
        'd': 0.51        # density g/cm³
    }
    
    contaminants, charge_results = calculate_and_run_contaminants(
        target_A=238,
        target_q=35,
        beam_energy=17.0,  # at stripper
        target_material=lithium_target,
        tolerance_percent=1.0,
        output_prefix="U238_analysis"
    )
