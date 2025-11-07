"""
Contaminant Calculator for Heavy Ion Accelerators
Calculates possible contaminants based on A/q
Created from CSD_select.py GUI - command-line version by Jinyu Wan [wan@frib.msu.edu]
"""

import json
import argparse
import csv
from pathlib import Path

# Default species list to search for contaminants
# This is a simplified list for the convenience of viewing results.
# Users can modify this list as needed.
DEFAULT_SPECIES = ["H", "C", "N", "O", "Ne", "Ca", "Ag", "Xe", "Ar", "Si", "W", "U"]

# Species database with mass numbers and possible charge states
species_database = {
    "H": {"mass_numbers": [1, 2, 3], "possible_charges": range(1, 2)},
    "He": {"mass_numbers": [3, 4], "possible_charges": range(1, 3)},
    "Li": {"mass_numbers": [6, 7], "possible_charges": range(1, 4)},
    "Be": {"mass_numbers": [9], "possible_charges": range(1, 5)},
    "B": {"mass_numbers": [10, 11], "possible_charges": range(1, 6)},
    "C": {"mass_numbers": [12, 13, 14], "possible_charges": range(1, 7)},
    "N": {"mass_numbers": [14, 15], "possible_charges": range(1, 8)},
    "O": {"mass_numbers": [16, 17, 18], "possible_charges": range(1, 9)},
    "F": {"mass_numbers": [19], "possible_charges": range(1, 10)},
    "Ne": {"mass_numbers": [20, 21, 22], "possible_charges": range(1, 11)},
    "Na": {"mass_numbers": [23], "possible_charges": range(1, 12)},
    "Mg": {"mass_numbers": [24, 25, 26], "possible_charges": range(1, 13)},
    "Al": {"mass_numbers": [27], "possible_charges": range(1, 14)},
    "Si": {"mass_numbers": [28, 29, 30], "possible_charges": range(1, 15)},
    "P": {"mass_numbers": [31], "possible_charges": range(1, 16)},
    "S": {"mass_numbers": [32, 33, 34, 36], "possible_charges": range(1, 17)},
    "Cl": {"mass_numbers": [35, 37], "possible_charges": range(1, 18)},
    "Ar": {"mass_numbers": [36, 38, 40], "possible_charges": range(1, 19)},
    "K": {"mass_numbers": [39, 40, 41], "possible_charges": range(1, 20)},
    "Ca": {"mass_numbers": [40, 42, 43, 44, 46, 48], "possible_charges": range(1, 21)},
    "Sc": {"mass_numbers": [45], "possible_charges": range(1, 22)},
    "Ti": {"mass_numbers": [46, 47, 48, 49, 50], "possible_charges": range(1, 23)},
    "V": {"mass_numbers": [50, 51], "possible_charges": range(1, 24)},
    "Cr": {"mass_numbers": [50, 52, 53, 54], "possible_charges": range(1, 25)},
    "Mn": {"mass_numbers": [55], "possible_charges": range(1, 26)},
    "Fe": {"mass_numbers": [54, 56, 57, 58], "possible_charges": range(1, 27)},
    "Co": {"mass_numbers": [59], "possible_charges": range(1, 28)},
    "Ni": {"mass_numbers": [58, 60, 61, 62, 64], "possible_charges": range(1, 29)},
    "Cu": {"mass_numbers": [63, 65], "possible_charges": range(1, 30)},
    "Zn": {"mass_numbers": [64, 66, 67, 68, 70], "possible_charges": range(1, 31)},
    "Ga": {"mass_numbers": [69, 71], "possible_charges": range(1, 32)},
    "Ge": {"mass_numbers": [70, 72, 73, 74, 76], "possible_charges": range(1, 33)},
    "As": {"mass_numbers": [75], "possible_charges": range(1, 34)},
    "Se": {"mass_numbers": [74, 76, 77, 78, 80, 82], "possible_charges": range(1, 35)},
    "Br": {"mass_numbers": [79, 81], "possible_charges": range(1, 36)},
    "Kr": {"mass_numbers": [78, 80, 82, 83, 84, 86], "possible_charges": range(1, 37)},
    "Rb": {"mass_numbers": [85, 87], "possible_charges": range(1, 38)},
    "Sr": {"mass_numbers": [84, 86, 87, 88], "possible_charges": range(1, 39)},
    "Y": {"mass_numbers": [89], "possible_charges": range(1, 40)},
    "Zr": {"mass_numbers": [90, 91, 92, 94, 96], "possible_charges": range(1, 41)},
    "Nb": {"mass_numbers": [93], "possible_charges": range(1, 42)},
    "Mo": {"mass_numbers": [92, 94, 95, 96, 97, 98, 100], "possible_charges": range(1, 43)},
    "Tc": {"mass_numbers": [97, 98, 99], "possible_charges": range(1, 44)},
    "Ru": {"mass_numbers": [96, 98, 99, 100, 101, 102, 104], "possible_charges": range(1, 45)},
    "Rh": {"mass_numbers": [103], "possible_charges": range(1, 46)},
    "Pd": {"mass_numbers": [102, 104, 105, 106, 108, 110], "possible_charges": range(1, 47)},
    "Ag": {"mass_numbers": [107, 109], "possible_charges": range(1, 48)},
    "Cd": {"mass_numbers": [106, 108, 110, 111, 112, 113, 114, 116], "possible_charges": range(1, 49)},
    "In": {"mass_numbers": [113, 115], "possible_charges": range(1, 50)},
    "Sn": {"mass_numbers": [112, 114, 115, 116, 117, 118, 119, 120, 122, 124], "possible_charges": range(1, 51)},
    "Sb": {"mass_numbers": [121, 123], "possible_charges": range(1, 52)},
    "Te": {"mass_numbers": [120, 122, 123, 124, 125, 126, 128, 130], "possible_charges": range(1, 53)},
    "I": {"mass_numbers": [127], "possible_charges": range(1, 54)},
    "Xe": {"mass_numbers": [124, 126, 128, 129, 130, 131, 132, 134, 136], "possible_charges": range(1, 55)},
    "Cs": {"mass_numbers": [133], "possible_charges": range(1, 56)},
    "Ba": {"mass_numbers": [130, 132, 134, 135, 136, 137, 138], "possible_charges": range(1, 57)},
    "La": {"mass_numbers": [138, 139], "possible_charges": range(1, 58)},
    "Ce": {"mass_numbers": [136, 138, 140, 142], "possible_charges": range(1, 59)},
    "Pr": {"mass_numbers": [141], "possible_charges": range(1, 60)},
    "Nd": {"mass_numbers": [142, 143, 144, 145, 146, 148, 150], "possible_charges": range(1, 61)},
    "Pm": {"mass_numbers": [145, 147], "possible_charges": range(1, 62)},
    "Sm": {"mass_numbers": [144, 147, 148, 149, 150, 152, 154], "possible_charges": range(1, 63)},
    "Eu": {"mass_numbers": [151, 153], "possible_charges": range(1, 64)},
    "Gd": {"mass_numbers": [152, 154, 155, 156, 157, 158, 160], "possible_charges": range(1, 65)},
    "Tb": {"mass_numbers": [159], "possible_charges": range(1, 66)},
    "Dy": {"mass_numbers": [156, 158, 160, 161, 162, 163, 164], "possible_charges": range(1, 67)},
    "Ho": {"mass_numbers": [165], "possible_charges": range(1, 68)},
    "Er": {"mass_numbers": [162, 164, 166, 167, 168, 170], "possible_charges": range(1, 69)},
    "Tm": {"mass_numbers": [169], "possible_charges": range(1, 70)},
    "Yb": {"mass_numbers": [168, 170, 171, 172, 173, 174, 176], "possible_charges": range(1, 71)},
    "Lu": {"mass_numbers": [175, 176], "possible_charges": range(1, 72)},
    "Hf": {"mass_numbers": [174, 176, 177, 178, 179, 180], "possible_charges": range(1, 73)},
    "Ta": {"mass_numbers": [180, 181], "possible_charges": range(1, 74)},
    "W": {"mass_numbers": [180, 182, 183, 184, 186], "possible_charges": range(1, 75)},
    "Re": {"mass_numbers": [185, 187], "possible_charges": range(1, 76)},
    "Os": {"mass_numbers": [184, 186, 187, 188, 189, 190, 192], "possible_charges": range(1, 77)},
    "Ir": {"mass_numbers": [191, 193], "possible_charges": range(1, 78)},
    "Pt": {"mass_numbers": [190, 192, 194, 195, 196, 198], "possible_charges": range(1, 79)},
    "Au": {"mass_numbers": [197], "possible_charges": range(1, 80)},
    "Hg": {"mass_numbers": [196, 198, 199, 200, 201, 202, 204], "possible_charges": range(1, 81)},
    "Tl": {"mass_numbers": [203, 205], "possible_charges": range(1, 82)},
    "Pb": {"mass_numbers": [204, 206, 207, 208], "possible_charges": range(1, 83)},
    "Bi": {"mass_numbers": [209], "possible_charges": range(1, 84)},
    "Po": {"mass_numbers": [209, 210], "possible_charges": range(1, 85)},
    "At": {"mass_numbers": [210, 211], "possible_charges": range(1, 86)},
    "Rn": {"mass_numbers": [211, 220, 222], "possible_charges": range(1, 87)},
    "Fr": {"mass_numbers": [223], "possible_charges": range(1, 88)},
    "Ra": {"mass_numbers": [223, 224, 226, 228], "possible_charges": range(1, 89)},
    "Ac": {"mass_numbers": [227], "possible_charges": range(1, 90)},
    "Th": {"mass_numbers": [230, 232], "possible_charges": range(1, 91)},
    "Pa": {"mass_numbers": [231], "possible_charges": range(1, 92)},
    "U": {"mass_numbers": [233, 234, 235, 236, 238], "possible_charges": range(1, 93)},
    "Np": {"mass_numbers": [236, 237], "possible_charges": range(1, 94)},
    "Pu": {"mass_numbers": [238, 239, 240, 241, 242, 244], "possible_charges": range(1, 95)},
    "Am": {"mass_numbers": [241, 243], "possible_charges": range(1, 96)},
    "Cm": {"mass_numbers": [243, 244, 245, 246, 247, 248], "possible_charges": range(1, 97)},
    "Bk": {"mass_numbers": [247, 249], "possible_charges": range(1, 98)},
    "Cf": {"mass_numbers": [249, 250, 251, 252], "possible_charges": range(1, 99)},
    "Es": {"mass_numbers": [252], "possible_charges": range(1, 100)},
    "Fm": {"mass_numbers": [257], "possible_charges": range(1, 101)},
    "Md": {"mass_numbers": [258, 260], "possible_charges": range(1, 102)},
    "No": {"mass_numbers": [259], "possible_charges": range(1, 103)},
    "Lr": {"mass_numbers": [262], "possible_charges": range(1, 104)},
    "Rf": {"mass_numbers": [267], "possible_charges": range(1, 105)},
    "Db": {"mass_numbers": [268], "possible_charges": range(1, 106)},
    "Sg": {"mass_numbers": [271], "possible_charges": range(1, 107)},
    "Bh": {"mass_numbers": [272], "possible_charges": range(1, 108)},
    "Hs": {"mass_numbers": [270], "possible_charges": range(1, 109)},
    "Mt": {"mass_numbers": [276], "possible_charges": range(1, 110)},
    "Ds": {"mass_numbers": [281], "possible_charges": range(1, 111)},
    "Rg": {"mass_numbers": [280], "possible_charges": range(1, 112)},
    "Cn": {"mass_numbers": [285], "possible_charges": range(1, 113)},
    "Nh": {"mass_numbers": [284], "possible_charges": range(1, 114)},
    "Fl": {"mass_numbers": [289], "possible_charges": range(1, 115)},
    "Mc": {"mass_numbers": [288], "possible_charges": range(1, 116)},
    "Lv": {"mass_numbers": [293], "possible_charges": range(1, 117)},
    "Ts": {"mass_numbers": [292], "possible_charges": range(1, 118)},
    "Og": {"mass_numbers": [294], "possible_charges": range(1, 119)},
}


def calculate_A_over_q(mass_number, charge):
    """Calculate A/q ratio (beam rigidity)"""
    return mass_number / charge


def find_contaminants(target_A, target_q, tolerance_percent=1.0, species_list=None, 
                      custom_charge_ranges=None):
    """
    Find potential contaminants within tolerance of target beam rigidity
    
    Parameters:
    -----------
    target_A : int
        Target mass number
    target_q : int
        Target charge state
    tolerance_percent : float
        Tolerance percentage (default 1.0%)
    species_list : list
        List of species symbols to check (default: DEFAULT_SPECIES)
    custom_charge_ranges : dict
        Custom charge ranges per species (optional)
    
    Returns:
    --------
    list of dict : Contaminants within tolerance
    """
    if species_list is None:
        species_list = DEFAULT_SPECIES
    
    target_Aq = calculate_A_over_q(target_A, target_q)
    contaminants = []
    
    for sp in species_list:
        if sp not in species_database:
            print(f"Warning: Species '{sp}' not in database, skipping...")
            continue
            
        mass_numbers = species_database[sp]["mass_numbers"]
        
        # Use custom charge range if provided, otherwise use database default
        if custom_charge_ranges and sp in custom_charge_ranges:
            charge_range = custom_charge_ranges[sp]
        else:
            charge_range = species_database[sp]["possible_charges"]
        
        for A_sp in mass_numbers:
            for q_sp in charge_range:
                aq_sp = calculate_A_over_q(A_sp, q_sp)
                rel_diff = abs(aq_sp - target_Aq) / target_Aq
                
                if rel_diff <= tolerance_percent / 100.0:
                    contaminants.append({
                        "species": sp,
                        "isotope": f"{sp}-{A_sp}",
                        "mass_number": A_sp,
                        "charge": q_sp,
                        "A/q": aq_sp,
                        "rel_diff_%": rel_diff * 100,
                        "abs_diff": abs(aq_sp - target_Aq)
                    })
    
    # Sort by absolute difference from target
    contaminants.sort(key=lambda x: x["abs_diff"])
    
    return contaminants, target_Aq


def calculate_contaminants(target_A, target_q, tolerance_percent=1.0, species_list=None, 
                          custom_charge_ranges=None, verbose=True):
    """
    Main callable function to calculate potential contaminants
    
    Parameters:
    -----------
    target_A : int
        Target mass number
    target_q : int
        Target charge state
    tolerance_percent : float
        Tolerance percentage (default 1.0%)
    species_list : list
        List of species symbols to check (default: DEFAULT_SPECIES)
    custom_charge_ranges : dict
        Custom charge ranges per species (optional)
    verbose : bool
        If True, print results to console (default: True)
    
    Returns:
    --------
    dict : Results dictionary containing:
        - 'target_A': target mass number
        - 'target_q': target charge state
        - 'target_Aq': target A/q ratio
        - 'tolerance_percent': tolerance used
        - 'contaminants': list of contaminant dictionaries
        - 'num_contaminants': number of contaminants found
    """
    contaminants, target_Aq = find_contaminants(
        target_A, target_q, tolerance_percent, species_list, custom_charge_ranges
    )
    
    results = {
        'target_A': target_A,
        'target_q': target_q,
        'target_Aq': target_Aq,
        'tolerance_percent': tolerance_percent,
        'contaminants': contaminants,
        'num_contaminants': len(contaminants)
    }
    
    if verbose:
        print_results(contaminants, target_Aq, target_A, target_q, tolerance_percent)
    
    return results


def print_results(contaminants, target_Aq, target_A, target_q, tolerance_percent):
    """Print contaminant results in a formatted table"""
    print(f"\nTarget Beam: A={target_A}, q={target_q}, A/q={target_Aq:.4f}")
    print(f"Tolerance: ±{tolerance_percent:.2f}%")
    print(f"A/q Range: {target_Aq*(1-tolerance_percent/100):.4f} - {target_Aq*(1+tolerance_percent/100):.4f}")
    print("\n" + "-"*80)
    
    if not contaminants:
        print("\nNo contaminants found within specified tolerance.")
        return
    
    print(f"\nFound {len(contaminants)} potential contaminant(s):\n")
    print(f"{'Isotope':<12} {'Charge':<8} {'A/q':<10} {'Diff (%)':<12} {'Abs Diff':<10}")
    print("-"*80)
    
    for cont in contaminants:
        print(f"{cont['isotope']:<12} {cont['charge']:<8} {cont['A/q']:<10.4f} "
              f"{cont['rel_diff_%']:<12.4f} {cont['abs_diff']:<10.6f}")

def save_results_to_csv(results, filename="contaminants_output.csv"):
    """
    Save results to a CSV file with clear column headers
    
    Parameters:
    -----------
    results : dict
        Results dictionary from calculate_contaminants()
    filename : str
        Output filename (default: "contaminants_output.csv")
    """
    contaminants = results['contaminants']
    target_Aq = results['target_Aq']
    target_A = results['target_A']
    target_q = results['target_q']
    tolerance_percent = results['tolerance_percent']
    
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header information
        writer.writerow(['Target Beam Information'])
        writer.writerow(['Target Mass Number (A)', target_A])
        writer.writerow(['Target Charge State (q)', target_q])
        writer.writerow(['Target A/q Ratio', f'{target_Aq:.6f}'])
        writer.writerow(['Tolerance (%)', f'{tolerance_percent:.2f}'])
        writer.writerow(['A/q Range Min', f'{target_Aq*(1-tolerance_percent/100):.6f}'])
        writer.writerow(['A/q Range Max', f'{target_Aq*(1+tolerance_percent/100):.6f}'])
        writer.writerow(['Number of Contaminants Found', len(contaminants)])
        writer.writerow([])  # Empty row for separation
        
        if not contaminants:
            writer.writerow(['No contaminants found within specified tolerance'])
            print(f"Results saved to: {filename}")
            return
        
        # Write contaminants table header
        writer.writerow(['Contaminant Data'])
        writer.writerow(['Species', 'Isotope', 'Mass Number (A)', 'Charge State (q)', 
                        'A/q Ratio', 'Relative Difference (%)', 'Absolute Difference'])
        
        # Write contaminant data
        for cont in contaminants:
            writer.writerow([
                cont['species'],
                cont['isotope'],
                cont['mass_number'],
                cont['charge'],
                f"{cont['A/q']:.6f}",
                f"{cont['rel_diff_%']:.6f}",
                f"{cont['abs_diff']:.6f}"
            ])
    
    print(f"Results saved to: {filename}")


def save_results_to_file(results, filename="contaminants_output.txt"):
    """
    Save results to a text file
    
    Parameters:
    -----------
    results : dict
        Results dictionary from calculate_contaminants()
    filename : str
        Output filename (default: "contaminants_output.txt")
    """
    contaminants = results['contaminants']
    target_Aq = results['target_Aq']
    target_A = results['target_A']
    target_q = results['target_q']
    tolerance_percent = results['tolerance_percent']
    
    with open(filename, 'w') as f:
        f.write(f"\nTarget Beam: A={target_A}, q={target_q}, A/q={target_Aq:.4f}\n")
        f.write(f"Tolerance: ±{tolerance_percent:.2f}%\n")
        f.write(f"A/q Range: {target_Aq*(1-tolerance_percent/100):.4f} - {target_Aq*(1+tolerance_percent/100):.4f}\n")
        f.write("\n" + "-"*80 + "\n")
        
        if not contaminants:
            f.write("\nNo contaminants found within specified tolerance.\n")
            return
        
        f.write(f"\nFound {len(contaminants)} potential contaminant(s):\n\n")
        f.write(f"{'Isotope':<12} {'Charge':<8} {'A/q':<10} {'Diff (%)':<12} {'Abs Diff':<10}\n")
        f.write("-"*80 + "\n")
        
        for cont in contaminants:
            f.write(f"{cont['isotope']:<12} {cont['charge']:<8} {cont['A/q']:<10.4f} "
                   f"{cont['rel_diff_%']:<12.4f} {cont['abs_diff']:<10.6f}\n")
        
        f.write("="*80 + "\n")
    
    print(f"Results saved to: {filename}")


# Command-line interface (optional)
def main():
    parser = argparse.ArgumentParser(
        description="Calculate potential contaminants in heavy ion accelerator based on beam rigidity (A/q)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with U-238 at charge state 34
  python contaminant_calculator.py --target-A 238 --target-q 34

  # With custom tolerance
  python contaminant_calculator.py --target-A 238 --target-q 34 --tolerance 0.5

  # Specify custom species to check
  python contaminant_calculator.py --target-A 238 --target-q 34 --species U Xe W

  # Check all available species
  python contaminant_calculator.py --target-A 238 --target-q 34 --all-species

  # Save results to CSV file (recommended)
  python contaminant_calculator.py --target-A 238 --target-q 34 --output results.csv

  # Save results to text file
  python contaminant_calculator.py --target-A 238 --target-q 34 --output results.txt
        """
    )
    
    parser.add_argument('--target-A', type=int, required=True,
                       help='Target mass number (A)')
    parser.add_argument('--target-q', type=int, required=True,
                       help='Target charge state (q)')
    parser.add_argument('--tolerance', type=float, default=1.0,
                       help='Tolerance in percent (default: 1.0%%)')
    parser.add_argument('--species', nargs='+', default=None,
                       help='List of species to check (e.g., U Xe W). Default: FRIB common species')
    parser.add_argument('--all-species', action='store_true',
                       help='Check all available species in database')
    parser.add_argument('--output', '-o', type=str, default=None,
                       help='Save results to file. Use .csv extension for CSV format (recommended) or .txt for text format')
    parser.add_argument('--list-species', action='store_true',
                       help='List all available species and exit')
    
    args = parser.parse_args()
    
    # Handle list-species flag
    if args.list_species:
        print("\nAvailable species in database:")
        print("-" * 40)
        all_species = sorted(species_database.keys())
        for i, sp in enumerate(all_species, 1):
            print(f"{sp:>3}", end="  ")
            if i % 10 == 0:
                print()
        print("\n")
        print(f"\nDefault FRIB species: {', '.join(DEFAULT_SPECIES)}")
        return
    
    # Determine species list
    if args.all_species:
        species_list = list(species_database.keys())
        print(f"Checking all {len(species_list)} species in database...")
    elif args.species:
        species_list = args.species
    else:
        species_list = DEFAULT_SPECIES
    
    # Calculate contaminants using the main function
    results = calculate_contaminants(
        args.target_A, 
        args.target_q, 
        args.tolerance,
        species_list,
        verbose=True
    )
    
    # Save to file if requested
    if args.output:
        # Determine file format based on extension
        if args.output.lower().endswith('.csv'):
            save_results_to_csv(results, args.output)
        else:
            save_results_to_file(results, args.output)


if __name__ == "__main__":
    main()
