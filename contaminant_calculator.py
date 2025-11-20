"""
Contaminant Calculator for Heavy Ion Accelerators
Calculates possible contaminants based on A/q
Created from CSD_select.py by Jinyu Wan [wan@frib.msu.edu]
"""

import json
import argparse
import csv
from pathlib import Path

# Default species list to search for contaminants
DEFAULT_SPECIES = ["H", "He", "C", "N", "O", "Ne", "Si", "Al", "Ar", "Ca", "Ni", "Zn", "Se", "Kr", "Ag", "Xe", "Tm", "W", "Pt", "Hg", "Bi", "U"]

# Species database with isotopes: {mass_amu: abundance_%}
# Only naturally occurring isotopes with measurable abundances are included
species_database = {
    "H": {
        "isotopes": {1.007825: 99.9885, 2.014102: 0.0115},
        "possible_charges": range(1, 2)
    },
    "He": {
        "isotopes": {3.016029: 0.000137, 4.002603: 99.999863},
        "possible_charges": range(1, 3)
    },
    "Li": {
        "isotopes": {6.015122: 7.59, 7.016004: 92.41},
        "possible_charges": range(1, 4)
    },
    "Be": {
        "isotopes": {9.012182: 100},
        "possible_charges": range(1, 5)
    },
    "B": {
        "isotopes": {10.012937: 19.9, 11.009305: 80.1},
        "possible_charges": range(1, 6)
    },
    "C": {
        "isotopes": {12.000000: 98.93, 13.003355: 1.07},
        "possible_charges": range(1, 7)
    },
    "N": {
        "isotopes": {14.003074: 99.632, 15.000109: 0.368},
        "possible_charges": range(1, 8)
    },
    "O": {
        "isotopes": {15.994915: 99.757, 16.999132: 0.038, 17.999160: 0.205},
        "possible_charges": range(1, 9)
    },
    "F": {
        "isotopes": {18.998403: 100},
        "possible_charges": range(1, 10)
    },
    "Ne": {
        "isotopes": {19.992440: 90.48, 20.993847: 0.27, 21.991386: 9.25},
        "possible_charges": range(1, 11)
    },
    "Na": {
        "isotopes": {22.989770: 100},
        "possible_charges": range(1, 12)
    },
    "Mg": {
        "isotopes": {23.985042: 78.99, 24.985837: 10.00, 25.982593: 11.01},
        "possible_charges": range(1, 13)
    },
    "Al": {
        "isotopes": {26.981538: 100},
        "possible_charges": range(1, 14)
    },
    "Si": {
        "isotopes": {27.976927: 92.2297, 28.976495: 4.6832, 29.973770: 3.0872},
        "possible_charges": range(1, 15)
    },
    "P": {
        "isotopes": {30.973762: 100},
        "possible_charges": range(1, 16)
    },
    "S": {
        "isotopes": {31.972071: 94.93, 32.971458: 0.76, 33.967867: 4.29, 35.967081: 0.02},
        "possible_charges": range(1, 17)
    },
    "Cl": {
        "isotopes": {34.968853: 75.78, 36.965903: 24.22},
        "possible_charges": range(1, 18)
    },
    "Ar": {
        "isotopes": {35.967546: 0.3365, 37.962732: 0.0632, 39.962383: 99.6003},
        "possible_charges": range(1, 19)
    },
    "K": {
        "isotopes": {38.963707: 93.2581, 39.963999: 0.0117, 40.961826: 6.7302},
        "possible_charges": range(1, 20)
    },
    "Ca": {
        "isotopes": {39.962591: 96.941, 41.958618: 0.647, 42.958767: 0.135, 43.955481: 2.086, 45.953693: 0.004, 47.952534: 0.187},
        "possible_charges": range(1, 21)
    },
    "Sc": {
        "isotopes": {44.955910: 100},
        "possible_charges": range(1, 22)  # Z = 21
    },
    "Ti": {
        "isotopes": {
            45.952629: 8.25,
            46.951764: 7.44,
            47.947947: 73.72,
            48.947871: 5.41,
            49.944792: 5.18
        },
        "possible_charges": range(1, 23)  # Z = 22
    },
    "V": {
        "isotopes": {49.947163: 0.250, 50.943964: 99.750},
        "possible_charges": range(1, 24)  # Z = 23
    },
    "Cr": {
        "isotopes": {
            49.946050: 4.345,
            51.940512: 83.789,
            52.940654: 9.501,
            53.938885: 2.365
        },
        "possible_charges": range(1, 25)  # Z = 24
    },
    "Mn": {
        "isotopes": {54.938050: 100},
        "possible_charges": range(1, 26)  # Z = 25
    },
    "Fe": {
        "isotopes": {
            53.939615: 5.845,
            55.934942: 91.754,
            56.935399: 2.119,
            57.933280: 0.282
        },
        "possible_charges": range(1, 27)  # Z = 26
    },
    "Co": {
        "isotopes": {58.933200: 100},
        "possible_charges": range(1, 28)  # Z = 27
    },
    "Ni": {
        "isotopes": {57.935348: 68.0769, 59.930791: 26.2231, 60.931060: 1.1399, 61.928349: 3.6345, 63.927970: 0.9256},
        "possible_charges": range(1, 29)
    },
    "Zn": {
        "isotopes": {63.929147: 48.63, 65.926037: 27.90, 66.927131: 4.10, 67.924848: 18.75, 69.925325: 0.62},
        "possible_charges": range(1, 31)
    },
    "Se": {
        "isotopes": {73.922477: 0.89, 75.919214: 9.37, 76.919915: 7.63, 77.917310: 23.77, 79.916522: 49.61, 81.916700: 8.73},
        "possible_charges": range(1, 35)
    },
    "Kr": {
        "isotopes": {77.920386: 0.35, 79.916378: 2.28, 81.913485: 11.58, 82.914136: 11.49, 83.911507: 57.00, 85.910610: 17.30},
        "possible_charges": range(1, 37)
    },
    "Rb": {
        "isotopes": {84.911789: 72.17, 86.909183: 27.83},
        "possible_charges": range(1, 38)  # Z = 37
    },
    "Sr": {
        "isotopes": {
            83.913425: 0.56,
            85.909262: 9.86,
            86.908879: 7.00,
            87.905614: 82.58
        },
        "possible_charges": range(1, 39)  # Z = 38
    },
    "Y": {
        "isotopes": {88.905848: 100},
        "possible_charges": range(1, 40)  # Z = 39
    },
    "Zr": {
        "isotopes": {
            89.904704: 51.45,
            90.905645: 11.22,
            91.905040: 17.15,
            93.906316: 17.38,
            95.908276: 2.80
        },
        "possible_charges": range(1, 41)  # Z = 40
    },
    "Nb": {
        "isotopes": {92.906378: 100},
        "possible_charges": range(1, 42)  # Z = 41
    },
    "Mo": {
        "isotopes": {
            91.906810: 14.84,
            93.905088: 9.25,
            94.905841: 15.92,
            95.904679: 16.68,
            96.906021: 9.55,
            97.905408: 24.13,
            99.907477: 9.63
        },
        "possible_charges": range(1, 43)  # Z = 42
    },
    # Tc has only "*" abundance in the reference → skipped
    "Ru": {
        "isotopes": {
            95.907598: 5.54,
            97.905287: 1.87,
            98.905939: 12.76,
            99.904220: 12.60,
            100.905582: 17.06,
            101.904350: 31.55,
            103.905430: 18.62
        },
        "possible_charges": range(1, 45)  # Z = 44
    },
    "Rh": {
        "isotopes": {102.905504: 100},
        "possible_charges": range(1, 46)  # Z = 45
    },
    "Pd": {
        "isotopes": {
            101.905608: 1.02,
            103.904035: 11.14,
            104.905084: 22.33,
            105.903483: 27.33,
            107.903894: 26.46,
            109.905152: 11.72
        },
        "possible_charges": range(1, 47)  # Z = 46
    },
    "Ag": {
        "isotopes": {106.905093: 51.839, 108.904756: 48.161},
        "possible_charges": range(1, 48)
    },
    "Cd": {
        "isotopes": {
            105.906458: 1.25,
            107.904183: 0.89,
            109.903006: 12.49,
            110.904182: 12.80,
            111.902757: 24.13,
            112.904401: 12.22,
            113.903358: 28.73,
            115.904755: 7.49
        },
        "possible_charges": range(1, 49)  # Z = 48
    },
    "In": {
        "isotopes": {112.904061: 4.29, 114.903878: 95.71},
        "possible_charges": range(1, 50)  # Z = 49
    },
    "Sn": {
        "isotopes": {
            111.904821: 0.97,
            113.902782: 0.66,
            114.903346: 0.34,
            115.901744: 14.54,
            116.902954: 7.68,
            117.901606: 24.22,
            118.903309: 8.59,
            119.902197: 32.58,
            121.903440: 4.63,
            123.905275: 5.79
        },
        "possible_charges": range(1, 51)  # Z = 50
    },
    "Sb": {
        "isotopes": {120.903818: 57.21, 122.904216: 42.79},
        "possible_charges": range(1, 52)  # Z = 51
    },
    "Te": {
        "isotopes": {
            119.904020: 0.09,
            121.903047: 2.55,
            122.904273: 0.89,
            123.902819: 4.74,
            124.904425: 7.07,
            125.903306: 18.84,
            127.904461: 31.74,
            129.906223: 34.08
        },
        "possible_charges": range(1, 53)  # Z = 52
    },
    "I": {
        "isotopes": {126.904468: 100},
        "possible_charges": range(1, 54)  # Z = 53
    },
    "Xe": {
        "isotopes": {123.905896: 0.09, 125.904269: 0.09, 127.903530: 1.92, 128.904779: 26.44, 129.903508: 4.08, 130.905082: 21.18, 131.904154: 26.89, 133.905395: 10.44, 135.907220: 8.87},
        "possible_charges": range(1, 55)
    },
    "Cs": {
        "isotopes": {132.905447: 100},
        "possible_charges": range(1, 56)  # Z = 55
    },
    "Ba": {
        "isotopes": {
            129.906310: 0.106,
            131.905056: 0.101,
            133.904503: 2.417,
            134.905683: 6.592,
            135.904570: 7.854,
            136.905821: 11.232,
            137.905241: 71.698
        },
        "possible_charges": range(1, 57)  # Z = 56
    },
    "La": {
        "isotopes": {137.907107: 0.090, 138.906348: 99.910},
        "possible_charges": range(1, 58)  # Z = 57
    },
    "Ce": {
        "isotopes": {
            135.907144: 0.185,
            137.905986: 0.251,
            139.905434: 88.450,
            141.909240: 11.114
        },
        "possible_charges": range(1, 59)  # Z = 58
    },
    "Pr": {
        "isotopes": {140.907648: 100},
        "possible_charges": range(1, 60)  # Z = 59
    },
    "Nd": {
        "isotopes": {
            141.907719: 27.2,
            142.909810: 12.2,
            143.910083: 23.8,
            144.912569: 8.3,
            145.913112: 17.2,
            147.916889: 5.7,
            149.920887: 5.6
        },
        "possible_charges": range(1, 61)  # Z = 60
    },
    # Pm only has "*" abundance → skipped
    "Sm": {
        "isotopes": {
            143.911995: 3.07,
            146.914893: 14.99,
            147.914818: 11.24,
            148.917180: 13.82,
            149.917271: 7.38,
            151.919728: 26.75,
            153.922205: 22.75
        },
        "possible_charges": range(1, 63)  # Z = 62
    },
    "Eu": {
        "isotopes": {150.919846: 47.81, 152.921226: 52.19},
        "possible_charges": range(1, 64)  # Z = 63
    },
    "Gd": {
        "isotopes": {
            151.919788: 0.20,
            153.920862: 2.18,
            154.922619: 14.80,
            155.922120: 20.47,
            156.923957: 15.65,
            157.924101: 24.84,
            159.927051: 21.86
        },
        "possible_charges": range(1, 65)  # Z = 64
    },
    "Tb": {
        "isotopes": {158.925343: 100},
        "possible_charges": range(1, 66)  # Z = 65
    },
    "Dy": {
        "isotopes": {
            155.924278: 0.06,
            157.924405: 0.10,
            159.925194: 2.34,
            160.926930: 18.91,
            161.926795: 25.51,
            162.928728: 24.90,
            163.929171: 28.18
        },
        "possible_charges": range(1, 67)  # Z = 66
    },
    "Ho": {
        "isotopes": {164.930319: 100},
        "possible_charges": range(1, 68)  # Z = 67
    },
    "Er": {
        "isotopes": {
            161.928775: 0.14,
            163.929197: 1.61,
            165.930290: 33.61,
            166.932045: 22.93,
            167.932368: 26.78,
            169.935460: 14.93
        },
        "possible_charges": range(1, 69)  # Z = 68
    },
    "Tm": {
        "isotopes": {168.934211: 100},
        "possible_charges": range(1, 70)
    },
    "Yb": {
        "isotopes": {
            167.933894: 0.13,
            169.934759: 3.04,
            170.936322: 14.28,
            171.936378: 21.83,
            172.938207: 16.13,
            173.938858: 31.83,
            175.942568: 12.76
        },
        "possible_charges": range(1, 71)  # Z = 70
    },
    "Lu": {
        "isotopes": {174.940768: 97.41, 175.942682: 2.59},
        "possible_charges": range(1, 72)  # Z = 71
    },
    "Hf": {
        "isotopes": {
            173.940040: 0.16,
            175.941402: 5.26,
            176.943220: 18.60,
            177.943698: 27.28,
            178.945815: 13.62,
            179.946549: 35.08
        },
        "possible_charges": range(1, 73)  # Z = 72
    },
    "Ta": {
        "isotopes": {179.947466: 0.012, 180.947996: 99.988},
        "possible_charges": range(1, 74)  # Z = 73
    },
    "W": {
        "isotopes": {179.946706: 0.12, 181.948206: 26.50, 182.950224: 14.31, 183.950933: 30.64, 185.954362: 28.43},
        "possible_charges": range(1, 75)
    },
    "Re": {
        "isotopes": {184.952956: 37.40, 186.955751: 62.60},
        "possible_charges": range(1, 76)  # Z = 75
    },
    "Os": {
        "isotopes": {
            183.952491: 0.02,
            185.953838: 1.59,
            186.955748: 1.96,
            187.955836: 13.24,
            188.958145: 16.15,
            189.958445: 26.26,
            191.961479: 40.78
        },
        "possible_charges": range(1, 77)  # Z = 76
    },
    "Ir": {
        "isotopes": {190.960591: 37.3, 192.962924: 62.7},
        "possible_charges": range(1, 78)  # Z = 77
    },
    "Pt": {
        "isotopes": {189.959932: 0.014, 191.961035: 0.782, 193.962664: 32.967, 194.964774: 33.832, 195.964935: 25.242, 197.967876: 7.163},
        "possible_charges": range(1, 79)
    },
    "Au": {
        "isotopes": {196.966552: 100},
        "possible_charges": range(1, 80)  # Z = 79
    },
    "Hg": {
        "isotopes": {195.965815: 0.15, 197.966752: 9.97, 198.968262: 16.87, 199.968309: 23.10, 200.970285: 13.18, 201.970626: 29.86, 203.973476: 6.87},
        "possible_charges": range(1, 81)
    },
    "Tl": {
        "isotopes": {202.972329: 29.524, 204.974412: 70.476},
        "possible_charges": range(1, 82)  # Z = 81
    },
    "Pb": {
        "isotopes": {
            203.973029: 1.4,
            205.974449: 24.1,
            206.975881: 22.1,
            207.976636: 52.4
        },
        "possible_charges": range(1, 83)  # Z = 82
    },
    "Bi": {
        "isotopes": {208.980383: 100},
        "possible_charges": range(1, 84)
    },
    # Po–Ac in the table are purely radiogenic (*), so left out
    "Th": {
        "isotopes": {232.038050: 100},
        "possible_charges": range(1, 91)  # Z = 90
    },
    "Pa": {
        "isotopes": {231.035879: 100},
        "possible_charges": range(1, 92)  # Z = 91
    },
    "U": {
        "isotopes": {234.040946: 0.0055, 235.043923: 0.7200, 238.050783: 99.2745},
        "possible_charges": range(1, 93)
    },
}


def calculate_A_over_q(mass_amu, charge):
    """Calculate A/q ratio (beam rigidity)"""
    return mass_amu / charge



def find_contaminants(target_mass_amu, target_q, tolerance_percent=1.0, species_list=None, 
                      custom_charge_ranges=None, use_exact_mass=False):
    """
    Find potential contaminants within tolerance of target beam rigidity
    
    Parameters:
    -----------
    target_mass_amu : float
        Target mass in atomic mass units (can use accurate mass or mass number)
    target_q : int
        Target charge state
    tolerance_percent : float
        Tolerance percentage (default 1.0%)
    species_list : list
        List of species symbols to check (default: DEFAULT_SPECIES)
    custom_charge_ranges : dict
        Custom charge ranges per species (optional)
    use_exact_mass : bool
        If True, use exact atomic mass. If False, use rounded mass number. (default: False)
    
    Returns:
    --------
    list of dict : Contaminants within tolerance
    """
    if species_list is None:
        species_list = DEFAULT_SPECIES
    
    target_Aq = calculate_A_over_q(target_mass_amu, target_q)
    contaminants = []
    
    for sp in species_list:
        if sp not in species_database:
            print(f"Warning: Species '{sp}' not in database, skipping...")
            continue
        
        # Get isotopes dict with {mass_amu: abundance_%}
        isotopes = species_database[sp]["isotopes"]
        
        # Use custom charge range if provided, otherwise use database default
        if custom_charge_ranges and sp in custom_charge_ranges:
            charge_range = custom_charge_ranges[sp]
        else:
            charge_range = species_database[sp]["possible_charges"]
        
        for mass_amu, abundance in isotopes.items():
            mass_number = int(round(mass_amu))  # For display purposes
            
            # Determine mass to use for calculation
            calc_mass = mass_amu if use_exact_mass else mass_number
            
            for q_sp in charge_range:
                aq_sp = calculate_A_over_q(calc_mass, q_sp)
                rel_diff = abs(aq_sp - target_Aq) / target_Aq
                
                if rel_diff <= tolerance_percent / 100.0:
                    contaminants.append({
                        "species": sp,
                        "isotope": f"{sp}-{mass_number}",
                        "mass_number": mass_number,
                        "mass_amu": mass_amu,  # Store accurate mass
                        "charge": q_sp,
                        "A/q": aq_sp,
                        "rel_diff_%": rel_diff * 100,
                        "abs_diff": abs(aq_sp - target_Aq),
                        "abundance_%": abundance  # Natural abundance percentage
                    })
    
    # Sort by absolute difference from target
    contaminants.sort(key=lambda x: x["abs_diff"])
    
    return contaminants, target_Aq


def calculate_contaminants(target_A, target_q, tolerance_percent=1.0, species_list=None, 
                          custom_charge_ranges=None, verbose=True, use_exact_mass=False):
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
    use_exact_mass : bool
        If True, use exact atomic mass. If False, use rounded mass number. (default: False)
    
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
        target_A, target_q, tolerance_percent, species_list, custom_charge_ranges, use_exact_mass
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
        writer.writerow(['Species', 'Isotope', 'Mass Number (A)', 'Mass (amu)', 'Charge State (q)', 
                        'A/q Ratio', 'Relative Difference (%)', 'Absolute Difference', 'Natural Abundance (%)'])
        
        # Write contaminant data
        for cont in contaminants:
            writer.writerow([
                cont['species'],
                cont['isotope'],
                cont['mass_number'],
                f"{cont['mass_amu']:.6f}",
                cont['charge'],
                f"{cont['A/q']:.6f}",
                f"{cont['rel_diff_%']:.6f}",
                f"{cont['abs_diff']:.6f}",
                f"{cont.get('abundance_%', 0):.4f}"
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

  # Use exact atomic mass instead of rounded mass number
  python contaminant_calculator.py --target-A 238 --target-q 34 --use-exact-mass
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
    parser.add_argument('--use-exact-mass', action='store_true',
                       help='Use exact atomic mass for calculation instead of rounded mass number')
    
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
        verbose=True,
        use_exact_mass=args.use_exact_mass
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
