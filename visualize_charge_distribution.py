"""
Visualize charge state distributions before and after stripping.
Shows contaminants and main beam species.
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def get_element_symbol(Z):
    """Convert atomic number to element symbol"""
    symbols = [
        "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
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
    if 0 <= Z < len(symbols):
        return symbols[Z]
    return f"Z{Z}"


def get_Z_from_symbol(symbol):
    """Convert element symbol to atomic number"""
    symbols = [
        "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
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
    try:
        return symbols.index(symbol)
    except ValueError:
        return 0


def plot_charge_distribution(contaminants_csv, charge_dist_csv, output_image="charge_states.png"):
    """
    Plot charge state distribution before and after stripping.
    Uses bar charts with labels at x-ticks.
    
    Args:
        contaminants_csv: CSV file with contaminant list
        charge_dist_csv: CSV file with ETACHA charge distributions
        output_image: Output PNG file name
    """
    before_data = []
    target_A = None
    target_q = None
    target_Aq = None
    
    with open(contaminants_csv, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'Target Mass Number (A)' in line:
                target_A = int(line.split(',')[1].strip())
            elif 'Target Charge State (q)' in line:
                target_q = int(line.split(',')[1].strip())
            elif 'Target A/q Ratio' in line:
                target_Aq = float(line.split(',')[1].strip())
        
        # Find where data starts
        data_start_idx = 0
        for i, line in enumerate(lines):
            if 'Contaminant Data' in line:
                data_start_idx = i + 1
                break
        
        # Read contaminant data
        if data_start_idx > 0:
            import io
            data_section = ''.join(lines[data_start_idx:])
            reader = csv.DictReader(io.StringIO(data_section))
            
            for row in reader:
                if row.get('Mass Number (A)') and row.get('Charge State (q)'):
                    A = int(row['Mass Number (A)'])
                    q = int(row['Charge State (q)'])
                    isotope = row['Isotope']
                    Aq = float(row['A/q Ratio'])
                    species = row.get('Species', '')
                    
                    Z = get_Z_from_symbol(species)
                    
                    # Check if this is the target beam
                    is_target = (A == target_A and q == target_q)
                    
                    before_data.append({
                        'A': A,
                        'q': q,
                        'Z': Z,  
                        'isotope': isotope,
                        'Aq': Aq,
                        'is_target': is_target
                    })
    
    target_in_list = any(d['is_target'] for d in before_data)
    if not target_in_list and target_A and target_q:
        print("⚠ Target beam not found in contaminants list, adding it...")
        est_Z = int(target_A / 2.5)
        element = get_element_symbol(est_Z)
        before_data.insert(0, {
            'A': target_A,
            'q': target_q,
            'Z': est_Z,
            'isotope': f"{element}-{target_A}",
            'Aq': target_Aq,
            'is_target': True
        })
    
    # Read charge distributions after stripping from CSV 
    after_data = []
    
    with open(charge_dist_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            Ab = int(row['Ab'])
            Zb = int(row['Zb'])
            Qb = int(row['Qb'])  # Initial charge state
            
            # Find corresponding isotope info - match by A, Z, and q
            isotope = None
            is_target = False
            for item in before_data:
                if item['A'] == Ab and item['q'] == Qb and item.get('Z', 0) == Zb:
                    isotope = item['isotope']
                    is_target = item['is_target']
                    break
            
            if isotope is None:
                for item in before_data:
                    if item['A'] == Ab and item['q'] == Qb:
                        isotope = item['isotope']
                        is_target = item['is_target']
                        break
            
            charges = []
            intensities = []
            
            for key, value in row.items():
                if key.startswith('q') and key[1:].isdigit() and value:
                    try:
                        charge = int(key[1:])
                        intensity = float(value)
                        if intensity > 0.01:  # Only keep significant peaks
                            charges.append(charge)
                            intensities.append(intensity)
                    except ValueError:
                        continue
            
            # Normalize intensities. This makes the highest peak = 1
            if intensities:
                max_int = max(intensities)
                intensities = [i / max_int for i in intensities]
            
            if charges:
                after_data.append({
                    'isotope': isotope or f"A={Ab}",
                    'A': Ab,
                    'Z': Zb,
                    'q_before': Qb,
                    'charges': charges,
                    'intensities': intensities,
                    'is_target': is_target
                })
    
    def format_species(isotope, charge):
        parts = isotope.split('-')
        if len(parts) == 2:
            element = parts[0]
            mass = parts[1]
            return f"{mass}{element}{charge}+"
        return f"{isotope}{charge}+"
    
    def format_species_no_charge(isotope):
        parts = isotope.split('-')
        if len(parts) == 2:
            element = parts[0]
            mass = parts[1]
            return f"{mass}{element}"
        return isotope
    
    # Create figure
    fig = plt.figure(figsize=(16, 8))
    # Top panel: Before stripping 
    ax1 = plt.subplot(2, 1, 1)
    ax1.set_title('Before stripping', fontsize=14)
    
    before_data_sorted = sorted(before_data, key=lambda x: x['Aq'])
    
    from matplotlib import cm
    import matplotlib.patches as mpatches
    from collections import defaultdict
    
    n_species = len(before_data_sorted)
    colors = cm.rainbow(np.linspace(0, 1, n_species))
    
    # Create marker types matching the second panel
    marker_types = ['o', '^', 'v', '<', '>', 's', 'D', 'p', '*', 'h']
    species_markers = {}
    species_colors = {}
    
    for i, data in enumerate(before_data_sorted):
        isotope = data['isotope']
        species_markers[isotope] = marker_types[i % len(marker_types)]
        species_colors[isotope] = colors[i]
    
    # Group species by A/q to handle overlaps
    aq_groups = defaultdict(list)
    for i, data in enumerate(before_data_sorted):
        aq_rounded = round(data['Aq'], 6)  
        aq_groups[aq_rounded].append((i, data))
    
    # Plot markers at each A/q position
    legend_handles = []
    plotted_species = set()
    
    for aq_value, species_list in sorted(aq_groups.items()):
        n_species_at_aq = len(species_list)
        
        # For multiple species at same A/q, stack them vertically
        base_height = 0.5
        spacing = 0.12
        
        for idx, (color_idx, data) in enumerate(species_list):
            isotope = data['isotope']
            color = species_colors[isotope]
            marker = species_markers[isotope]
            species_label = format_species_no_charge(isotope)
            
            # Calculate y position: stack multiple species
            if n_species_at_aq == 1:
                y_pos = base_height
            else:
                # Center the stack around base_height
                total_span = (n_species_at_aq - 1) * spacing
                y_pos = base_height - total_span/2 + idx * spacing
            
            # Draw stem from axis to marker
            ax1.plot([aq_value, aq_value], [0, y_pos], 
                    color=color, alpha=0.5, linewidth=2, zorder=1)
            
            # Draw marker
            marker_size = 180 if data['is_target'] else 100
            marker_edge = 2.5 if data['is_target'] else 1
            ax1.scatter(aq_value, y_pos, s=marker_size, c=[color], 
                       marker=marker, edgecolor='black', linewidth=marker_edge, 
                       alpha=0.9, zorder=3)
            
            # Create legend entry (only once per species)
            if isotope not in plotted_species:
                from matplotlib.lines import Line2D
                if data['is_target']:
                    handle = Line2D([0], [0], marker=marker, color='w', 
                                  markerfacecolor=color, markersize=12,
                                  markeredgecolor='black', markeredgewidth=2,
                                  label=f"{species_label} (Primary beam)")
                else:
                    handle = Line2D([0], [0], marker=marker, color='w', 
                                  markerfacecolor=color, markersize=10,
                                  markeredgecolor='black', markeredgewidth=0.5,
                                  label=species_label)
                legend_handles.append(handle)
                plotted_species.add(isotope)
    
    ax1.set_xlabel('A/q', fontsize=12)
    ax1.set_ylabel('Species', fontsize=12)
    ax1.set_ylim([0, 1.0])
    ax1.set_yticks([])  
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    # ax1.spines['left'].set_visible(False)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    ax1.legend(handles=legend_handles, loc='upper left', ncol=5, fontsize=9, 
              framealpha=0.9)
    
    # Bottom panel: After stripping
    ax2 = plt.subplot(2, 1, 2)
    ax2.set_title('After stripping (Maximum peak normalized to 1)', fontsize=14)
    ax2.set_ylabel('Relative Intensity', fontsize=12)
    
    # Use the same colors and markers as the top panel
    legend_handles_2 = []
    plotted_species_2 = set()
    
    for species in after_data:
        isotope = species['isotope']
        is_target = species['is_target']
        color = species_colors.get(isotope, 'gray')
        marker = species_markers.get(isotope, 'o')
        charges = species['charges']
        intensities = species['intensities']
        
        for q, intensity in zip(charges, intensities):
            if intensity < 0.01:  # Skip very weak peaks
                continue
            
            Aq = species['A'] / q
            
            # Create smooth Gaussian around this charge state
            sigma = 0.008
            x_range = np.linspace(Aq - 3*sigma, Aq + 3*sigma, 100)
            gaussian = intensity * np.exp(-0.5 * ((x_range - Aq) / sigma) ** 2)
            
            ax2.plot(x_range, gaussian, color=color, alpha=0.6, linewidth=2)
            ax2.fill_between(x_range, gaussian, alpha=0.2, color=color)
            
            ax2.plot(Aq, intensity, marker=marker, color=color, markersize=8, 
                    markeredgecolor='black', markeredgewidth=1.5)
        
        if isotope not in plotted_species_2:
            label = format_species_no_charge(isotope)
            from matplotlib.lines import Line2D
            if is_target:
                handle = Line2D([0], [0], marker=marker, color='w', 
                              markerfacecolor=color, markersize=10,
                              markeredgecolor='black', markeredgewidth=2,
                              label=f"{label} (Primary beam)")
            else:
                handle = Line2D([0], [0], marker=marker, color='w', 
                              markerfacecolor=color, markersize=10,
                              markeredgecolor='black', markeredgewidth=0.5,
                              label=label)
            legend_handles_2.append(handle)
            plotted_species_2.add(isotope)
    
    ax2.set_xlabel('A/q', fontsize=12)
    ax2.set_ylim([0, 1.4])
    ax2.grid(axis='y', alpha=0.3)
    ax2.legend(handles=legend_handles_2, loc='upper left', ncol=5, fontsize=9, 
              framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(output_image, dpi=150, bbox_inches='tight')
    print(f"Picture saved: {output_image}")
    if target_A and target_q:
        print(f"  Target beam: A={target_A} q={target_q} (A/q={target_Aq:.3f})")
    print(f"  Total species before stripping: {len(before_data)}")
    total_charge_states = sum(len(species['charges']) for species in after_data)
    print(f"  Total charge states after stripping: {total_charge_states}")
    plt.show()


if __name__ == "__main__":
    contaminants_file = "u238_analysis_list.csv"
    charge_dist_file = "u238_analysis_charge_distributions.csv"
    
    if Path(contaminants_file).exists() and Path(charge_dist_file).exists():
        plot_charge_distribution(contaminants_file, charge_dist_file, "charge_states.png")
    else:
        print("Example files not found. Run calculate_and_run_contaminants first!")
        print("\nUsage:")
        print("  from visualize_charge_distribution import plot_charge_distribution")
        print("  plot_charge_distribution('my_analysis_list.csv', 'my_analysis_charge_distributions.csv')")
