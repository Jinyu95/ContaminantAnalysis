"""
Simple beam tracking for contaminant analysis
Reads charge distributions from CSV and tracks through FODO lattice
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pyJuTrack as jt

# Constants
C_LIGHT = 299792458.0  # m/s
AMU_TO_EV = 931.494102e6  # eV/u

# Element names
ELEMENTS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf'
}

def get_element_name(Z):
    """Get element symbol from atomic number"""
    return ELEMENTS.get(Z, f'Z{Z}')


def get_species_name(A, Z):
    """Get species name like 238U or 182W"""
    elem = get_element_name(Z)
    return f"{A}{elem}"


def get_element_type(elem):
    """Classify pyJuTrack element type from its string representation."""
    elem_str = str(elem)
    if 'DRIFT' in elem_str:
        return 'DRIFT'
    if 'KQUAD' in elem_str or 'QUAD' in elem_str:
        return 'KQUAD'
    if 'SBEND' in elem_str or 'RBEND' in elem_str or 'BEND' in elem_str:
        return 'BEND'
    if 'KSEXT' in elem_str or 'SEXT' in elem_str:
        return 'KSEXT'
    if 'KOCT' in elem_str or 'OCT' in elem_str:
        return 'KOCT'
    if 'RFCA' in elem_str or 'CAV' in elem_str:
        return 'RFCA'
    return 'OTHER'


def plot_lattice_layout(ax, lattice, width=0.25):
    """
    Plot lattice layout in a simplified flat view (no bending)
    
    Parameters:
    -----------
    ax : matplotlib axis
        Axis to plot on
    lattice : jt.Lattice
        JuTrack lattice object
    width : float
        Height of elements in plot
    """
    pos = 0.0 
    
    try:
        elements = lattice.line if hasattr(lattice, 'line') else list(lattice)
    except:
        elements = []
    
    for elem in elements:
        L = 0.0
        if hasattr(elem, 'length'):
            L = elem.length
        elif hasattr(elem, 'len'):
            L = elem.len
        
        if L == 0:
            continue
            
        elem_type = get_element_type(elem)
        elem_name = ''
        
        elem_str = str(elem)
        
        if hasattr(elem, 'name'):
            elem_name = elem.name
        elif '"' in elem_str:
            # Extract name from string like 'DRIFT{Float64}("D1", ...'
            start = elem_str.find('"') + 1
            end = elem_str.find('"', start)
            if start > 0 and end > start:
                elem_name = elem_str[start:end]
        
        # Quadrupoles
        if 'QUAD' in elem_type or 'KQUAD' in elem_type:
            k1 = elem.k1 if hasattr(elem, 'k1') else 0.0
            color = 'lightblue'
            
            if k1 >= 0:  # Focusing in x (above beamline)
                rect = Rectangle((pos, 0), L, width, 
                                edgecolor='black', facecolor=color, linewidth=0.8)
            else:  # Defocusing in x (below beamline)
                rect = Rectangle((pos, -width), L, width,
                                edgecolor='black', facecolor=color, linewidth=0.8)
            ax.add_patch(rect)
        
        # Sextupoles
        elif 'SEXT' in elem_type or 'KSEXT' in elem_type:
            color = 'blue'
            rect = Rectangle((pos, -width*0.7/2), L, width*0.7,
                            edgecolor='black', facecolor=color, linewidth=0.8)
            ax.add_patch(rect)
        
        # Octupoles
        elif 'OCT' in elem_type or 'KOCT' in elem_type:
            color = 'green'
            rect = Rectangle((pos, -width*0.4/2), L, width*0.4,
                            edgecolor='black', facecolor=color, linewidth=0.8)
            ax.add_patch(rect)
        
        # Dipoles (SBEND, RBEND, etc.)
        elif 'BEND' in elem_type or 'SBEND' in elem_type:
            color = 'orange'
            rect = Rectangle((pos, -width/2), L, width,
                            edgecolor='black', facecolor=color, linewidth=0.8)
            ax.add_patch(rect)
        
        # RF Cavities
        elif 'RFCA' in elem_type or 'CAV' in elem_type or 'DRIFT' in elem_name.upper():
            if 'CAV' in elem_name.upper():
                color = 'red'
                rect = Rectangle((pos, -width/4), L, width/2,
                                edgecolor='black', facecolor=color, 
                                hatch='//', linewidth=0.8)
                ax.add_patch(rect)
        
        pos += L
    
    # Draw beamline
    total_length = pos if pos > 0 else lattice.total_length()
    ax.plot([0, total_length], [0, 0], 'k-', linewidth=0.8, alpha=0.5)
    
    ax.set_xlim(0, total_length)
    ax.set_ylim(-width*1.5, width*1.5)
    ax.set_aspect('auto')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

def plot_loss_analysis(loss_locations, loss_summary, lattice, output_prefix):
    """
    Create visualizations of particle loss analysis
    
    Parameters:
    -----------
    loss_locations : list of dict
        Detailed loss information at each element
    loss_summary : list of dict
        Summary of transmission for each charge state
    lattice : jt.Lattice
        Lattice object
    output_prefix : str
        Prefix for output files
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    
    fig = plt.subplots(2, 2, figsize=(14, 10))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    ax = axes[0, 0]
    df_summary = pd.DataFrame(loss_summary)
    
    species_list = df_summary['species'].unique()
    
    n_species = len(species_list)
    if n_species <= 10:
        cmap = cm.get_cmap('tab10', n_species)
    else:
        cmap = cm.get_cmap('tab20', n_species)
    species_colors = {sp: cmap(i) for i, sp in enumerate(species_list)}
    
    df_filtered = df_summary[(df_summary['transmission'] < 0.999) | (df_summary['fraction'] > 1.0)]
    
    if len(df_filtered) == 0:
        df_filtered = df_summary.nlargest(min(10, len(df_summary)), 'fraction')
    
    x_positions = list(range(len(df_filtered)))
    labels = [f"{row['species']}{int(row['q'])}+" for _, row in df_filtered.iterrows()]
    colors_list = [species_colors[row['species']] for _, row in df_filtered.iterrows()]
    
    bars = ax.bar(x_positions, df_filtered['transmission']*100, color=colors_list, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=90, ha='center', fontsize=8)
    ax.set_ylabel('Transmission (%)', fontsize=12)
    ax.set_title(f'Transmission by Charge State (showing {len(df_filtered)} with losses or >1% fraction)', 
                 fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim([0, 105])
    
    ax = axes[0, 1]
    ax.bar(x_positions, df_filtered['loss_percent'], color=colors_list, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=90, ha='center', fontsize=8)
    ax.set_ylabel('Loss (%)', fontsize=12)
    ax.set_title(f'Particle Loss by Charge State ({len(df_filtered)} states)', fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    ax = axes[1, 0]
    df_locs = pd.DataFrame(loss_locations)
    
    if len(df_locs) > 0:
        loss_by_s = df_locs.groupby('s_center').agg({
            'n_lost': 'sum',
            'element_type': 'first'
        }).reset_index()
        
        elem_type_colors = {
            'DRIFT': '0.8',
            'KQUAD': 'blue',
            'BEND': 'red',
            'KSEXT': 'green',
            'KOCT': 'purple',
            'RFCA': 'orange',
            'OTHER': 'gray'
        }
        
        colors = [elem_type_colors.get(et, 'gray') for et in loss_by_s['element_type']]
        
        ax.bar(loss_by_s['s_center'], loss_by_s['n_lost'], 
               width=0.05, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.set_xlabel('s (m)', fontsize=12)
        ax.set_ylabel('Particles Lost', fontsize=12)
        ax.set_title('Loss Distribution Along Lattice', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=elem_type_colors[et], label=et, alpha=0.7)
                          for et in set(loss_by_s['element_type'])]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    # Cumulative loss along lattice
    ax = axes[1, 1]
    
    if len(df_locs) > 0:
        # Plot cumulative loss for each species/charge
        # Only plot significant losses (>1% loss) to avoid clutter
        plotted_count = 0
        max_plots = 20  # Limit number of lines to avoid overcrowding
        
        for species_q in df_locs[['species', 'q']].drop_duplicates().values:
            species, q = species_q
            data = df_locs[(df_locs['species'] == species) & (df_locs['q'] == q)]
            data = data.sort_values('s_center')
            
            total_loss_fraction = data['loss_fraction'].sum()
            
            # Only plot if significant loss (>1%) and under max_plots limit
            if total_loss_fraction > 0.01 and plotted_count < max_plots:
                cumulative_loss = np.cumsum(data['loss_fraction'].values) * 100
                s_positions = data['s_center'].values
                
                color = species_colors.get(species, 'gray')
                ax.plot(s_positions, cumulative_loss, 'o-', 
                       label=f"{species}{int(q)}+", color=color, linewidth=2, markersize=3)
                plotted_count += 1
        
        ax.set_xlabel('s (m)', fontsize=12)
        ax.set_ylabel('Cumulative Loss (%)', fontsize=12)
        ax.set_title('Cumulative Loss Along Lattice (>1% loss only)', fontsize=13, fontweight='bold')
        ax.legend(loc='upper left', fontsize=8, ncol=3, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 105])
    
    plt.tight_layout()
    output_file = f"{output_prefix}_loss_analysis.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Loss analysis plot saved: {output_file}")
    plt.close(fig)


def analyze_particle_loss_from_tracking(particles_initial_xy, rout, lattice, aperture_arr, 
                                        spos, species_name, q, fraction):
    """
    Analyze particle losses using already-tracked coordinates from jt.linepass
    
    Parameters:
    -----------
    particles_initial_xy : ndarray
        Initial particle x, px coordinates (N x 2) for tracking particle IDs
    rout : list of ndarray
        List of particle coordinates at each element from jt.linepass
    lattice : jt.Lattice
        Lattice object
    aperture_arr : ndarray
        Aperture at each element (m)
    spos : list
        s-positions of elements
    species_name : str
        Name of species
    q : int
        Charge state
    fraction : float
        Fraction of this charge state
        
    Returns:
    --------
    dict : Loss analysis results
    """
    n_particles = particles_initial_xy.shape[0]
    
    alive_mask = np.ones(n_particles, dtype=bool)
    
    try:
        elements = lattice.line if hasattr(lattice, 'line') else list(lattice)
    except:
        elements = []
    
    loss_info = []
    current_s = 0.0
    
    # Get initial particle coordinates (entrance to first element)
    particles_entrance = particles_initial_xy  # Initial coordinates (x, px)
    
    for elem_idx, elem in enumerate(elements):
        if hasattr(elem, 'length'):
            L = elem.length
        elif hasattr(elem, 'len'):
            L = elem.len
        else:
            L = 0.0
        
        elem_type = get_element_type(elem)
        elem_name = str(elem).split('(')[0] if '(' in str(elem) else str(elem)
        
        s_start = current_s
        s_end = current_s + L
        s_center = (s_start + s_end) / 2.0
        
        # Get aperture at this element
        if elem_idx < len(aperture_arr):
            aperture = aperture_arr[elem_idx]
        else:
            aperture = aperture_arr[-1]
        
        # Check particles at both element entrance and exit
        # rout[elem_idx] corresponds to position at element exit (after element elem_idx)
        
        # Check at entrance (use previous element's exit or initial coordinates)
        if elem_idx == 0:
            # First element: use initial coordinates
            coords_entrance = np.column_stack([particles_entrance[:, 0], 
                                               np.zeros(len(particles_entrance)),
                                               np.zeros(len(particles_entrance)),
                                               np.zeros(len(particles_entrance)),
                                               np.zeros(len(particles_entrance)),
                                               np.zeros(len(particles_entrance))])
        elif elem_idx > 0 and (elem_idx - 1) < len(rout):
            # Use previous element's exit as this element's entrance
            coords_entrance = rout[elem_idx - 1]
        else:
            coords_entrance = None
        
        # Check entrance aperture
        if coords_entrance is not None:
            x_entrance = coords_entrance[alive_mask, 0]
            y_entrance = coords_entrance[alive_mask, 2]
            r_entrance = np.sqrt(x_entrance**2 + y_entrance**2)
            
            # Particles that exceed aperture at entrance
            hit_aperture_entrance = r_entrance > aperture
            n_lost_entrance = np.sum(hit_aperture_entrance)
            
            if n_lost_entrance > 0:
                # Record loss at entrance
                loss_info.append({
                    'element_index': elem_idx,
                    'element_type': elem_type,
                    'element_name': elem_name + '_entrance',
                    's_start': s_start,
                    's_end': s_end,
                    's_center': s_start,  # Entrance position
                    'length': L,
                    'aperture': aperture,
                    'n_lost': n_lost_entrance,
                    'loss_fraction': n_lost_entrance / n_particles
                })
                
                # Update alive mask
                alive_indices = np.where(alive_mask)[0]
                lost_global_indices = alive_indices[hit_aperture_entrance]
                alive_mask[lost_global_indices] = False
        
        # Check at exit
        if elem_idx < len(rout):
            coords_exit = rout[elem_idx]
            
            # Check which alive particles hit the aperture at element exit
            if np.any(alive_mask):
                x_exit = coords_exit[alive_mask, 0]
                y_exit = coords_exit[alive_mask, 2]
                r_exit = np.sqrt(x_exit**2 + y_exit**2)
                
                # Particles that exceed aperture at exit
                hit_aperture_exit = r_exit > aperture
                n_lost_exit = np.sum(hit_aperture_exit)
                
                if n_lost_exit > 0:
                    # Record loss at exit
                    loss_info.append({
                        'element_index': elem_idx,
                        'element_type': elem_type,
                        'element_name': elem_name + '_exit',
                        's_start': s_start,
                        's_end': s_end,
                        's_center': s_end,  # Exit position
                        'length': L,
                        'aperture': aperture,
                        'n_lost': n_lost_exit,
                        'loss_fraction': n_lost_exit / n_particles
                    })
                    
                    # Update alive mask
                    alive_indices = np.where(alive_mask)[0]
                    lost_global_indices = alive_indices[hit_aperture_exit]
                    alive_mask[lost_global_indices] = False
        
        current_s = s_end
    
    n_survived = np.sum(alive_mask)
    transmission = n_survived / n_particles
    
    return {
        'species': species_name,
        'q': q,
        'fraction': fraction,
        'n_initial': n_particles,
        'n_survived': n_survived,
        'transmission': transmission,
        'loss_locations': loss_info
    }


def get_aperture_from_element(elem):
    """
    Extract aperture from element's RApertures attribute.
    Returns the first value in RApertures list, or default if not available.
    
    Parameters:
    -----------
    elem : jt element
        JuTrack element
    
    Returns:
    --------
    float : Aperture radius in meters
    """
    if hasattr(elem, 'RApertures'):
        raper = elem.RApertures
        return float(raper[0])
    return 0.5  # Default 50 cm


def track_contaminants(csv_file, lattice, aperture, main_A, main_q, 
                      output_prefix="tracking", emit_x=1e-6, emit_y=1e-6,
                      beta_twiss_x=10.0, beta_twiss_y=10.0,
                      alpha_twiss_x=0.0, alpha_twiss_y=0.0,
                      x_offset=0.0, y_offset=0.0, plot_envelope=True,
                      s_plot_range_m=1.0, plot_ylim_mm=None, analyze_loss=True,
                      normalized_emittance=True, main_energy_MeV_u=None):
    """
    Track all beams from charge distribution CSV through lattice
    
    Parameters:
    -----------
    csv_file : str
        Charge distribution CSV from ETACHA
    lattice : jt.Lattice
        JuTrack lattice object to track through
    aperture : float or array-like
        Beam pipe radius in meters. Can be:
        - Single float: constant aperture for all elements
        - Array: aperture at each element position (length must match lattice)
    main_A : int
        Mass number of the main beam that accelerator is tuned for
    main_q : int
        Charge state of the main beam that accelerator is tuned for
    output_prefix : str
        Prefix for output files
    emit_x, emit_y : float
        Emittance in m·rad. If normalized_emittance=True, this is normalized rms emittance.
        If normalized_emittance=False, this is geometric emittance.
    beta_twiss_x, beta_twiss_y : float
        Twiss beta functions in meters
    alpha_twiss_x, alpha_twiss_y : float
        Twiss alpha parameters (dimensionless)
    x_offset : float
        Initial horizontal offset (m)
    y_offset : float
        Initial vertical offset (m)
    plot_envelope : bool
        Whether to plot envelope evolution along lattice
    s_plot_range_m : float
        Range of s to plot (meters)
    plot_ylim_mm : float, optional
        Y-axis limit for envelope plots in mm. If None, uses adaptive limits based on data.
    analyze_loss : bool
        Whether to perform detailed particle loss analysis
    normalized_emittance : bool
        If True, emit_x and emit_y are normalized rms emittances (ε_n = beta·gamma·ε_geometric)
        If False, emit_x and emit_y are geometric emittances
    main_energy_MeV_u : float, optional
        Kinetic energy of the main beam in MeV/u. If None, assumes main beam has same energy as contaminants.
    """
    # Read charge distributions calculated by ETACHA4
    df = pd.read_csv(csv_file)
    
    spos = jt.spos(lattice)
    n_elements = len(spos)
    
    # Extract apertures from lattice elements
    # If aperture parameter is provided (backward compatibility), use it
    # Otherwise, extract from element RApertures
    if aperture is None:
        # Extract apertures from lattice elements
        try:
            elements = lattice.line if hasattr(lattice, 'line') else list(lattice)
        except:
            elements = []
        
        aperture_arr = np.array([get_aperture_from_element(elem) for elem in elements])
        print(f"Lattice length: {lattice.total_length():.2f} m")
        print(f"Aperture: from element RApertures, range [{np.min(aperture_arr)*1000:.1f}, "
              f"{np.max(aperture_arr)*1000:.1f}] mm")
    elif isinstance(aperture, (int, float)):
        # Constant aperture
        aperture_arr = np.ones(n_elements) * aperture
        print(f"Lattice length: {lattice.total_length():.2f} m")
        print(f"Aperture: {aperture*1000:.1f} mm (constant)")
    else:
        # Variable aperture array
        aperture_arr = np.array(aperture)
        if len(aperture_arr) != n_elements:
            raise ValueError(f"Aperture array length ({len(aperture_arr)}) must match "
                           f"number of lattice elements ({n_elements})")
        print(f"Lattice length: {lattice.total_length():.2f} m")
        print(f"Aperture: variable, range [{np.min(aperture_arr)*1000:.1f}, "
              f"{np.max(aperture_arr)*1000:.1f}] mm")
    
    print(f"Initial offset: x={x_offset*1000:.1f} mm, y={y_offset*1000:.1f} mm")
    if main_energy_MeV_u is not None:
        print(f"Main beam: A={main_A}, q={main_q}, A/q={main_A/main_q:.4f}, E={main_energy_MeV_u} MeV/u")
    else:
        print(f"Main beam: A={main_A}, q={main_q}, A/q={main_A/main_q:.4f}, E=variable (from CSV)")
    
    results = []
    envelope_data = []
    
    # Calculate main beam Brho (reference for dp/p calculation)
    # Note: We assume main beam has same energy per nucleon as the beams being tracked
    # This will be updated for each species based on their actual energy
    
    # Process each beam
    for idx, row in df.iterrows():
        A = int(row['Ab'])
        Z = int(row['Zb'])
        q_initial = int(row['Qb'])
        E = float(row['EnergyResidual'])                # MeV/u
        abundance = float(row.get('abundance_%', 0))    # Natural abundance percentage
        
        # Find charge states with combined ratios (abundance × charge_fraction)
        charge_states = {}  # {q: (charge_fraction, combined_ratio)}
        for col in df.columns:
            if col.startswith('q') and col[1:].lstrip('-').isdigit() and not col.endswith('_ratio'):
                q_str = col[1:]
                if not q_str.lstrip('-').isdigit():
                    continue
                q = int(q_str)
                charge_frac = row[col]
                
                ratio_col = f"q{q}_ratio"
                if ratio_col in df.columns and pd.notna(row[ratio_col]):
                    combined_ratio = float(row[ratio_col])
                else:
                    # Fallback: calculate it manually if ratio column doesn't exist
                    combined_ratio = abundance * charge_frac / 100.0 if pd.notna(charge_frac) else 0
                
                # Filter: only track if combined ratio >= 1%
                if pd.notna(charge_frac) and charge_frac > 0 and combined_ratio >= 1.0:
                    charge_states[q] = (float(charge_frac), combined_ratio)
        
        if not charge_states:
            continue
        
        species_name = get_species_name(A, Z)
        
        print(f"\n{species_name}: A={A}, Z={Z}, q_initial={q_initial}, E={E} MeV/u, Abundance={abundance:.4f}%")
        
        # Calculate reference Brho based on main beam settings
        # The accelerator is tuned for main_A/main_q at energy E (or main_energy_MeV_u if provided)
        m0_main = main_A * AMU_TO_EV
        
        if main_energy_MeV_u is not None:
            E_main_u = main_energy_MeV_u
        else:
            E_main_u = E
            
        E_kinetic_main = E_main_u * main_A * 1e6  # eV
        E_total_main = E_kinetic_main + m0_main
        gamma_main = E_total_main / m0_main
        beta_main = np.sqrt(1.0 - 1.0/gamma_main**2)
        p_ev_main = beta_main * gamma_main * m0_main  # momentum in eV/c
        brho_main = p_ev_main / (main_q * C_LIGHT * 1.602176634e-19)
        
        # Calculate momentum for this species (same energy per nucleon)
        m0 = A * AMU_TO_EV
        E_kinetic = E * A * 1e6  # eV
        E_total = E_kinetic + m0 
        gamma = E_total / m0
        beta = np.sqrt(1.0 - 1.0/gamma**2)
        p_ev = beta * gamma * m0 
        
        # Track each charge state
        for q, (charge_frac, combined_ratio) in sorted(charge_states.items()):
            # Calculate Brho for this charge state
            brho = p_ev / (q * C_LIGHT * 1.602176634e-19)
            
            # dp/p relative to main beam (not nominal charge state)
            dp_p = (brho - brho_main) / brho_main
            
            # Check if this is the main beam
            is_main_beam = (A == main_A and q == main_q)
            
            if normalized_emittance:
                beta_rel = beta  
                gamma_rel = gamma  
                emit_x_geom = emit_x / (beta_rel * gamma_rel)
                emit_y_geom = emit_y / (beta_rel * gamma_rel)
            else:
                emit_x_geom = emit_x
                emit_y_geom = emit_y
            
            n_particles = 10000
            
            # Calculate Twiss gamma
            gamma_twiss_x = (1 + alpha_twiss_x**2) / beta_twiss_x
            gamma_twiss_y = (1 + alpha_twiss_y**2) / beta_twiss_y
            
            # RMS values using GEOMETRIC emittance
            sigma_x = np.sqrt(emit_x_geom * beta_twiss_x)
            sigma_px = np.sqrt(emit_x_geom * gamma_twiss_x)
            sigma_y = np.sqrt(emit_y_geom * beta_twiss_y)
            sigma_py = np.sqrt(emit_y_geom * gamma_twiss_y)
            
            # Generate uncorrelated Gaussian distributions
            x_uncorr = np.random.randn(n_particles)
            px_uncorr = np.random.randn(n_particles)
            y_uncorr = np.random.randn(n_particles)
            py_uncorr = np.random.randn(n_particles)
            
            # Clip at 2.5 sigma to avoid extreme outliers
            x_uncorr = np.clip(x_uncorr, -2.5, 2.5)
            px_uncorr = np.clip(px_uncorr, -2.5, 2.5)
            y_uncorr = np.clip(y_uncorr, -2.5, 2.5)
            py_uncorr = np.clip(py_uncorr, -2.5, 2.5)
            
            # Apply Twiss correlation: x' = x + α*px (in normalized coordinates)
            # In phase space: px_new = px - α*x/β
            x = x_uncorr * sigma_x + x_offset
            px = (px_uncorr * sigma_px - alpha_twiss_x * x_uncorr * sigma_x / beta_twiss_x)
            
            y = y_uncorr * sigma_y + y_offset
            py = (py_uncorr * sigma_py - alpha_twiss_y * y_uncorr * sigma_y / beta_twiss_y)
            
            z = np.zeros(n_particles)
            
            # Set dp/p
            dp = np.ones(n_particles) * dp_p
            
            particles = np.column_stack([x, px, y, py, z, dp])
            
            # Track with reference points to get envelope along lattice
            beam = jt.Beam(particles.copy(), energy=E_total, mass=m0)
            
            # Get all element positions
            spos = jt.spos(lattice)
            n_elements = len(spos)
            # Julia uses 1-based indexing
            refpts = list(range(1, n_elements + 1))
            
            rout = jt.linepass(lattice, beam, refpts=refpts)
            
            # Perform detailed loss analysis using the tracked coordinates
            loss_result = None
            if analyze_loss:
                loss_result = analyze_particle_loss_from_tracking(
                    particles[:, :2],  # Initial x, px (for particle IDs)
                    rout, lattice, aperture_arr, spos,
                    species_name, q, charge_frac
                )
            
            # Calculate initial envelope and centroid at s=0
            x_initial = particles[:, 0]
            y_initial = particles[:, 2]
            x_rms_initial = np.std(x_initial)
            y_rms_initial = np.std(y_initial)
            x_mean_initial = np.mean(x_initial)
            y_mean_initial = np.mean(y_initial)
            r_rms_initial = np.sqrt(x_rms_initial**2 + y_rms_initial**2)
            clearance_initial = aperture_arr[0] - r_rms_initial
            
            # Initialize arrays with initial values
            x_env_arr = [x_rms_initial]
            y_env_arr = [y_rms_initial]
            x_mean_arr = [x_mean_initial]
            y_mean_arr = [y_mean_initial]
            r_env_arr = [r_rms_initial]
            clearance_arr = [clearance_initial]
            s_arr = [spos[0]]  # s=0
            
            # rout is a list of particle coordinates at each refpt
            # Calculate envelope and centroid at each position after initial
            for i, coords in enumerate(rout):
                # Calculate RMS beam size (envelope) and centroid
                x_arr = coords[:, 0]
                y_arr = coords[:, 2]
                
                # RMS values (standard deviation) and mean (centroid)
                x_rms = np.std(x_arr)
                y_rms = np.std(y_arr)
                x_mean = np.mean(x_arr)
                y_mean = np.mean(y_arr)
                r_rms = np.sqrt(x_rms**2 + y_rms**2)
                
                # Calculate clearance at this position (i+1 because we already added initial)
                clearance_at_i = aperture_arr[i+1] - r_rms
                
                x_env_arr.append(x_rms)
                y_env_arr.append(y_rms)
                x_mean_arr.append(x_mean)
                y_mean_arr.append(y_mean)
                r_env_arr.append(r_rms)
                clearance_arr.append(clearance_at_i)
                s_arr.append(spos[i+1])
            
            # Find minimum clearance and maximum envelope
            min_clearance_idx = np.argmin(clearance_arr)
            min_clearance = clearance_arr[min_clearance_idx]
            min_clearance_s = s_arr[min_clearance_idx]
            
            max_r = r_env_arr[min_clearance_idx]
            max_x = x_env_arr[min_clearance_idx]
            max_y = y_env_arr[min_clearance_idx]
            
            hit = min_clearance < 0
            status = "HIT!" if hit else "OK"
            main_label = " [MAIN BEAM]" if is_main_beam else ""
            
            print(f"  q={q:2d} (charge:{charge_frac:5.2f}%, ratio:{combined_ratio:6.3f}%)")
            
            # Print loss analysis results if available
            if loss_result:
                transmission_pct = loss_result['transmission'] * 100
                n_lost_total = loss_result['n_initial'] - loss_result['n_survived']
                print(f"    Loss analysis: {transmission_pct:.2f}% transmitted, "
                      f"{n_lost_total}/{loss_result['n_initial']} particles lost")
            
            results.append({
                'species': species_name,
                'A': A, 'Z': Z, 'q_initial': q_initial, 'q': q,
                'energy_MeV_u': E, 
                'abundance_%': abundance,  # Natural abundance
                'charge_fraction_%': charge_frac,  # Charge state fraction
                'combined_ratio_%': combined_ratio,  # abundance × charge_fraction
                'brho_Tm': brho, 'dp_p': dp_p,
                'max_x_mm': max_x*1000, 'max_y_mm': max_y*1000, 
                'max_r_mm': max_r*1000,
                'min_clearance_mm': min_clearance*1000, 
                'min_clearance_s_m': min_clearance_s,
                'hit_pipe': hit,
                'transmission': loss_result['transmission'] if loss_result else None,
                'n_particles_lost': (loss_result['n_initial'] - loss_result['n_survived']) if loss_result else None,
                'loss_analysis': loss_result
            })
            
            # Store envelope data for plotting
            if plot_envelope:
                # Extract aperture values at tracked positions only
                aperture_at_tracked = [aperture_arr[i] for i in range(len(s_arr))]
                envelope_data.append({
                    'species': species_name,
                    'A': A, 'Z': Z, 'q': q, 'is_main': is_main_beam, 
                    'abundance': abundance,  # Natural abundance
                    'charge_fraction': charge_frac,  # Charge state fraction
                    'combined_ratio': combined_ratio,  # abundance × charge_fraction
                    's': s_arr, 
                    'x_env': x_env_arr, 'y_env': y_env_arr,
                    'x_mean': x_mean_arr, 'y_mean': y_mean_arr,
                    'aperture': aperture_at_tracked,
                    'loss_result': loss_result  # Add loss analysis results
                })
    
    # Save results
    df_results = pd.DataFrame(results)
    output_csv = f"{output_prefix}_results.csv"
    df_results.to_csv(output_csv, index=False)
    print(f"\nResults saved: {output_csv}")
    
    # Summary
    n_hit = df_results['hit_pipe'].sum()
    print(f"\nSummary: {len(df_results)} charge states tracked, {n_hit} hit pipe")
    
    # Generate loss analysis report and plots if requested
    if analyze_loss:
        print("\n" + "="*60)
        print("PARTICLE LOSS ANALYSIS")
        print("="*60)
        
        # Collect all loss data
        loss_summary = []
        all_loss_locations = []
        
        for result in results:
            if result.get('loss_analysis'):
                loss_data = result['loss_analysis']
                species = loss_data['species']
                q = loss_data['q']
                frac = loss_data['fraction']
                trans = loss_data['transmission']
                
                loss_summary.append({
                    'species': species,
                    'q': q,
                    'fraction': frac,
                    'transmission': trans,
                    'loss_percent': (1 - trans) * 100
                })
                
                # Collect loss locations for this charge state
                for loss_loc in loss_data['loss_locations']:
                    all_loss_locations.append({
                        'species': species,
                        'q': q,
                        'fraction': frac,
                        **loss_loc
                    })
                
                print(f"\n{species}{q}+ (fraction={frac:.1f}%):")
                print(f"  Transmission: {trans*100:.2f}%")
                print(f"  Loss: {(1-trans)*100:.2f}%")
                
                if loss_data['loss_locations']:
                    print(f"  Major loss locations:")
                    # Sort by number of particles lost
                    sorted_losses = sorted(loss_data['loss_locations'], 
                                         key=lambda x: x['n_lost'], reverse=True)
                    for i, loc in enumerate(sorted_losses[:5]):  # Top 5
                        print(f"    {i+1}. {loc['element_type']} @ s={loc['s_center']:.2f}m: "
                              f"{loc['n_lost']} particles ({loc['loss_fraction']*100:.2f}%)")
        
        # Save loss summary to CSV
        if loss_summary:
            df_loss_summary = pd.DataFrame(loss_summary)
            loss_csv = f"{output_prefix}_loss_summary.csv"
            df_loss_summary.to_csv(loss_csv, index=False)
            print(f"\nLoss summary saved: {loss_csv}")
        
        # Save detailed loss locations
        if all_loss_locations:
            df_loss_locs = pd.DataFrame(all_loss_locations)
            loss_detail_csv = f"{output_prefix}_loss_locations.csv"
            df_loss_locs.to_csv(loss_detail_csv, index=False)
            print(f"Loss locations saved: {loss_detail_csv}")
            
            # Create loss visualization
            # plot_loss_analysis(all_loss_locations, loss_summary, lattice, output_prefix)
    
    # Plot envelope evolution along lattice
    if plot_envelope and envelope_data:
        # Extract element apertures once (s_start, s_end, aperture) for plotting
        element_apertures = []
        if lattice is not None:
            try:
                elements = lattice.line if hasattr(lattice, 'line') else list(lattice)
            except:
                elements = []
            
            current_s = 0.0
            for elem in elements:
                if hasattr(elem, 'length'):
                    L = elem.length
                elif hasattr(elem, 'len'):
                    L = elem.len
                else:
                    L = 0.0

                # Skip zero-length elements
                if L == 0.0:
                    continue

                s_start = current_s
                s_end = current_s + L
                
                # Get aperture from element's RApertures (first value)
                elem_aperture = get_aperture_from_element(elem)
                
                element_apertures.append({
                    's_start': s_start,
                    's_end': s_end,
                    'aperture': elem_aperture
                })
                
                current_s = s_end
        
        plot_envelope_evolution(envelope_data, aperture_arr, element_apertures, output_prefix, 
                                lattice=lattice, s_plot_range_m=s_plot_range_m, plot_ylim_mm=plot_ylim_mm)
        plot_all_species_envelope(envelope_data, aperture_arr, element_apertures, output_prefix, 
                                  lattice=lattice, s_plot_range_m=s_plot_range_m, plot_ylim_mm=plot_ylim_mm)
    
    return df_results


def plot_envelope_evolution(envelope_data, aperture_arr, element_apertures, output_prefix, 
                           lattice=None, s_plot_range_m=1.0, plot_ylim_mm=None):
    """
    Plot beam envelope evolution along the lattice
    
    Parameters:
    -----------
    envelope_data : list
        List of dicts with envelope data for each charge state
    aperture_arr : ndarray
        Aperture at each tracked position (m)
    element_apertures : list of dict
        Pre-calculated element apertures with keys: 's_start', 's_end', 'aperture'
    output_prefix : str
        Prefix for output file
    lattice : jt.Lattice, optional
        Lattice object for plotting layout
    s_plot_range_m : float or None
        Range of s to plot (meters). If None, uses full data range without setting xlim.
    plot_ylim_mm : float, optional
        Y-axis limit for envelope plots in mm. If None, uses adaptive limits based on data.
    """
    # Find the main beam data (to plot on all species plots)
    main_beam_data = None
    for data in envelope_data:
        if data['is_main']:
            main_beam_data = data
            break
    
    # Group by species
    beams = {}
    for data in envelope_data:
        species = data['species']
        if species not in beams:
            beams[species] = []
        beams[species].append(data)
    
    # Create plots for each species
    for species, charge_states in beams.items():
        # Create figure with 4 rows if lattice provided, else 3 rows (X, Y, Loss)
        if lattice:
            fig, axes = plt.subplots(4, 1, figsize=(16, 9), 
                                    gridspec_kw={'height_ratios': [3, 3, 2, 0.8], 'hspace': 0.3})
            ax_x = axes[0]
            ax_y = axes[1]
            ax_loss = axes[2]
            ax_lattice = axes[3]
        else:
            fig, axes = plt.subplots(3, 1, figsize=(16, 8),
                                    gridspec_kw={'height_ratios': [3, 3, 2], 'hspace': 0.3})
            ax_x = axes[0]
            ax_y = axes[1]
            ax_loss = axes[2]
        
        # Sort by charge state
        charge_states = sorted(charge_states, key=lambda x: x['q'])

        # Get aperture array from first charge state (all have same aperture)
        s_arr = np.array(charge_states[0]['s'])
        
        # Determine if aperture is constant or variable
        is_constant_aperture = len(np.unique(aperture_arr)) == 1
        
        # Collect all envelope data to determine adaptive y-limits
        all_x_upper = []
        all_x_lower = []
        all_y_upper = []
        all_y_lower = []
        
        for cs in charge_states:
            x_upper = (np.array(cs['x_mean']) + np.array(cs['x_env'])) * 1000
            x_lower = (np.array(cs['x_mean']) - np.array(cs['x_env'])) * 1000
            y_upper = (np.array(cs['y_mean']) + np.array(cs['y_env'])) * 1000
            y_lower = (np.array(cs['y_mean']) - np.array(cs['y_env'])) * 1000
            all_x_upper.extend(x_upper)
            all_x_lower.extend(x_lower)
            all_y_upper.extend(y_upper)
            all_y_lower.extend(y_lower)
        
        # Add main beam data if different species
        if main_beam_data and main_beam_data['species'] != species:
            x_upper = (np.array(main_beam_data['x_mean']) + np.array(main_beam_data['x_env'])) * 1000
            x_lower = (np.array(main_beam_data['x_mean']) - np.array(main_beam_data['x_env'])) * 1000
            y_upper = (np.array(main_beam_data['y_mean']) + np.array(main_beam_data['y_env'])) * 1000
            y_lower = (np.array(main_beam_data['y_mean']) - np.array(main_beam_data['y_env'])) * 1000
            all_x_upper.extend(x_upper)
            all_x_lower.extend(x_lower)
            all_y_upper.extend(y_upper)
            all_y_lower.extend(y_lower)
        
        # Calculate adaptive limits with margin
        x_max_envelope = max(abs(max(all_x_upper)), abs(min(all_x_lower)))
        y_max_envelope = max(abs(max(all_y_upper)), abs(min(all_y_lower)))
        
        # Use aperture as reference, but ensure envelope is visible
        if is_constant_aperture:
            aperture_mm = aperture_arr[0] * 1000
        else:
            aperture_mm = np.max(aperture_arr) * 1000
        
        # Use user-specified limits if provided, otherwise adaptive
        if plot_ylim_mm is not None:
            x_plot_lim = plot_ylim_mm
            y_plot_lim = plot_ylim_mm
        else:
            x_plot_lim = max(x_max_envelope * 1.3, aperture_mm * 1.1)
            y_plot_lim = max(y_max_envelope * 1.3, aperture_mm * 1.1)
        
        # Plot 1: X envelope
        ax = ax_x
        
        # First plot the main beam (if not in current species)
        if main_beam_data and main_beam_data['species'] != species:
            main_species = main_beam_data['species']
            main_q = main_beam_data['q']
            label_env = f"{main_species}{main_q}+ (MAIN)"
            # Plot centroid with label for legend
            ax.plot(main_beam_data['s'], np.array(main_beam_data['x_mean'])*1000, 
                   linewidth=2, color='k', linestyle='-', zorder=10, alpha=0.8, label=label_env)
            # Plot envelope (centroid ± σ) without label
            x_upper = (np.array(main_beam_data['x_mean']) + np.array(main_beam_data['x_env'])) * 1000
            x_lower = (np.array(main_beam_data['x_mean']) - np.array(main_beam_data['x_env'])) * 1000
            ax.fill_between(main_beam_data['s'], x_lower, x_upper, 
                           color='k', alpha=0.2, zorder=10)
        
        # Then plot this species' charge states
        for idx, cs in enumerate(charge_states):
            q = cs['q']
            is_main = cs['is_main']
            abundance = cs['abundance']
            combined_ratio = cs['combined_ratio']
            charge_fraction = cs['charge_fraction']
            # Assign consistent color for this charge state
            color = plt.cm.tab10(idx % 10)
            
            if is_main:
                # Main beam from this species
                label = f"{species}{q}+ (MAIN, {charge_fraction:.3f}%)"
                # Plot centroid with label
                ax.plot(cs['s'], np.array(cs['x_mean'])*1000, 
                       linewidth=2, color='k', zorder=10, alpha=0.8, label=label)
                # Plot envelope (centroid ± σ) without label
                x_upper = (np.array(cs['x_mean']) + np.array(cs['x_env'])) * 1000
                x_lower = (np.array(cs['x_mean']) - np.array(cs['x_env'])) * 1000
                ax.fill_between(cs['s'], x_lower, x_upper, 
                               color='k', alpha=0.2, zorder=10)
            else:
                # Contaminants
                label = f"{species}{q}+ ({charge_fraction:.3f}%)"
                # Plot centroid with label
                ax.plot(cs['s'], np.array(cs['x_mean'])*1000, 
                       linewidth=1.2, alpha=0.7, label=label, color=color)
                # Plot envelope (centroid ± σ) without label
                x_upper = (np.array(cs['x_mean']) + np.array(cs['x_env'])) * 1000
                x_lower = (np.array(cs['x_mean']) - np.array(cs['x_env'])) * 1000
                ax.fill_between(cs['s'], x_lower, x_upper, 
                               alpha=0.15, color=color)
        
        # Aperture shadow for X (outside pipe, element-by-element)
        # Filter elements to visible range if s_plot_range_m is set
        visible_elements = element_apertures
        if s_plot_range_m is not None:
            visible_elements = [e for e in element_apertures if e['s_start'] < s_plot_range_m]
        
        for elem_ap in visible_elements:
            s_start = elem_ap['s_start']
            s_end = elem_ap['s_end']
            ap_mm = elem_ap['aperture'] * 1000.0
            
            # Shadow above aperture
            ax.fill_between(
                [s_start, s_end],
                ap_mm, x_plot_lim,
                color='0.3',
                alpha=0.8,
                edgecolor=None,
                linewidth=0,
                zorder=-5
            )
            
            # Shadow below aperture
            ax.fill_between(
                [s_start, s_end],
                -x_plot_lim, -ap_mm,
                color='0.3',
                alpha=0.8,
                edgecolor=None,
                linewidth=0,
                zorder=-5
            )
            
            # Draw aperture boundary lines
            ax.plot([s_start, s_end], [ap_mm, ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')
            ax.plot([s_start, s_end], [-ap_mm, -ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')

        # ax.set_xlabel('s (m)', fontsize=13)
        ax.set_ylabel('X (mm)', fontsize=13)
        # Get abundance from first charge state (all have same abundance for this species)
        species_abundance = charge_states[0]['abundance']
        ax.set_title(f'{species} (Abundance: {species_abundance:.4f}%)', fontsize=14, fontweight='bold')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=11, ncol=1, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        if s_plot_range_m is not None:
            ax.set_xlim([0, s_plot_range_m])
        ax.set_ylim([-x_plot_lim, x_plot_lim])
        
        # Plot 2: Y envelope
        ax = ax_y
        
        # First plot the main beam (if not in current species)
        if main_beam_data and main_beam_data['species'] != species:
            main_species = main_beam_data['species']
            main_q = main_beam_data['q']
            label_env = f"{main_species}{main_q}+ (MAIN)"
            # Plot centroid with label
            ax.plot(main_beam_data['s'], np.array(main_beam_data['y_mean'])*1000, 
                   linewidth=2, color='k', linestyle='-', zorder=10, alpha=0.8, label=label_env)
            # Plot envelope (centroid ± σ) without label
            y_upper = (np.array(main_beam_data['y_mean']) + np.array(main_beam_data['y_env'])) * 1000
            y_lower = (np.array(main_beam_data['y_mean']) - np.array(main_beam_data['y_env'])) * 1000
            ax.fill_between(main_beam_data['s'], y_lower, y_upper, 
                           color='k', alpha=0.2, zorder=10)
        
        # Then plot this species' charge states
        for idx, cs in enumerate(charge_states):
            q = cs['q']
            is_main = cs['is_main']
            combined_ratio = cs['combined_ratio']
            charge_fraction = cs['charge_fraction']
            # Assign consistent color for this charge state (same as X plot)
            color = plt.cm.tab10(idx % 10)
            
            if is_main:
                # Main beam from this species
                label = f"{species}{q}+ (MAIN, {charge_fraction:.3f}%)"
                # Plot centroid with label
                ax.plot(cs['s'], np.array(cs['y_mean'])*1000, 
                       linewidth=2, color='black', zorder=10, alpha=0.8, label=label)
                # Plot envelope (centroid ± σ) without label
                y_upper = (np.array(cs['y_mean']) + np.array(cs['y_env'])) * 1000
                y_lower = (np.array(cs['y_mean']) - np.array(cs['y_env'])) * 1000
                ax.fill_between(cs['s'], y_lower, y_upper, 
                               color='black', alpha=0.2, zorder=10)
            else:
                # Contaminants
                label = f"{species}{q}+ ({charge_fraction:.3f}%)"
                # Plot centroid with label
                ax.plot(cs['s'], np.array(cs['y_mean'])*1000, 
                       linewidth=1.2, alpha=0.7, label=label, color=color)
                # Plot envelope (centroid ± σ) without label
                y_upper = (np.array(cs['y_mean']) + np.array(cs['y_env'])) * 1000
                y_lower = (np.array(cs['y_mean']) - np.array(cs['y_env'])) * 1000
                ax.fill_between(cs['s'], y_lower, y_upper, 
                               alpha=0.15, color=color)
        
        # Aperture shadow for Y (element-by-element)
        # Filter elements to visible range if s_plot_range_m is set
        visible_elements = element_apertures
        if s_plot_range_m is not None:
            visible_elements = [e for e in element_apertures if e['s_start'] < s_plot_range_m]
        
        for elem_ap in visible_elements:
            s_start = elem_ap['s_start']
            s_end = elem_ap['s_end']
            ap_mm = elem_ap['aperture'] * 1000.0
            
            # Shadow above aperture
            ax.fill_between(
                [s_start, s_end],
                ap_mm, y_plot_lim,
                color='0.3',
                alpha=0.8,
                edgecolor=None,
                linewidth=0,
                zorder=-5
            )
            
            # Shadow below aperture
            ax.fill_between(
                [s_start, s_end],
                -y_plot_lim, -ap_mm,
                color='0.3',
                alpha=0.8,
                edgecolor=None,
                linewidth=0,
                zorder=-5
            )
            
            # Draw aperture boundary lines
            ax.plot([s_start, s_end], [ap_mm, ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')
            ax.plot([s_start, s_end], [-ap_mm, -ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')

        # ax.set_xlabel('s (m)', fontsize=13)
        ax.set_ylabel('Y (mm)', fontsize=13)
      #   ax.set_title('Vertical Envelope', fontsize=14, fontweight='bold')
        # ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=9, ncol=1, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        if s_plot_range_m is not None:
            ax.set_xlim([0, s_plot_range_m])
        ax.set_ylim([-y_plot_lim, y_plot_lim])
        
        # Plot cumulative particle loss for each charge state
        ax = ax_loss
        for idx, data in enumerate(charge_states):
            q = data['q']
            combined_ratio = data['combined_ratio']
            loss_result = data.get('loss_result')
            transmission = loss_result['transmission'] if loss_result else 1.0
            # Use same color as envelope plots
            is_main = data['is_main']
            if is_main:
                color = 'black'
            else:
                color = plt.cm.tab10(idx % 10)
            
            if loss_result and loss_result['loss_locations']:
                # Build cumulative loss curve from loss_locations
                loss_locations = loss_result['loss_locations']
                n_initial = loss_result['n_initial']
                
                # Sort loss locations by s position (use s_center for plotting)
                loss_s = [loc['s_center'] for loc in loss_locations]
                loss_counts = [loc['n_lost'] for loc in loss_locations]
                
                # Create cumulative loss array
                if loss_s:
                    # Combine s positions and losses
                    s_loss_pairs = sorted(zip(loss_s, loss_counts))
                    s_sorted = [0] + [s for s, _ in s_loss_pairs]
                    cumulative_loss_pct = [0]
                    total_lost = 0
                    
                    for s, count in s_loss_pairs:
                        total_lost += count
                        cumulative_loss_pct.append(total_lost / n_initial * 100)
                    
                    # Extend line to end of beamline (last tracked position)
                    s_max = np.max(s_arr)
                    if s_sorted[-1] < s_max:
                        s_sorted.append(s_max)
                        cumulative_loss_pct.append(cumulative_loss_pct[-1])  # Maintain final loss value
                    
                    # Plot cumulative loss as step function
                    label = f"{species}{q}+, transmission:{transmission*100:.2f}%"
                    ax.plot(s_sorted, cumulative_loss_pct, '-', color=color, 
                           label=label, linewidth=2.5, drawstyle='steps-post')
                else:
                    # No losses - plot flat line at 0
                    label = f"{species}{q}+, transmission:{transmission*100:.2f}%"
                    ax.plot([0, s_plot_range_m if s_plot_range_m else np.max(s_arr)], [0, 0], '-', color=color, 
                           label=label, linewidth=2, alpha=0.5)
            else:
                # No losses - plot flat line at 0
                label = f"{species}{q}+, transmission:{transmission*100:.2f}%"
                ax.plot([0, s_plot_range_m if s_plot_range_m else np.max(s_arr)], [0, 0], '-', color=color,
                       label=label, linewidth=2, alpha=0.5)
        
        # ax.set_xlabel('s (m)', fontsize=13)
        ax.set_ylabel('Cumulative Loss (%)', fontsize=13)
        if s_plot_range_m is not None:
            ax.set_xlim([0, s_plot_range_m])
        ax.set_ylim([0, 105])
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.01), fontsize=11, framealpha=0.9)
        
        # Plot lattice layout if provided (at the bottom)
        if lattice:
            plot_lattice_layout(ax_lattice, lattice, width=0.25)
            if s_plot_range_m is not None:
                ax_lattice.set_xlim(0, s_plot_range_m)
            ax_lattice.set_xlabel('s (m)', fontsize=13)
        
        plt.tight_layout()
        # Clean filename - replace special characters
        species_clean = species.replace('-', '')
        output_file = f"{output_prefix}_envelope_{species_clean}.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Envelope plot saved: {output_file}")
        plt.close()


def plot_all_species_envelope(envelope_data, aperture_arr, element_apertures, output_prefix, 
                             lattice=None, s_plot_range_m=1.0, plot_ylim_mm=None):
    """
    Plot all species and charge states together in a single comprehensive plot
    
    Parameters:
    -----------
    envelope_data : list
        List of dicts with envelope data for each charge state
    aperture_arr : ndarray
        Aperture at each tracked position (m)
    element_apertures : list of dict
        Pre-calculated element apertures with keys: 's_start', 's_end', 'aperture'
    output_prefix : str
        Prefix for output file
    lattice : jt.Lattice, optional
        Lattice object for plotting layout
    s_plot_range_m : float or None
        Range of s to plot (meters). If None, uses full data range without setting xlim.
    plot_ylim_mm : float, optional
        Y-axis limit for envelope plots in mm. If None, uses adaptive limits based on data.
    """
    import matplotlib.cm as cm
    from itertools import cycle
    
    # Create figure with 3 rows if lattice provided, else 2 rows
    if lattice:
        fig, axes = plt.subplots(3, 1, figsize=(13, 7), 
                                gridspec_kw={'height_ratios': [4, 4, 0.8], 'hspace': 0.3})
        ax_x = axes[0]
        ax_y = axes[1]
        ax_lattice = axes[2]
    else:
        fig, axes = plt.subplots(2, 1, figsize=(13, 5))
        ax_x = axes[0]
        ax_y = axes[1]
    
    # Get aperture array from first charge state
    s_arr = np.array(envelope_data[0]['s'])
    
    # Determine if aperture is constant or variable
    is_constant_aperture = len(np.unique(aperture_arr)) == 1
    
    # Group by species to assign colors
    species_dict = {}
    for data in envelope_data:
        species = data['species']
        if species not in species_dict:
            species_dict[species] = []
        species_dict[species].append(data)
    
    # Generate color palette for species
    n_species = len(species_dict)
    colors = cm.tab10(np.linspace(0, 1, min(10, n_species)))
    if n_species > 10:
        colors = cm.tab20(np.linspace(0, 1, n_species))
    
    color_map = {species: colors[i] for i, species in enumerate(species_dict.keys())}
    
    # Track which species have been added to legend
    species_in_legend = set()
    
    # Collect all envelope data for adaptive y-limits
    all_x_upper = []
    all_x_lower = []
    all_y_upper = []
    all_y_lower = []
    
    for data in envelope_data:
        x_upper = (np.array(data['x_mean']) + np.array(data['x_env'])) * 1000
        x_lower = (np.array(data['x_mean']) - np.array(data['x_env'])) * 1000
        y_upper = (np.array(data['y_mean']) + np.array(data['y_env'])) * 1000
        y_lower = (np.array(data['y_mean']) - np.array(data['y_env'])) * 1000
        all_x_upper.extend(x_upper)
        all_x_lower.extend(x_lower)
        all_y_upper.extend(y_upper)
        all_y_lower.extend(y_lower)
    
    # Calculate adaptive limits with margin
    x_max_envelope = max(abs(max(all_x_upper)), abs(min(all_x_lower)))
    y_max_envelope = max(abs(max(all_y_upper)), abs(min(all_y_lower)))
    
    # Use aperture as reference, but ensure envelope is visible
    if is_constant_aperture:
        aperture_mm = aperture_arr[0] * 1000
    else:
        aperture_mm = np.max(aperture_arr) * 1000
    
    # Use user-specified limits if provided, otherwise adaptive
    if plot_ylim_mm is not None:
        x_plot_lim = plot_ylim_mm
        y_plot_lim = plot_ylim_mm
    else:
        x_plot_lim = max(x_max_envelope * 1.3, aperture_mm * 1.1)
        y_plot_lim = max(y_max_envelope * 1.3, aperture_mm * 1.1)
    
    # Plot 1: X envelope
    ax = ax_x
    
    for data in envelope_data:
        species = data['species']
        q = data['q']
        is_main = data['is_main']
        abundance = data['abundance']
        combined_ratio = data['combined_ratio']
        
        color = color_map[species]
        
        # Main beam gets thick line, contaminants get thin line
        if is_main:
            centroid_linewidth = 2.5
            centroid_alpha = 0.9
            envelope_alpha = 0.25
            color = 'k'  # Main beam in black
            label = f"{species}{q}+ (MAIN, abund:{abundance:.2f}%)"
        else:
            centroid_linewidth = 1.0
            centroid_alpha = 0.6
            envelope_alpha = 0.1
            # Only add species to legend once, for the first charge state
            if species in species_in_legend:
                label = None  # Don't add to legend
            else:
                label = f"{species} (abund:{abundance:.2f}%)"
                species_in_legend.add(species)
        
        # Plot centroid with label for legend
        ax.plot(data['s'], np.array(data['x_mean'])*1000, 
               color=color, linewidth=centroid_linewidth, alpha=centroid_alpha, label=label)
        
        # Plot envelope (centroid ± σ) without label
        x_upper = (np.array(data['x_mean']) + np.array(data['x_env'])) * 1000
        x_lower = (np.array(data['x_mean']) - np.array(data['x_env'])) * 1000
        ax.fill_between(data['s'], x_lower, x_upper, 
                       color=color, alpha=envelope_alpha)
    
    # Aperture shadow for X (outside pipe, element-by-element)
    # Filter elements to visible range if s_plot_range_m is set
    visible_elements = element_apertures
    if s_plot_range_m is not None:
        visible_elements = [e for e in element_apertures if e['s_start'] < s_plot_range_m]
    
    for elem_ap in visible_elements:
        s_start = elem_ap['s_start']
        s_end = elem_ap['s_end']
        ap_mm = elem_ap['aperture'] * 1000.0
        
        # Shadow above aperture
        ax.fill_between(
            [s_start, s_end],
            ap_mm, x_plot_lim,
            color='0.3',
            alpha=0.8,
            edgecolor=None,
            linewidth=0,
            zorder=-5
        )
        
        # Shadow below aperture
        ax.fill_between(
            [s_start, s_end],
            -x_plot_lim, -ap_mm,
            color='0.3',
            alpha=0.8,
            edgecolor=None,
            linewidth=0,
            zorder=-5
        )
        
        # Draw aperture boundary lines
        ax.plot([s_start, s_end], [ap_mm, ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')
        ax.plot([s_start, s_end], [-ap_mm, -ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')
    
    # ax.set_xlabel('s (m)', fontsize=13)
    ax.set_ylabel('X (mm)', fontsize=13)
    
    # Smart legend placement: keep legend compact and above plot to avoid covering data
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), fontsize=9,
              ncol=min(6, len(ax.get_legend_handles_labels()[0])), framealpha=0.9)
    
    ax.grid(True, alpha=0.3)
    if s_plot_range_m is not None:
        ax.set_xlim([0, s_plot_range_m])
    ax.set_ylim([-x_plot_lim, x_plot_lim])
    
    # Reset for Y plot
    species_in_legend = set()
    
    # Plot 2: Y envelope
    ax = ax_y
    
    for data in envelope_data:
        species = data['species']
        q = data['q']
        is_main = data['is_main']
        abundance = data['abundance']
        
        color = color_map[species]
        
        # Main beam gets thick line, contaminants get thin line
        if is_main:
            centroid_linewidth = 2.5
            centroid_alpha = 0.9
            envelope_alpha = 0.25
            color = 'k'  # Main beam in black
            label = f"{species}{q}+ (MAIN, abund:{abundance:.2f}%)"
        else:
            centroid_linewidth = 1.0
            centroid_alpha = 0.6
            envelope_alpha = 0.1
            # Only add species to legend once
            if species in species_in_legend:
                label = None
            else:
                label = f"{species} (abund:{abundance:.2f}%)"
                species_in_legend.add(species)
        
        # Plot centroid with label for legend
        ax.plot(data['s'], np.array(data['y_mean'])*1000, 
               color=color, linewidth=centroid_linewidth, alpha=centroid_alpha, label=label)
        
        # Plot envelope (centroid ± σ) without label
        y_upper = (np.array(data['y_mean']) + np.array(data['y_env'])) * 1000
        y_lower = (np.array(data['y_mean']) - np.array(data['y_env'])) * 1000
        ax.fill_between(data['s'], y_lower, y_upper, 
                       color=color, alpha=envelope_alpha)
    
    # Aperture shadow for Y (element-by-element)
    # Filter elements to visible range if s_plot_range_m is set
    visible_elements = element_apertures
    if s_plot_range_m is not None:
        visible_elements = [e for e in element_apertures if e['s_start'] < s_plot_range_m]
    
    for elem_ap in visible_elements:
        s_start = elem_ap['s_start']
        s_end = elem_ap['s_end']
        ap_mm = elem_ap['aperture'] * 1000.0
        
        # Shadow above aperture
        ax.fill_between(
            [s_start, s_end],
            ap_mm, y_plot_lim,
            color='0.3',
            alpha=0.8,
            edgecolor=None,
            linewidth=0,
            zorder=-5
        )
        
        # Shadow below aperture
        ax.fill_between(
            [s_start, s_end],
            -y_plot_lim, -ap_mm,
            color='0.3',
            alpha=0.8,
            edgecolor=None,
            linewidth=0,
            zorder=-5
        )
        
        # Draw aperture boundary lines
        ax.plot([s_start, s_end], [ap_mm, ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')
        ax.plot([s_start, s_end], [-ap_mm, -ap_mm], color='black', linewidth=2.5, zorder=10, alpha=1.0, solid_capstyle='butt')
    
    # ax.set_xlabel('s (m)', fontsize=13)
    ax.set_ylabel('Y (mm)', fontsize=13)
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), fontsize=9,
    #           ncol=min(6, len(ax.get_legend_handles_labels()[0])), framealpha=0.9)
    ax.grid(True, alpha=0.3)
    if s_plot_range_m is not None:
        ax.set_xlim([0, s_plot_range_m])
    ax.set_ylim([-y_plot_lim, y_plot_lim])
    
    # Plot lattice layout if provided (at the bottom)
    if lattice:
        plot_lattice_layout(ax_lattice, lattice, width=0.25)
        if s_plot_range_m is not None:
            ax_lattice.set_xlim(0, s_plot_range_m)
        ax_lattice.set_xlabel('s (m)', fontsize=13)
    
    plt.tight_layout()
    
    output_file = f"{output_prefix}_envelope_all_species.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"All species envelope plot saved: {output_file}")
    plt.close()


if __name__ == "__main__":
    # Define main beam that accelerator is tuned for
    main_A = 238  
    main_q = 75
    main_energy_MeV = 16.55  # MeV/u
    
    # Create lattice
    import stripper2CSS_lattice
    lattice = stripper2CSS_lattice.create_stripper2CSS_lattice()
    
    print(f"  Lattice length: {lattice.total_length():.2f} m")
    
    # Beam parameters
    emit_n_x = 0.15 * 1e-6  # m·rad (normalized rms emittance)
    emit_n_y = 0.15 * 1e-6  # m·rad (normalized rms emittance)
    beta_twiss_x = 0.2  # m
    beta_twiss_y = 0.2  # m
    alpha_twiss_x = 0.0  
    alpha_twiss_y = 0.0  
    
    print(f"\nBeam parameters:")
    print(f"  Normalized rms emittance: εₙ = {emit_n_x*1e6/np.pi:.2f}π mm·mrad")
    print(f"  Beta functions: βx = {beta_twiss_x} m, βy = {beta_twiss_y} m")
    print(f"  Alphax = {alpha_twiss_x}, Alphay = {alpha_twiss_y}")
    
    # Track contaminants through the lattice
    csv_file = "U238_analysis_charge_distributions.csv" 

    # make dir
    import os
    output_dir = "U238_fs1d_tracking_results"
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, "U238_fs1d_tracking")

    track_contaminants(csv_file, 
                      lattice=lattice,
                      aperture=None,  # Apertures from element RApertures
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
                      plot_ylim_mm=50.0,  
                      main_energy_MeV_u=main_energy_MeV)  