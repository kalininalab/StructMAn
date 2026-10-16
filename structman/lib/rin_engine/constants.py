"""Codes and thresholds of the RIN engine - the single place where its criteria are defined.

ENGINE_VERSION goes into every RIN file. Change it whenever a threshold, a typing rule or the
meaning of a stored value changes: the RIN DB builder recomputes every entry written by another
version.
"""

ENGINE_VERSION = '1.0.0'

# ---------------------------------------------------------------------------------------------
# Interaction families (edge types). The values are stored in the RIN files: append, never renumber.

CONTACT = 0          # Voronota-LT contact area > 0; score = area in A^2
HBOND = 1            # D-H...A
WEAK_HBOND = 2       # C-H...O
POLAR = 3            # donor-acceptor proximity without angle terms
WEAK_POLAR = 4       # C(H)-O proximity without angle terms
IONIC = 5            # cationic - anionic atoms
PIPI_PARALLEL = 6    # aromatic ring stacking, parallel (face to face, displaced)
PIPI_T = 7           # aromatic ring stacking, perpendicular (T-shaped, edge to face)
PICATION = 8         # aromatic ring - cation
PIHYDROGEN = 9       # X-H...aromatic ring
DISULFIDE = 10       # S-S
METAL = 11           # metal - N/O/S coordination
HALOGEN = 12         # C-X...A and C-X...aromatic ring (X = Cl, Br, I)
HYDROPHOBIC = 13     # hydrophobic atom - hydrophobic atom

FAMILY_NAMES = (
    'contact', 'hbond', 'weak_hbond', 'polar', 'weak_polar', 'ionic', 'pipi_parallel', 'pipi_t',
    'pication', 'pihydrogen', 'disulfide', 'metal', 'halogen', 'hydrophobic',
)

# ---------------------------------------------------------------------------------------------
# Atom classes, the side of an edge an interaction uses (legacy RINerator: mc / sc / ligand).

MAIN_CHAIN = 0       # amino-acid backbone: N, CA, C, O, OXT and their hydrogens
SIDE_CHAIN = 1       # the rest of an amino acid
NUCLEIC_BACKBONE = 2 # phosphate and sugar of a nucleotide
NUCLEIC_BASE = 3     # the base of a nucleotide
OTHER = 4            # ligands, metals, ions, sugars, ...

CLASS_NAMES = ('mc', 'sc', 'nbb', 'nba', 'lig')

AMINO_ACID_BACKBONE = frozenset(('N', 'CA', 'C', 'O', 'OXT'))
NUCLEOTIDE_BACKBONE = frozenset((
    'P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P',
    "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
))

# ---------------------------------------------------------------------------------------------
# Node kinds (residues). Partner categories of the runtime profiles (ligand, metal, ion, DNA, RNA,
# chain types) are derived from these and from the entity types at runtime.

KIND_PROTEIN = 0     # residue of a polypeptide entity, modified residues included
KIND_NUCLEOTIDE = 1  # residue of a DNA / RNA / hybrid entity
KIND_LIGAND = 2      # non-polymer with more than one heavy atom
KIND_METAL = 3       # single metal atom (ZN, MG, FE, ...)
KIND_ION = 4         # single non-metal atom (CL, BR, IOD, ...)
KIND_SUGAR = 5       # residue of a branched (oligosaccharide) entity
KIND_OTHER = 6       # residue of any other polymer type

KIND_NAMES = ('protein', 'nucleotide', 'ligand', 'metal', 'ion', 'sugar', 'other')

# Upper case, as biotite writes elements.
METAL_ELEMENTS = frozenset((
    'LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU',
    'ZN', 'GA', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',
    'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
    'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'FR', 'RA',
    'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM',
))
HALOGEN_DONOR_ELEMENTS = frozenset(('CL', 'BR', 'I'))
WATER_NAMES = frozenset(('HOH', 'DOD', 'WAT', 'H2O'))

# The 20 amino acids whose protonation hydride estimates; modified residues keep CCD charges.
CANONICAL_AMINO_ACIDS = frozenset((
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
))

# Charged groups of the canonical amino acids, by atom name. Both carboxylate oxygens and all
# guanidinium nitrogens count, whichever atom carries the formal charge. Histidine is neutral at
# pH 7 and only counts when its estimated charge is positive.
AMINO_ACID_CATIONIC = {'LYS': ('NZ',), 'ARG': ('NE', 'NH1', 'NH2')}
AMINO_ACID_ANIONIC = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}

# Bondi van der Waals radii (A) of halogen bond partners.
BONDI_RADII = {'O': 1.52, 'N': 1.55, 'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'I': 1.98}

# ---------------------------------------------------------------------------------------------
# Thresholds. Distances in A, angles in degrees.

PH = 7.0

THRESHOLDS = {
    # Atoms within this distance of a kept chain take part; everything else is cut away first.
    # Above the largest interaction distance (pi-pi 6.5).
    'relevant_radius': 7.0,
    # Atom pairs closer in the bond graph than this are not interactions (1-2, 1-3, 1-4).
    'max_excluded_bond_separation': 3,
    # Voronota-LT probe radius.
    'probe_radius': 1.4,
    # Hydrogen bond, Baker & Hubbard criteria as in biotite.structure.hbond; RING strict uses
    # D...A <= 3.5, which these imply for linear bonds.
    'hbond_h_acceptor': 2.5,
    'hbond_angle': 120.0,
    # Weak hydrogen bond C-H...O, Arpeggio (arpeggio/core/config.py CONTACT_TYPES['weak hbond']).
    'weak_hbond_donor_acceptor': 3.6,
    'weak_hbond_angle': 130.0,
    # Polar and weak polar contacts, Arpeggio ('polar distance', 'weak polar distance').
    'polar': 3.5,
    'weak_polar': 3.5,
    # Ionic, RING strict and Arpeggio.
    'ionic': 4.0,
    # pi-pi stacking, RING strict = biotite.structure.find_stacking_interactions defaults.
    'stacking_centroid': 6.5,
    'stacking_plane_tolerance': 30.0,
    'stacking_shift_tolerance': 30.0,
    # pi-cation, RING strict = biotite.structure.find_pi_cation_interactions defaults.
    'pication': 5.0,
    'pication_angle_tolerance': 30.0,
    # pi-hydrogen X-H...ring centroid (Steiner & Koellner 2001, J Mol Biol 305:535; RING 4).
    'pihydrogen_h_centroid': 3.0,
    'pihydrogen_normal_angle': 30.0,
    'pihydrogen_xh_centroid_angle': 120.0,
    # Disulfide, RING.
    'disulfide': 2.5,
    # Metal coordination, Arpeggio ('metal').
    'metal': 2.8,
    # Halogen bond, Arpeggio: theta1 = C-X...A, theta2 = X...A-Y; X...A within the Bondi vdW sum.
    # X...ring centroid for pi acceptors (RING 4 allows pi systems as acceptors).
    'halogen_theta1': 120.0,
    'halogen_theta2_min': 70.0,
    'halogen_theta2_max': 170.0,
    'halogen_pi_centroid': 4.0,
    # Hydrophobic, Arpeggio.
    'hydrophobic': 4.5,
}

# Radius of the neighbour search for all atom-atom families.
MAX_ATOM_PAIR_DISTANCE = max(
    THRESHOLDS['polar'], THRESHOLDS['weak_polar'], THRESHOLDS['ionic'], THRESHOLDS['hydrophobic'],
    THRESHOLDS['metal'], THRESHOLDS['disulfide'], THRESHOLDS['weak_hbond_donor_acceptor'],
    2 * max(BONDI_RADII.values()),
)
