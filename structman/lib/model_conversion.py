"""PDB format models (AlphaFold DB, AlphaFold dimers, complex models) as mmCIF.

The model files carry little more than ATOM records. The mmCIF written here adds what the mmCIF side
of StructMAn reads:

- `_entity` / `_entity_poly`: chains with the same polymer sequence share an entity, typed by
  residue names (CCD amino acids -> polypeptide(L), nucleotides -> polydeoxyribonucleotide /
  polyribonucleotide); every other residue name gets a non-polymer entity. The RIN engine decides
  what is a protein residue from these.
- `_pdbx_poly_seq_scheme`: every polymer residue, all observed (models have no gaps), which
  mmcifParser reads the sequences from.
- `_atom_site` with `label_seq_id` = position in the chain (1..n); the author numbering, which
  for the complex models starts somewhere in the protein, stays in auth_seq_id. B-factors
  (pLDDT in AlphaFold models, 0-100) stay in B_iso_or_equiv.
"""

from __future__ import annotations

import io

import numpy as np

import biotite.structure as struc
import biotite.structure.info as info
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

from structman.lib import mmcifParser

_AMINO_ACIDS = frozenset(info.amino_acid_names())
_NUCLEOTIDES = frozenset(info.nucleotide_names())


def _poly_type(res_names: list[str]) -> str | None:
    if res_names and all(name in _AMINO_ACIDS for name in res_names):
        return 'polypeptide(L)'
    if res_names and all(name in _NUCLEOTIDES for name in res_names):
        return 'polydeoxyribonucleotide' if all(name.startswith('D') for name in res_names) else 'polyribonucleotide'
    return None


def _category(columns: dict) -> pdbx.CIFCategory:
    return pdbx.CIFCategory({name: [str(value) for value in values] for name, values in columns.items()})


def read_pdb_atoms(pdb_text: str) -> struc.AtomArray:
    return pdb.PDBFile.read(io.StringIO(pdb_text)).get_structure(model = 1, altloc = 'first', extra_fields = ['b_factor', 'occupancy'])


def atoms_to_mmcif(atoms: struc.AtomArray, entry_id: str) -> tuple[str, dict[str, str]]:
    """(mmCIF text, {chain: one letter sequence of its polymer}) of a model."""
    starts = struc.get_residue_starts(atoms)
    residue_of_atom = np.zeros(atoms.array_length(), dtype = np.int64)
    residue_of_atom[starts[1:]] = 1
    residue_of_atom = np.cumsum(residue_of_atom)

    # Chains in order of appearance; polymer = the chain's residues, if they all are amino acids
    # (or all nucleotides). Anything else in the chain is a non-polymer of its own.
    chain_order = list(dict.fromkeys(atoms.chain_id[starts].tolist()))
    entity_of_residue = np.empty(len(starts), dtype = object)
    seq_id_of_residue = np.full(len(starts), '.', dtype = object)
    entities: dict = {}            # key -> entity id
    entity_rows = []               # (id, type)
    poly_rows = []                 # (entity id, poly type, sequence, strands)
    poly_strands: dict[str, list[str]] = {}
    scheme_rows = []
    chains: dict[str, str] = {}

    for chain in chain_order:
        residues = [index for index, start in enumerate(starts) if atoms.chain_id[start] == chain]
        names = [atoms.res_name[starts[index]] for index in residues]
        polymer = [index for index, name in zip(residues, names) if name in _AMINO_ACIDS or name in _NUCLEOTIDES]
        poly_names = [atoms.res_name[starts[index]] for index in polymer]
        poly_type = _poly_type(poly_names) if len(polymer) > 1 else None
        if poly_type is not None:
            key = ('polymer', poly_type, tuple(poly_names))
            if key not in entities:
                entities[key] = str(len(entities) + 1)
                entity_rows.append((entities[key], 'polymer'))
                one_letter = ''.join(mmcifParser._ONE_LETTER.get(name, 'X') for name in poly_names)
                poly_rows.append([entities[key], poly_type, one_letter])
                poly_strands[entities[key]] = []
            entity_id = entities[key]
            poly_strands[entity_id].append(chain)
            for seq_id, index in enumerate(polymer, start = 1):
                entity_of_residue[index] = entity_id
                seq_id_of_residue[index] = str(seq_id)
                start = starts[index]
                name = atoms.res_name[start]
                ins_code = atoms.ins_code[start] or '.'
                scheme_rows.append((chain, entity_id, seq_id, name, seq_id, int(atoms.res_id[start]), int(atoms.res_id[start]),
                                    name, name, chain, ins_code, 'n'))
            if poly_type == 'polypeptide(L)':
                chains[chain] = ''.join(mmcifParser._ONE_LETTER.get(name, 'X') for name in poly_names)
            polymer_set = set(polymer)
        else:
            polymer_set = set()
        for index in residues:
            if index in polymer_set:
                continue
            key = ('non-polymer', atoms.res_name[starts[index]])
            if key not in entities:
                entities[key] = str(len(entities) + 1)
                entity_rows.append((entities[key], 'non-polymer'))
            entity_of_residue[index] = entities[key]

    atoms = atoms.copy()
    atoms.set_annotation('label_entity_id', entity_of_residue[residue_of_atom].astype(str))
    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, atoms, data_block = entry_id)
    block = cif.block
    atom_site = block['atom_site']
    atom_site['label_seq_id'] = pdbx.CIFColumn(pdbx.CIFData(seq_id_of_residue[residue_of_atom].astype(str).tolist()))

    block['entry'] = _category({'id': [entry_id]})
    block['entity'] = _category({'id': [row[0] for row in entity_rows], 'type': [row[1] for row in entity_rows]})
    if poly_rows:
        block['entity_poly'] = _category({
            'entity_id': [row[0] for row in poly_rows],
            'type': [row[1] for row in poly_rows],
            'pdbx_seq_one_letter_code_can': [row[2] for row in poly_rows],
            'pdbx_strand_id': [','.join(poly_strands[row[0]]) for row in poly_rows],
        })
        names = ('asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num', 'pdb_seq_num', 'auth_seq_num',
                 'pdb_mon_id', 'auth_mon_id', 'pdb_strand_id', 'pdb_ins_code', 'hetero')
        block['pdbx_poly_seq_scheme'] = _category({name: [row[pos] for row in scheme_rows] for pos, name in enumerate(names)})

    # atom_site last, as in the PDB's files: readers that stop at it see all other categories.
    del block['atom_site']
    block['atom_site'] = atom_site

    out = io.StringIO()
    cif.write(out)
    return out.getvalue(), chains


def pdb_to_mmcif(pdb_text: str, entry_id: str) -> tuple[str, dict[str, str]]:
    return atoms_to_mmcif(read_pdb_atoms(pdb_text), entry_id)
