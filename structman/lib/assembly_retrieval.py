"""PDB-id inputs: which assembly chains they mean, and their sequences from the mmCIF assemblies.

Accepted input forms (case of the PDB id does not matter, chain names are case sensitive):

    1ABC        the protein chains of the first assembly (assembly 1 as a rule)
    1ABC-2      the protein chains of assembly 2
    1ABC:A      chain A of the first assembly that has it
    1ABC-2:A    chain A of assembly 2

and the extended ids `PDB_0000xxxx` the same way. Every input resolves to assembly chain ids
`<PDBID>-<n>:<chain>` (mmcifParser.make_chain_id), the ids of the search DB and the RIN DB. A
chain is matched by name first and then by the name of the deposited chain it was made from, so
`A` also finds `A-2` when assembly n only holds a transformed copy of A.

Sequences are the observed residues of the kept chain instance, read with the same code as the
search DB (mmcifParser.get_observed_protein_chains), so a PDB input chain is identical to its own
search DB entry. Each comes with {author residue id ('52', '52A'): sequence position}, the residue
ids that PDB style mutations ('1ABC:A  K52AE') name.

Files come from the local mirror (<pdb_path>/data/assemblies/mmCIF/divided); an id the mirror
does not have is looked up among its obsolete entries (data/structures/obsolete/mmCIF) and replaced
by its successor, and what is still missing is downloaded from the RCSB.
"""

from __future__ import annotations

import gzip
import io
import os
import re
import urllib.error
import urllib.request

from structman.lib import mmcifParser

RCSB_ASSEMBLY_URL = 'https://files.rcsb.org/download/{pdb_id}-assembly{nr}.cif.gz'
DOWNLOAD_TIMEOUT = 60
# Online, assemblies are tried 1, 2, ... until one is missing; no entry has nearly this many.
MAX_ONLINE_ASSEMBLIES = 64
EXPAND = '-'

_INPUT_RE = re.compile(r'^(?P<pdb>[0-9][A-Za-z0-9]{3}|[Pp][Dd][Bb]_[0-9A-Za-z]{8})(?:-(?P<nr>[1-9][0-9]*))?(?::(?P<chain>[^:]+))?$')


# ---------------------------------------------------------------------------------------------
# Input ids

def parse_pdb_input(text: str) -> tuple[str, int | None, str | None] | None:
    """(PDB id upper case, assembly number or None, chain or None), None if text is no PDB input."""
    match = _INPUT_RE.match(text.strip())
    if match is None:
        return None
    chain = match.group('chain')
    if chain == EXPAND:
        chain = None
    elif chain is not None and not mmcifParser.is_valid_chain_name(chain):
        return None
    return match.group('pdb').upper(), int(match.group('nr')) if match.group('nr') else None, chain


def input_key(pdb_id: str, assembly_nr: int | None, chain: str | None) -> str:
    """The key a PDB input is filed under until it is resolved: 1ABC:A, 1ABC:-, 1ABC-2:-, 1ABC-2:A."""
    structure = pdb_id if assembly_nr is None else mmcifParser.make_assembly_id(pdb_id, assembly_nr)
    return f'{structure}:{chain if chain is not None else EXPAND}'


# ---------------------------------------------------------------------------------------------
# Files

def _log(config, text: str, level: int = 2) -> None:
    if getattr(config, 'verbosity', 0) >= level:
        logger = getattr(config, 'logger', None)
        if logger is not None:
            logger.info(text)
        else:
            print(text)


def _urlopen(url: str):
    return urllib.request.urlopen(urllib.request.Request(url), timeout = DOWNLOAD_TIMEOUT)


def _wait_for_connection(config, url: str) -> None:
    from structman.lib.sdsc import sdsc_utils
    sdsc_utils.connection_sleep_cycle(getattr(config, 'verbosity', 0), url)


def local_assembly_numbers(pdb_path: str | None, pdb_id: str) -> list[int]:
    """The assembly numbers of the entry present in the local mirror, ascending."""
    if not pdb_path:
        return []
    folder = os.path.dirname(mmcifParser.assembly_file_path(pdb_path, pdb_id, 1))
    try:
        names = os.listdir(folder)
    except OSError:
        return []
    prefix = f'{pdb_id.lower()}-assembly'
    numbers = []
    for name in names:
        if not name.startswith(prefix):
            continue
        parsed = mmcifParser.parse_assembly_filename(name)
        if parsed is not None and parsed[0].upper() == pdb_id.upper():
            numbers.append(parsed[1])
    return sorted(numbers)


def obsolete_successor(pdb_path: str | None, pdb_id: str) -> str | None:
    if not pdb_path:
        return None
    path = mmcifParser.obsolete_file_path(pdb_path, pdb_id)
    if not os.path.isfile(path):
        return None
    try:
        with mmcifParser.open_cif_text(path) as f:
            return mmcifParser.read_obsolete_replacement(f)
    except (OSError, ValueError):
        return None


def download_assembly(pdb_id: str, assembly_nr: int, config) -> str | None:
    """The text of an assembly file from the RCSB; None if it does not exist there."""
    url = RCSB_ASSEMBLY_URL.format(pdb_id = pdb_id.upper(), nr = assembly_nr)
    _wait_for_connection(config, url)
    _log(config, f'Assembly {pdb_id}-{assembly_nr} is not in the local PDB, downloading {url}')
    try:
        with _urlopen(url) as response:
            blob = response.read()
    except urllib.error.HTTPError as e:
        if e.code != 404:
            _log(config, f'Download of {url} failed: {e}', level = 1)
        return None
    except (urllib.error.URLError, OSError) as e:
        _log(config, f'Download of {url} failed: {e}', level = 1)
        return None
    try:
        return gzip.decompress(blob).decode('utf-8', errors = 'replace')
    except OSError:
        return blob.decode('utf-8', errors = 'replace')


class AssemblyReader:
    """Reads (and remembers) the protein chains of assemblies, from the mirror or the RCSB.

    One reader per batch of inputs, so that a file is read once however many inputs name it.
    """

    def __init__(self, config):
        self.config = config
        self.pdb_path = getattr(config, 'pdb_path', None) or None
        self._chains: dict[tuple[str, int], dict | None] = {}
        self._local_numbers: dict[str, list[int]] = {}
        self._successor: dict[str, str] = {}

    def local_numbers(self, pdb_id: str) -> list[int]:
        if pdb_id not in self._local_numbers:
            self._local_numbers[pdb_id] = local_assembly_numbers(self.pdb_path, pdb_id)
        return self._local_numbers[pdb_id]

    def current_id(self, pdb_id: str) -> str:
        """The id itself, or its successor when the mirror has it only as an obsolete entry."""
        if pdb_id in self._successor:
            return self._successor[pdb_id]
        current = pdb_id
        if self.pdb_path and not self.local_numbers(pdb_id):
            successor = obsolete_successor(self.pdb_path, pdb_id)
            if successor is not None and successor != pdb_id:
                _log(self.config, f'{pdb_id} is obsolete, replaced by {successor}', level = 1)
                current = successor
        self._successor[pdb_id] = current
        return current

    def chains(self, pdb_id: str, assembly_nr: int) -> dict[str, tuple[str, dict[str, int]]] | None:
        """{chain: (observed sequence, {residue id: position})} of the protein chains; None if the
        assembly does not exist."""
        key = (pdb_id, assembly_nr)
        if key in self._chains:
            return self._chains[key]
        result = None
        try:
            local = self.local_numbers(pdb_id)
            if assembly_nr in local:
                extract = mmcifParser.read_cif_categories(mmcifParser.assembly_file_path(self.pdb_path, pdb_id, assembly_nr))
            elif local:
                # The mirror has the entry, so it has all of its assemblies.
                extract = None
            else:
                text = download_assembly(pdb_id, assembly_nr, self.config)
                extract = None if text is None else mmcifParser.extract_cif_categories(
                    io.StringIO(text), mmcifParser.SEQUENCE_CATEGORIES, required = mmcifParser.REQUIRED_SEQUENCE_CATEGORIES)
            if extract is not None:
                result, _ = mmcifParser.get_observed_protein_residues(extract)
        except (OSError, ValueError) as e:
            _log(self.config, f'Could not read assembly {pdb_id}-{assembly_nr}: {type(e).__name__}: {e}', level = 1)
            result = None
        self._chains[key] = result
        return result

    def candidate_numbers(self, pdb_id: str):
        """Assembly numbers to try, in order: the local ones, else 1, 2, ... online."""
        local = self.local_numbers(pdb_id)
        if local:
            yield from local
            return
        for nr in range(1, MAX_ONLINE_ASSEMBLIES + 1):
            if self.chains(pdb_id, nr) is None:
                return
            yield nr


# ---------------------------------------------------------------------------------------------
# Resolving inputs

def _match_chain(chains, wanted: str) -> str | None:
    if wanted in chains:
        return wanted
    for chain in chains:
        if chain.split('-', 1)[0] == wanted:
            return chain
    return None


def resolve_inputs(keys, config, reader: AssemblyReader | None = None) -> tuple[dict[str, list[str]], dict[str, str]]:
    """({input key: [assembly chain ids]}, {input key: why it could not be resolved})."""
    reader = reader or AssemblyReader(config)
    resolved: dict[str, list[str]] = {}
    errors: dict[str, str] = {}
    for key in keys:
        parsed = parse_pdb_input(key)
        if parsed is None:
            errors[key] = 'not a PDB id'
            continue
        pdb_id, assembly_nr, chain = parsed
        current = reader.current_id(pdb_id)
        numbers = [assembly_nr] if assembly_nr is not None else reader.candidate_numbers(current)

        found: list[str] = []
        exact_only = True
        checked = []
        for nr in numbers:
            chains = reader.chains(current, nr)
            if chains is None:
                continue
            checked.append(nr)
            assembly_id = mmcifParser.make_assembly_id(current, nr)
            if chain is None:
                found = [mmcifParser.make_chain_id(assembly_id, name) for name in chains]
                break
            if chain in chains:
                found = [mmcifParser.make_chain_id(assembly_id, chain)]
                break
        if chain is not None and not found:
            # No assembly has the chain under its own name: a transformed copy of it then.
            exact_only = False
            for nr in checked:
                match = _match_chain(reader.chains(current, nr), chain)
                if match is not None:
                    found = [mmcifParser.make_chain_id(mmcifParser.make_assembly_id(current, nr), match)]
                    break

        if found:
            resolved[key] = found
            if not exact_only:
                _log(config, f'{key}: chain {chain} resolved to {found[0]}')
        elif not checked:
            what = f'assembly {assembly_nr} of {current}' if assembly_nr is not None else current
            errors[key] = f'{what} was found neither in the local PDB nor at the RCSB'
        elif chain is None:
            errors[key] = f'{current} has no protein chain in assembly {checked[0]}'
        else:
            errors[key] = f'no protein chain {chain} in the assemblies {", ".join(map(str, checked))} of {current}'
    return resolved, errors


# ---------------------------------------------------------------------------------------------
# Sequences

def get_sequences(prot_ids, config, reader: AssemblyReader | None = None) -> list[tuple[str, str, dict[str, int]]]:
    """(assembly chain id, observed sequence, {residue id: position}) for each id.

    An id whose assembly or chain cannot be read comes back with an empty sequence, which the
    sequence scan takes as a protein to remove.
    """
    reader = reader or AssemblyReader(config)
    results = []
    for prot_id in prot_ids:
        try:
            structure_id, chain = mmcifParser.split_chain_id(prot_id)
            pdb_id, assembly_nr = mmcifParser.split_assembly_id(structure_id)
        except ValueError:
            results.append((prot_id, '', {}))
            continue
        chains = reader.chains(pdb_id, assembly_nr)
        if chains is None or chain not in chains:
            _log(config, f'No sequence for {prot_id}', level = 3)
            results.append((prot_id, '', {}))
            continue
        seq, res_pos_map = chains[chain]
        results.append((prot_id, seq, res_pos_map))
    return results
