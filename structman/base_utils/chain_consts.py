"""Canonical ordering of chain identifiers.

Chain ids are single characters drawn from this 62 symbol alphabet: the live parse
path reads PDB column 22 (`line[21:22]`), and the multi-model renaming in
pdbParser/rin_utils/createCustomDb hands out replacement ids from this same order and
already refuses to go past its end. Fixing the order in one place lets chain-keyed
dicts be replaced by arrays indexed by chain number, and lets sets of chains be
carried as a bitmask, since 62 fits in a uint64.

Anything outside the alphabet - a blank chain id in particular - has no chain number
and is meant to be skipped by the parsers rather than mapped onto something.
"""

CHAIN_ID_ORDER: str = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz'

MAX_NUMBER_OF_CHAINS: int = len(CHAIN_ID_ORDER)

CHAIN_ID_NUMBERS: dict[str, int] = {chain_id: number for number, chain_id in enumerate(CHAIN_ID_ORDER)}


def chain_id_number(chain_id: str) -> int:
    """Chain number for a chain id, or -1 if it is not a usable chain identifier."""
    return CHAIN_ID_NUMBERS.get(chain_id, -1)


def is_usable_chain_id(chain_id: str) -> bool:
    return chain_id in CHAIN_ID_NUMBERS
