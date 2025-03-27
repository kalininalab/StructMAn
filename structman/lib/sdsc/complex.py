from structman.lib import pdbParser 
from structman.lib.sdsc.sdsc_utils import doomsday_protocol, is_alphafold_model, alphafold_model_id_to_file_path, Slotted_obj

class Complex(Slotted_obj):
    __slots__ = [
        'pdb_id',       'resolution',           'interaction_partners',
        'database_id',  'chains',               'chainlist',
        'stored',       'lig_profile',          'metal_profile',
        'ion_profile',  'chain_chain_profile',  'homomers',
        'atom_count',   'page',                 'IAmap',
        'interfaces'
    ]

    slot_mask = [
        False,  True,   False,
        False,  True,   False,
        False,  False,  False,
        False,  False,  False,
        False,  False,  False,
        True    
    ]

    def __init__(self, pdb_id = None, resolution=None, chains_str=None, lig_profile=None, lig_profile_str=None, metal_profile=None,
                 metal_profile_str=None, ion_profile=None, ion_profile_str=None, chain_chain_profile=None,
                 chain_chain_profile_str=None, stored=False, database_id=None, homomers=None, homomers_str=None, atom_count=None,
                 IAmap = None, interfaces = None):
        self.pdb_id = pdb_id
        self.chains = {}  # {One-letter-chain-id:chaintype} ; possible chaintypes: 'Protein','DNA','RNA,'Peptide'
        self.chainlist = None
        self.resolution = resolution
        self.database_id = database_id
        self.interaction_partners = []
        self.stored = stored
        self.lig_profile = lig_profile
        self.metal_profile = metal_profile
        self.ion_profile = ion_profile
        self.chain_chain_profile = chain_chain_profile
        if homomers is None:
            self.homomers = {}
        else:
            self.homomers = homomers  # Output of pdbParser.getInfo
        self.atom_count = atom_count
        self.page = None
        if IAmap is None:
            self.IAmap = {}
        else:
            self.IAmap = IAmap
        if interfaces is None:
            self.interfaces = {}
        else:
            self.interfaces = interfaces

        if homomers_str is not None:
            self.parseHomomersStr(homomers_str)

        if chains_str is not None:
            self.parseChainStr(chains_str)

        if lig_profile_str is not None:
            self.lig_profile = self.parseProfileStr(lig_profile_str)

        if metal_profile_str is not None:
            self.metal_profile = self.parseProfileStr(metal_profile_str)

        if ion_profile_str is not None:
            self.ion_profile = self.parseProfileStr(ion_profile_str)

        if chain_chain_profile_str is not None:
            self.chain_chain_profile = self.parseProfileStr(chain_chain_profile_str)

    def __serialize__(self) -> list[any]:
        serialized_complex: list[any] = []
        for slot_number, attribute_name in enumerate(self.__slots__):
            if self.slot_mask[slot_number]:
                serialized_complex.append(self.__getattribute__(attribute_name))
            else:
                serialized_complex.append(None)
        return serialized_complex

    def deconstruct(self):
        del self.chain_chain_profile
        del self.ion_profile
        del self.metal_profile
        del self.lig_profile
        for chain in self.interfaces:
            for chain2 in self.interfaces[chain]:
                self.interfaces[chain][chain2].deconstruct()
        del self.interfaces
        del self.chains
        doomsday_protocol(self)

    def parseHomomersStr(self, homomers_str):
        homomers = {}
        for chains in homomers_str.split(','):
            for chain in chains:
                homomers[chain] = []
                for chain2 in chains:
                    homomers[chain].append(chain2)
        self.homomers = homomers

    def get_homomers(self, chain: str) -> str:
        if chain not in self.homomers:
            return chain
        return self.homomers[chain]

    def getHomomersStr(self):
        return ','.join(self.getHomomersSets())

    def getHomomersSets(self):
        used_chains = set()
        parts = []
        for chain in self.homomers:
            if chain in used_chains:
                continue
            part = chain
            used_chains.add(chain)
            for chain2 in self.homomers[chain]:
                if chain2 == chain:
                    continue
                used_chains.add(chain2)
                part += chain2
            parts.append(part)
        return parts

    def getChainStr(self):
        chain_str_parts = []
        for chain in self.chainlist:
            chain_str_parts.append('%s:%s' % (chain, self.chains[chain]))
        return ','.join(chain_str_parts)

    def parseChainStr(self, chain_str):
        chainlist = []
        for chain_str_part in chain_str.split(','):
            if chain_str_part.count(':') < 1:
                continue
            chain, chaintype = chain_str_part.split(':')
            self.chains[chain] = chaintype
            chainlist.append(chain)
        self.chainlist = chainlist

    def parseProfileStr(self, profile_str):
        profile = {}
        if profile_str == '':
            return profile
        parts = profile_str.split(',')
        for part in parts:
            id_part, value_part = part.split(':')
            id_1, id_2 = id_part.split('_')
            deg, score = value_part.split('_')
            deg = int(deg)
            score = float(score)
            profile[(id_1, id_2)] = (deg, score)
        return profile

    def getProfileStr(self, profile):
        if profile is None:
            return ''
        profile_str_parts = []
        for id_1 in profile:
            for id_2 in profile[id_1].get_keys():
                deg, score = profile[id_1].get_item(id_2)
                profile_str_parts.append('%s_%s:%s_%s' % (id_1, id_2, str(deg), str(score)))
        return ','.join(profile_str_parts)

    def getLigProfileStr(self):
        return self.getProfileStr(self.lig_profile)

    def set_IAmap(self, IAmap):
        self.IAmap = IAmap

    def get_IAmap(self):
        return self.IAmap

    def set_interfaces(self, interfaces):
        self.interfaces = interfaces

    def get_interfaces(self):
        return self.interfaces

    def set_metal_profile(self, value):
        self.metal_profile = value

    def get_metal_profile(self):
        return self.metal_profile

    def getMetalProfileStr(self):
        return self.getProfileStr(self.metal_profile)

    def getIonProfileStr(self):
        return self.getProfileStr(self.ion_profile)

    def getChainChainProfileStr(self):
        return self.getProfileStr(self.chain_chain_profile)

    def set_chain_type_map(self, value, chainlist):
        self.chains = value
        self.chainlist = chainlist
        del_list = []
        for chain in self.homomers:
            if chain not in self.chains:
                del_list.append(chain)
            else:
                del_pos = []
                for pos, h_chain in enumerate(self.homomers[chain]):
                    if h_chain not in self.chains:
                        del_pos.append(pos)
                for pos in reversed(del_pos):
                    del self.homomers[chain][pos]
        for chain in del_list:
            del self.homomers[chain]

    def set_interaction_partners(self, value):
        self.interaction_partners = value

    def get_interaction_partners(self):
        return self.interaction_partners

    def countLigands(self, tchain):
        count = 0
        for iap in self.interaction_partners:
            ia_type = iap[0]
            if ia_type != "Ligand":
                continue
            chain = iap[3]
            if tchain != chain:
                continue
            count += 1
        return count

    def get_stored(self):
        return self.stored

    def get_chains(self, only_protein=False):
        if only_protein:
            chains = {}
            for chain in self.chains:
                if self.chains[chain] != 'Protein':
                    continue
                chains[chain] = self.chains[chain]
            return chains
        return self.chains

    def set_database_id(self, value):
        self.database_id = value

    def get_atom_count(self):
        return self.atom_count

    def set_atom_count(self, atom_count):
        self.atom_count = atom_count

    def print_interfaces(self):
        for chain_1 in self.interfaces:
            for chain_2 in self.interfaces[chain_1]:
                print(chain_1, chain_2, self.interfaces[chain_1][chain_2].interactions)
                print(self.interfaces[chain_1][chain_2].residues)

    def getPage(self, config, self_update=False):
        if self_update:
            if config.model_db_active:
                if is_alphafold_model(self.pdb_id):
                    model_path = alphafold_model_id_to_file_path(self.pdb_id, config)
                else:
                    model_path = None
            else:
                model_path = None
            parse_out = pdbParser.getStandardizedPdbFile(self.pdb_id, config.pdb_path, verbosity=config.verbosity, model_path = model_path)

            if parse_out is None:
                config.errorlog.add_error('pdbParser failed: %s %s %s' % (u_ac, pdb_id, chain))
                return

            (template_page, interaction_partners, chain_type_map, oligo, atom_count, chainlist, rare_residues) = parse_out

            print(f'In complex object getPage - self_update: interaction_partners: {interaction_partners}')

            self.page = template_page
            self.interaction_partners = interaction_partners
            self.chains = chain_type_map
            self.chainlist = chainlist

        if self.page is None:
            if config.model_db_active:
                if is_alphafold_model(self.pdb_id):
                    model_path = alphafold_model_id_to_file_path(self.pdb_id, config)
                else:
                    model_path = None
            else:
                model_path = None
            page, atom_count = pdbParser.standardParsePDB(self.pdb_id, config.pdb_path, model_path = model_path)
            self.page = page
        return self.page
