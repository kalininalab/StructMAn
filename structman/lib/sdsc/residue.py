# sdsc: structman datastructures and classes
from structman.lib.sdsc.consts.residues import METAL_ATOMS, ION_ATOMS
from structman.lib.sdsc.sdsc_utils import classify, doomsday_protocol, SparseArray, Slotted_obj
from typing import TypeVar, Generic

T = TypeVar('T')

class Residue_Map(Slotted_obj, Generic[T]):
    __slots__ = ['value_list', 'value_dict']

    def __init__(self):
        self.value_list = SparseArray()
        self.value_dict = {}

    #def __setitem__(self, key: int | str, value: any) -> None:
    #    return self.add_item(key, value)

    def add_item(self, key: int | str, value: any) -> None:
        if type(key) == int:
            self.add_intkey_item(key, value)
        else:
            try:
                int_key = int(key)
                self.add_intkey_item(int_key, value)
            except:
                self.add_strkey_item(key, value)

    def add_intkey_item(self, intkey, value):
        self.value_list.insert(intkey, value)

    def add_strkey_item(self, strkey, value):
        self.value_dict[strkey] = value

    #def __getitem__(self, key: int | str) -> any:
    #    return self.get_item(key)

    def get_item(self, key: int | str) -> any:
        if isinstance(key, int):
            return self.value_list.get(key)
        try:
            intkey = int(key)
            value = self.value_list.get(intkey)
            return value
        except:
            try:
                return self.value_dict[key]
            except:
                return None

    def get_keys(self) -> list[int | str]:
        return self.value_list.get_keys() + list(self.value_dict.keys())

    def contains(self, key: int | str) -> bool:
        return self.get_item(key) is not None

    def __len__(self):
        return len(self.value_list) + len(self.value_dict)

class Residue(Slotted_obj):
    __slots__ = [
        'res_num', 'aa', 'lig_dists',
        'chain_distances', 'RSA', 'relative_main_chain_acc',
        'relative_side_chain_acc', 'SSA', 'homomer_distances',
        'interaction_profile', 'centralities', 'modres',
        'b_factor', 'pLDDT', 'database_id',
        'stored', 'phi', 'psi',
        'intra_ssbond', 'inter_ssbond', 'ssbond_length',
        'intra_link', 'inter_link', 'link_length',
        'cis_conformation', 'cis_follower', 'inter_chain_median_kd',
        'inter_chain_dist_weighted_kd', 'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa',
        'intra_chain_median_kd', 'intra_chain_dist_weighted_kd', 'inter_chain_interactions_median',
        'inter_chain_interactions_dist_weighted', 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
        'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa', 'Class',
        'simpleClass'
        ]

    def __init__(self):
        pass

    def deconstruct(self):
        del self.res_num
        doomsday_protocol(self)

    def get_interaction_partners(self):
        if self.interaction_profile is None:
            return None, None
        return (self.interaction_profile.interacting_chains, self.interaction_profile.interacting_ligands)

    def get_rsa(self, splitted=False):
        if not splitted:
            return self.RSA
        else:
            return (self.RSA, self.relative_main_chain_acc, self.relative_side_chain_acc)

    def get_classification(self, config):
        if self.Class is None and self.simpleClass is None:
            if self.RSA is None:
                sc = None
            else:
                if self.RSA > config.surface_threshold:
                    sc = "Surface"
                else:
                    sc = "Core"
            
            if self.interaction_profile is None:
                config.errorlog.add_warning('Interaction profile is None, classification will fail')
            rin_class, rin_simple_class = classify(self.interaction_profile, sc)
            self.Class = rin_class
            self.simpleClass = rin_simple_class
        return self.Class, self.simpleClass

    def set_classification(self, Class, simpleClass):
        self.Class = Class
        self.simpleClass = simpleClass

    def checked_get_attribute(self, attribute_name: str):
        try:
            value = self.__getattribute__(attribute_name)
        except AttributeError:
            value = None
        return value

    def get_res_info(self):
        res_info = (self.checked_get_attribute('RSA'), self.checked_get_attribute('relative_main_chain_acc'), self.checked_get_attribute('relative_side_chain_acc'), self.checked_get_attribute('SSA'),
                    self.checked_get_attribute('interaction_profile'), self.checked_get_attribute('centralities'),
                    self.checked_get_attribute('phi'), self.checked_get_attribute('psi'), self.checked_get_attribute('intra_ssbond'), self.checked_get_attribute('inter_ssbond'), self.checked_get_attribute('ssbond_length'),
                    self.checked_get_attribute('intra_link'), self.checked_get_attribute('inter_link'), self.checked_get_attribute('link_length'), self.checked_get_attribute('cis_conformation'), self.checked_get_attribute('cis_follower'),
                    self.checked_get_attribute('inter_chain_median_kd'), self.checked_get_attribute('inter_chain_dist_weighted_kd'),
                    self.checked_get_attribute('inter_chain_median_rsa'), self.checked_get_attribute('inter_chain_dist_weighted_rsa'), self.checked_get_attribute('intra_chain_median_kd'),
                    self.checked_get_attribute('intra_chain_dist_weighted_kd'), self.checked_get_attribute('intra_chain_median_rsa'), self.checked_get_attribute('intra_chain_dist_weighted_rsa'),
                    self.checked_get_attribute('inter_chain_interactions_median'), self.checked_get_attribute('inter_chain_interactions_dist_weighted'),
                    self.checked_get_attribute('intra_chain_interactions_median'), self.checked_get_attribute('intra_chain_interactions_dist_weighted'),
                    self.checked_get_attribute('b_factor'), self.checked_get_attribute('pLDDT'), self.checked_get_attribute('modres'),
                    self.checked_get_attribute('lig_dists'), self.checked_get_attribute('chain_distances'), self.checked_get_attribute('homomer_distances'), self.checked_get_attribute('aa'))

        return res_info
    
    def set_res_info(self, res_data):
        (
            self.RSA, self.relative_main_chain_acc, self.relative_side_chain_acc,
            self.SSA, self.profile, self.centralities, 
            self.phi, self.psi, self.intra_ssbond,
            self.inter_ssbond, self.ssbond_length, self.intra_link,
            self.inter_link, self.link_length, self.cis_conformation,
            self.cis_follower, self.inter_chain_median_kd, self.inter_chain_dist_weighted_kd,
            self.inter_chain_median_rsa, self.inter_chain_dist_weighted_rsa, self.intra_chain_median_kd,
            self.intra_chain_dist_weighted_kd, self.intra_chain_median_rsa, self.intra_chain_dist_weighted_rsa,
            self.inter_chain_interactions_median, self.inter_chain_interactions_dist_weighted, self.intra_chain_interactions_median,
            self.intra_chain_interactions_dist_weighted, self.b_factor, self.pLDDT,
            self.modres, self.lig_dists, self.chain_distances,
            self.homomer_distances, self.aa
        ) = res_data
 