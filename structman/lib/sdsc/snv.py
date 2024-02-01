from structman.lib.sdsc.sdsc_utils import doomsday_protocol

class SNV:
    __slots__ = ['new_aa', 'database_id', 'stored', 'tags']

    def __init__(self, new_aa, tags=set(), database_id=None):
        self.new_aa = new_aa
        self.tags = tags.copy()
        self.database_id = database_id
        self.stored = (database_id is not None)

    def fuse(self, other_snv):
        self.tags = self.tags | other_snv.tags
        other_snv.deconstruct()
        del other_snv

    def copy(self):
        copy_snv = SNV(self.new_aa, database_id = self.database_id, tags = self.tags)
        return copy_snv

    def deconstruct(self):
        doomsday_protocol(self)
