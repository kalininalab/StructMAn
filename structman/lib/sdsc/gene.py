class Gene:
    __slots__ = [
                    'gene_id', 'proteins', 'database_id'
                ]

    def __init__(self, gene_id, proteins = None):
        self.gene_id = gene_id
        if proteins is not None:
            self.proteins = set(proteins)
        else:
            self.proteins = set([])
        self.database_id = None
