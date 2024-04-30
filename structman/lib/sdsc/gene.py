class Gene:
    __slots__ = [
                    'gene_id', 'proteins', 'database_id', 'gene_name'
                ]

    def __init__(self, gene_id, proteins = None, gene_name = None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        if proteins is not None:
            self.proteins = set(proteins)
        else:
            self.proteins = set([])
        self.database_id = None
