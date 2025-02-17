class Gene:
    __slots__ = [
                    'gene_id', 'proteins', 'database_id', 'gene_name'
                ]

    def __init__(self, gene_id, proteins = None, gene_name = None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        if proteins is not None:
            self.proteins = proteins
        else:
            self.proteins = {}
        self.database_id = None

    def print_content(self):
        print(f'\nPrint out of Gene: {self.gene_id} {self.gene_name}')
        for prot_id in self.proteins:
            print(f'{prot_id} -> {self.proteins[prot_id]}')