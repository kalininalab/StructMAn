from structman.lib.output import out_generator

def init_aggregated_interface_table(interface_file, obj_only = False):
    interface_output = out_generator.OutputGenerator()
    headers = [
        'Input Protein ID', 'Primary Protein ID', 'Uniprot-Ac', 'WT Amino Acid', 'Position',
        'Interface Number', 'Interface Structure Recommendation', 'Position Structure Recommendation'
    ]

    interface_output.add_headers(headers)
    if obj_only:
        return interface_output

    interface_f = open(interface_file, 'a')
    interface_f.write(interface_output.get_header())
    return interface_output, interface_f

def init_protein_protein_table(ppi_file, obj_only = False):
    ppi_output = out_generator.OutputGenerator()

    headers = [
        'Gene A', 'Input Protein ID A', 'Primary Protein ID A', 'Uniprot-Ac A', 'Interface Number A',
        'Gene B', 'Input Protein ID B', 'Primary Protein ID B', 'Uniprot-Ac B', 'Interface Number B',
        'Structure Recommendation', 'Interaction Score'
    ]

    ppi_output.add_headers(headers)
    if obj_only:
        return ppi_output

    ppi_f = open(ppi_file, 'a')
    ppi_f.write(ppi_output.get_header())
    return ppi_output, ppi_f

def init_position_position_table(pos_pos_file, obj_only = False):
    pos_pos_output = out_generator.OutputGenerator()

    headers = [
        'Input Protein ID A', 'Primary Protein ID A', 'Uniprot-Ac A', 'Interface Number A', 'WT Amino Acid A', 'Position A',
        'Input Protein ID B', 'Primary Protein ID B', 'Uniprot-Ac B', 'Interface Number B', 'WT Amino Acid B', 'Position B',
        'Structure Recommendation', 'Interaction Score'
    ]

    pos_pos_output.add_headers(headers)
    if obj_only:
        return pos_pos_output

    pos_pos_f = open(pos_pos_file, 'a')
    pos_pos_f.write(pos_pos_output.get_header())
    return pos_pos_output, pos_pos_f
