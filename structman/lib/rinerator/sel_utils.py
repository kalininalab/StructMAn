# Copyright 2010 Max-Planck-Institut Informatik
#
#    This file is part of RINerator
#
#    RINerator is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RINerator is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RINerator.  If not, see <http://www.gnu.org/licenses/>.


def chains_to_sel_list(comp_name, chains):
    # selection_lst = [[comp_name, comp_type, res_ranges],...]
    sel_comp = [comp_name, 'protein', []]
    for chain in chains:
        if chain == "":
            chain = " "
        sel_chain = [0, chain, [' ', '_', ' '], [' ', '_', ' ']]
        sel_comp[2].append(sel_chain)
    return [sel_comp]


def segs_to_sel_list(comp_name, segments):
    # selection_lst = [[comp_name, comp_type, res_ranges],...]
    sel_comp = [comp_name, 'protein', []]
    for segment in segments:
        [str_model_nr, str_chain, str_start_res, str_end_res] = segment
        model_nr = int(str_model_nr)
        (start_res_hflag, start_res, start_res_icode) = parse_str_res(str_start_res)
        (end_res_hflag, end_res, end_res_icode) = parse_str_res(str_end_res)
        sel_segment = [model_nr, str_chain, [start_res_hflag, start_res, start_res_icode], [end_res_hflag, end_res, end_res_icode]]
        sel_comp[2].append(sel_segment)
    return [sel_comp]


def ligands_to_sel_list(ligands):
    ligand_list = []
    for lig_info in ligands:
        if len(lig_info) != 3:
            continue
        # lig_info = [lig_name, chain, lig_index]
        chain = lig_info[1]
        if chain == "":
            chain = " "
        try:
            icode = lig_info[2].lstrip('0123456789')
            if len(icode) == 0:
                icode = ' '
                res = int(lig_info[2])
            else:
                res = int(lig_info[2][:-len(icode)])
            sel_component = [lig_info[0], 'ligand', [[0, chain, ['H_' + lig_info[0], res, icode], [None, None, None]]]]
        except ValueError:
            sel_component = [lig_info[0], 'ligand', [[0, chain, ['H_' + lig_info[0], '_', ' '], [None, None, None]]]]
        ligand_list.append(sel_component)
    return ligand_list


def parse_str_res(str_res):
    if str_res == '_':
        return (' ', str_res, ' ')
    elif str_res is None or str_res == 'None':
        return (None, None, None)
    else:
        try:
            res = int(str_res)
            return (' ', res, ' ')
        except ValueError:
            res_icode = str_res[-1]
            res = int(str_res[0:-1])
            return (' ', res, res_icode)
