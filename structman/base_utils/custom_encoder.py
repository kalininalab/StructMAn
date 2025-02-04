import structman.lib.sdsc as sdsc
import copy

try:
    from structguy.support_classes import CrossValidationSlice, Feature
except:
#    [e, f, g] = sys.exc_info()
#    g = traceback.format_exc()
#    errortext = '\n'.join([str(e), str(f), str(g)]) + '\n\n'
#    print(errortext)
    pass

def custom_encoder_complete(obj):
    complete = True

    if isinstance(obj, set):
        return {'__set__': True, 'as_list': list(obj)}

    if 'Residue' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_residue = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'interaction_profile' and not complete:
                serialized_residue.append(None)
            elif attribute_name == 'centralities' and not complete:
                serialized_residue.append(None)
            else:
                serialized_residue.append(obj.__getattribute__(attribute_name))
        return {'__residue__': True, 'as_list': serialized_residue}

    if 'StructureAnnotation' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_annotation = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'alignment' and not complete:
                serialized_annotation.append(None)
            else:
                serialized_annotation.append(obj.__getattribute__(attribute_name))
        return {'__structureannotation__': True, 'as_list': serialized_annotation}

    if 'Structure' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_structure = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'sequence' and not complete:
                serialized_structure.append(None)
            else:
                serialized_structure.append(obj.__getattribute__(attribute_name))
        return {'__structure__': True, 'as_list': serialized_structure}

    if 'Complex' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_complex = []
        for attribute_name in obj.__slots__:
            #Complex objects get slimmed down, when packed !!!!!!!!!
            if (attribute_name == 'chains' or attribute_name == 'resolution' or attribute_name == 'interfaces') and not complete:
                serialized_complex.append(obj.__getattribute__(attribute_name))
            else:
                serialized_complex.append(None)
        return {'__complex__': True, 'as_list': serialized_complex}

    if 'Protein' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_protein = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'sequence' and not complete:
                serialized_protein.append(None)
            else:
                serialized_protein.append(obj.__getattribute__(attribute_name))
        return {'__protein__': True, 'as_list': serialized_protein}

    if 'Position_Position_Interaction' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_pos_pos_i = []
        for attribute_name in obj.__slots__:
            serialized_pos_pos_i.append(obj.__getattribute__(attribute_name))
        return {'__pos_pos_i__': True, 'as_list': serialized_pos_pos_i}

    if 'Position' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_position = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'mappings' and not complete:
                serialized_position.append(None)
            #elif attribute_name == 'mut_aas':
            #    serialized_position.append(None)
            else:
                serialized_position.append(obj.__getattribute__(attribute_name))
        return {'__position__': True, 'as_list': serialized_position}

    if 'Aggregated_interface' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_agg_interface = []
        for attribute_name in obj.__slots__:
            serialized_agg_interface.append(obj.__getattribute__(attribute_name))
        return {'__aggregated_interface__': True, 'as_list': serialized_agg_interface}

    if 'Interface' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_interface = []
        for attribute_name in obj.__slots__:
            serialized_interface.append(obj.__getattribute__(attribute_name))
        return {'__interface__': True, 'as_list': serialized_interface}


    if 'Insertion' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_insertion = []
        for attribute_name in obj.__slots__:
            serialized_insertion.append(obj.__getattribute__(attribute_name))
        return {'__insertion__': True, 'as_list': serialized_insertion}


    if 'Deletion' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_deletion = []
        for attribute_name in obj.__slots__:
            serialized_deletion.append(obj.__getattribute__(attribute_name))
        return {'__deletion__': True, 'as_list': serialized_deletion}


    if 'Substitution' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_substitution = []
        for attribute_name in obj.__slots__:
            serialized_substitution.append(obj.__getattribute__(attribute_name))
        return {'__substitution__': True, 'as_list': serialized_substitution}

    if 'SNV' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_snv = []
        for attribute_name in obj.__slots__:
            serialized_snv.append(obj.__getattribute__(attribute_name))
        return {'__snv__': True, 'as_list': serialized_snv}

    if 'Microminer_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Microminer_features__': True, 'as_list': serialized_object}
    
    if 'Structural_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Structural_features__': True, 'as_list': serialized_object}
    
    if 'Integrated_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Integrated_features__': True, 'as_list': serialized_object}

    if 'RIN_based_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = obj.get_raw_list()
        return {'__RIN_based_features__': True, 'as_list': serialized_object}

    if 'CrossValidationSlice' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            #serialized_object.append(copy.deepcopy(obj.__getattribute__(attribute_name)))
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__CrossValidationSlice__': True, 'as_list': serialized_object}

    if 'Feature' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            #serialized_object.append(copy.deepcopy(obj.__getattribute__(attribute_name)))
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Feature__': True, 'as_list': serialized_object}

    return obj

def custom_encoder(obj):
    complete = False

    if isinstance(obj, set):
        return {'__set__': True, 'as_list': list(obj)}

    if 'Residue_Map' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Residue_Map__': True, 'as_list': serialized_object}

    if 'Residue' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_residue = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'interaction_profile' and not complete:
                serialized_residue.append(None)
            elif attribute_name == 'centralities' and not complete:
                serialized_residue.append(None)
            else:
                serialized_residue.append(obj.__getattribute__(attribute_name))
        return {'__residue__': True, 'as_list': serialized_residue}

    if 'StructureAnnotation' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_annotation = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'alignment' and not complete:
                serialized_annotation.append(None)
            else:
                serialized_annotation.append(obj.__getattribute__(attribute_name))
        return {'__structureannotation__': True, 'as_list': serialized_annotation}

    if 'Structure' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_structure = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'sequence' and not complete:
                serialized_structure.append(None)
            else:
                serialized_structure.append(obj.__getattribute__(attribute_name))
        return {'__structure__': True, 'as_list': serialized_structure}

    if 'Complex' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_complex = []
        for attribute_name in obj.__slots__:
            #Complex objects get slimmed down, when packed !!!!!!!!!
            if (attribute_name == 'chains' or attribute_name == 'resolution' or attribute_name == 'interfaces') and not complete:
                serialized_complex.append(obj.__getattribute__(attribute_name))
            else:
                serialized_complex.append(None)
        return {'__complex__': True, 'as_list': serialized_complex}

    if 'Protein' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_protein = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'sequence' and not complete:
                serialized_protein.append(None)
            else:
                serialized_protein.append(obj.__getattribute__(attribute_name))
        return {'__protein__': True, 'as_list': serialized_protein}

    if 'Position_Position_Interaction' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_pos_pos_i = []
        for attribute_name in obj.__slots__:
            serialized_pos_pos_i.append(obj.__getattribute__(attribute_name))
        return {'__pos_pos_i__': True, 'as_list': serialized_pos_pos_i}

    if 'Position' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_position = []
        for attribute_name in obj.__slots__:
            if attribute_name == 'mappings' and not complete:
                serialized_position.append(None)
            #elif attribute_name == 'mut_aas':
            #    serialized_position.append(None)
            else:
                serialized_position.append(obj.__getattribute__(attribute_name))
        return {'__position__': True, 'as_list': serialized_position}

    if 'Aggregated_interface' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_agg_interface = []
        for attribute_name in obj.__slots__:
            serialized_agg_interface.append(obj.__getattribute__(attribute_name))
        return {'__aggregated_interface__': True, 'as_list': serialized_agg_interface}

    if 'Interface' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_interface = []
        for attribute_name in obj.__slots__:
            serialized_interface.append(obj.__getattribute__(attribute_name))
        return {'__interface__': True, 'as_list': serialized_interface}


    if 'Insertion' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_insertion = []
        for attribute_name in obj.__slots__:
            serialized_insertion.append(obj.__getattribute__(attribute_name))
        return {'__insertion__': True, 'as_list': serialized_insertion}


    if 'Deletion' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_deletion = []
        for attribute_name in obj.__slots__:
            serialized_deletion.append(obj.__getattribute__(attribute_name))
        return {'__deletion__': True, 'as_list': serialized_deletion}


    if 'Substitution' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_substitution = []
        for attribute_name in obj.__slots__:
            serialized_substitution.append(obj.__getattribute__(attribute_name))
        return {'__substitution__': True, 'as_list': serialized_substitution}

    if 'SNV' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_snv = []
        for attribute_name in obj.__slots__:
            serialized_snv.append(obj.__getattribute__(attribute_name))
        return {'__snv__': True, 'as_list': serialized_snv}

    if 'Microminer_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Microminer_features__': True, 'as_list': serialized_object}
    
    if 'Structural_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Structural_features__': True, 'as_list': serialized_object}
    
    if 'Integrated_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Integrated_features__': True, 'as_list': serialized_object}

    if 'RIN_based_features' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = obj.get_raw_list()
        return {'__RIN_based_features__': True, 'as_list': serialized_object}

    if 'CrossValidationSlice' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            #serialized_object.append(copy.deepcopy(obj.__getattribute__(attribute_name)))
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__CrossValidationSlice__': True, 'as_list': serialized_object}

    if 'Feature' in str(type(obj)): #isinstance just won't work, don't know why
        serialized_object = []
        for attribute_name in obj.__slots__:
            #serialized_object.append(copy.deepcopy(obj.__getattribute__(attribute_name)))
            serialized_object.append(obj.__getattribute__(attribute_name))
        return {'__Feature__': True, 'as_list': serialized_object}




    return obj

def custom_decoder(obj):

    if '__set__' in obj:
        return set(obj['as_list'])


    if '__Residue_Map__' in obj:
        serialized_object = obj['as_list']
        rebuild_obj = sdsc.residue.Residue_Map()
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])
        return rebuild_obj    

    if '__residue__' in obj:
        serialized_residue = obj['as_list']
        res = sdsc.residue.Residue(0)
        for i, attribute_name in enumerate(res.__slots__):
            res.__setattr__(attribute_name, serialized_residue[i])
        return res

    if '__structureannotation__' in obj:
        serialized_annotation = obj['as_list']
        u_ac = serialized_annotation[0]
        pdb_id = serialized_annotation[1]
        chain = serialized_annotation[2]
        anno = sdsc.structure.StructureAnnotation(u_ac, pdb_id, chain)
        for i, attribute_name in enumerate(anno.__slots__[3:]):
            anno.__setattr__(attribute_name, serialized_annotation[i+3])
        return anno

    if '__structure__' in obj:
        serialized_structure = obj['as_list']
        pdb_id = serialized_structure[0]
        chain = serialized_structure[1]
        struct = sdsc.structure.Structure(pdb_id, chain)
        for i, attribute_name in enumerate(struct.__slots__[2:]):
            struct.__setattr__(attribute_name, serialized_structure[i+2])
        return struct

    if '__complex__' in obj:
        serialized_complex = obj['as_list']
        pdb_id = serialized_complex[0]
        compl = sdsc.complex.Complex(pdb_id)
        for i, attribute_name in enumerate(compl.__slots__[1:]):
            compl.__setattr__(attribute_name, serialized_complex[i+1])
        return compl

    if '__protein__' in obj:
        serialized_protein = obj['as_list']
        prot = sdsc.protein.Protein(None)
        for i, attribute_name in enumerate(prot.__slots__):
            prot.__setattr__(attribute_name, serialized_protein[i])
        return prot

    if '__pos_pos_i__' in obj:
        serialized_pos_pos_i = obj['as_list']
        pos_pos_i = sdsc.interface.Position_Position_Interaction(None, None, None, None)
        for i, attribute_name in enumerate(pos_pos_i.__slots__):
            pos_pos_i.__setattr__(attribute_name, serialized_pos_pos_i[i])
        return pos_pos_i

    if '__position__' in obj:
        serialized_position = obj['as_list']
        posi = sdsc.position.Position()
        for i, attribute_name in enumerate(posi.__slots__):
            try:
                posi.__setattr__(attribute_name, serialized_position[i])
            except:
                print('Warning: undefined slot in position object in custom_decoder: ', attribute_name, i, serialized_position)
                posi.__setattr__(attribute_name, None)
        posi.mappings = sdsc.mappings.Mappings()
        return posi

    if '__aggregated_interface__' in obj:
        serialized_agg_interface = obj['as_list']
        agg_inter = sdsc.interface.Aggregated_interface(None)
        for i, attribute_name in enumerate(agg_inter.__slots__):
            agg_inter.__setattr__(attribute_name, serialized_agg_interface[i])
        return agg_inter

    if '__interface__' in obj:
        serialized_interface = obj['as_list']
        inter = sdsc.interface.Interface(None, None)
        for i, attribute_name in enumerate(inter.__slots__):
            inter.__setattr__(attribute_name, serialized_interface[i])
        return inter


    if '__insertion__' in obj:
        serialized_insertion = obj['as_list']
        left_flank = serialized_insertion[0]
        right_flank = serialized_insertion[1]
        inserted_sequence = serialized_insertion[-1]
        insertion = sdsc.indel.Insertion(left_flank = left_flank, right_flank = right_flank, inserted_sequence = inserted_sequence)
        for i, attribute_name in enumerate(insertion.__slots__[2:-1]):
            insertion.__setattr__(attribute_name, serialized_insertion[i+2])
        return insertion


    if '__deletion__' in obj:
        serialized_deletion = obj['as_list']
        left_flank = serialized_deletion[0]
        right_flank = serialized_deletion[1]
        deletion = sdsc.indel.Deletion(left_flank = left_flank, right_flank = right_flank)
        for i, attribute_name in enumerate(deletion.__slots__[2:]):
            deletion.__setattr__(attribute_name, serialized_deletion[i+2])
        return deletion


    if '__substitution__' in obj:
        serialized_substitution = obj['as_list']
        left_flank = serialized_substitution[0]
        right_flank = serialized_substitution[1]
        inserted_sequence = serialized_substitution[-1]
        substitution = sdsc.indel.Substitution(left_flank = left_flank, right_flank = right_flank, inserted_sequence = inserted_sequence)
        for i, attribute_name in enumerate(substitution.__slots__[2:-1]):
            substitution.__setattr__(attribute_name, serialized_substitution[i+2])
        return substitution
    
    if '__snv__' in obj:
        serialized_snv = obj['as_list']
        snv_obj = sdsc.snv.SNV(None)
        for i, attribute_name in enumerate(snv_obj.__slots__):
            snv_obj.__setattr__(attribute_name, serialized_snv[i])
        return snv_obj
    
    if '__Microminer_features__' in obj:
        serialized_object = obj['as_list']
        rebuild_obj = sdsc.mappings.Microminer_features()
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])
        return rebuild_obj
    
    if '__Structural_features__' in obj:
        serialized_object = obj['as_list']
        rebuild_obj = sdsc.mappings.Structural_features()
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])
        return rebuild_obj

    if '__Integrated_features__' in obj:
        serialized_object = obj['as_list']
        rebuild_obj = sdsc.mappings.Integrated_features()
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])
        return rebuild_obj

    if '__RIN_based_features__' in obj:
        raw_rin_based_features = obj['as_list']
        rebuild_obj = sdsc.mappings.RIN_based_features()
        rebuild_obj.set_values(raw_rin_based_features)
        return rebuild_obj

    if '__CrossValidationSlice__' in obj:
        serialized_object = obj['as_list']
        rebuild_obj = CrossValidationSlice(raw_init =  True)
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])
        rebuild_obj.feature_names = list(rebuild_obj.feature_names)

        return rebuild_obj

    if '__Feature__' in obj:
        serialized_object = obj['as_list']
        rebuild_obj = Feature()
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])

        return rebuild_obj


    return obj