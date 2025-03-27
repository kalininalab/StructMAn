from importlib import import_module

try:
    from structguy.support_classes import CrossValidationSlice, Feature
except:
#    [e, f, g] = sys.exc_info()
#    g = traceback.format_exc()
#    errortext = '\n'.join([str(e), str(f), str(g)]) + '\n\n'
#    print(errortext)
    pass

#https://stackoverflow.com/questions/2020014/get-fully-qualified-class-name-of-an-object-in-python
def fullname(o):
    klass = o.__class__
    module = klass.__module__
    if module == 'builtins':
        return klass.__qualname__ # avoid outputs like 'builtins.str'
    return f'{module}.{klass.__qualname__}'

def custom_encoder(obj):

    if isinstance(obj, set):
        return {'__set__': True, 'as_list': list(obj)}

    try:
        serialized_object = obj.__serialize__()
        return {f'__{fullname(obj)}__': True, 'as_list': serialized_object}
    except AttributeError:
        pass

    return obj

def standard_decode(serialized_object: list[any], class_init: callable) -> any:
    rebuild_obj = class_init()
    try:
        slot_mask: list[bool] = rebuild_obj.slot_mask
        i = 0
        for slot_number, attribute_name in enumerate(rebuild_obj.__slots__):
            if slot_mask[slot_number]:
                rebuild_obj.__setattr__(attribute_name, serialized_object[i])
                i += 1

    except AttributeError:
        for i, attribute_name in enumerate(rebuild_obj.__slots__):
            rebuild_obj.__setattr__(attribute_name, serialized_object[i])
    return rebuild_obj


def custom_decoder(obj: dict[str, any]) -> any:
    if len(obj) != 2:
        return obj

    for key in obj:
        if key != 'as_list':
            dunder_class_name = key
        else:
            serialized_object = obj[key]

    if not isinstance(dunder_class_name, str):
        return obj

    if '__set__' == dunder_class_name:
        return set(serialized_object)

    if dunder_class_name[:2] == '__':
        full_class_name: str = dunder_class_name[2:-2]
        module_path, class_name = full_class_name.rsplit('.', 1)
        module = import_module(module_path)
        try:
            class_obj = getattr(module, class_name)
        except AttributeError as e:
            print(f'Error in decoding {full_class_name}: {e}')

        return standard_decode(serialized_object, class_obj)

    return obj