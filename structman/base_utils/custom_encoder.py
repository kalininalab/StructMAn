from importlib import import_module
import ast
import numpy as np
from numpy.lib import format as npy_format
from io import BytesIO
from ray import cloudpickle, ObjectRef

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

_NPY_MAGIC = b'\x93NUMPY'
_NPY_V1_HEADER = _NPY_MAGIC + b'\x01\x00'


def ndarray_to_npy_bytes(array: np.ndarray) -> bytes:
    """The .npy payload np.save would write, without going through np.save.

    np.save costs ~8us an array almost independently of its size - BytesIO plus the
    generic write path - where building the header and calling tobytes is ~1.4us. The
    bytes are identical to np.save's, so blobs written here stay readable by np.load and
    by any older reader.

    Anything the fast path cannot express byte-for-byte (object dtype, which np.save
    pickles, or a non C-contiguous array, which it stores in Fortran order) falls back to
    np.save.
    """
    if array.dtype.hasobject or not array.flags.c_contiguous:
        memfile = BytesIO()
        np.save(memfile, array, fix_imports=False)
        return memfile.getvalue()

    header = "{'descr': %s, 'fortran_order': False, 'shape': %s, }" % (
        repr(npy_format.dtype_to_descr(array.dtype)), repr(array.shape))
    # magic(6) + version(2) + header length(2) + header + '\n', padded to a multiple of 64
    padding = (64 - (10 + len(header) + 1) % 64) % 64
    header = f'{header}{" " * padding}\n'

    return (_NPY_V1_HEADER + len(header).to_bytes(2, 'little')
            + header.encode('latin1') + array.tobytes())


def parse_npy_header(header: str) -> dict:
    """Read numpy's .npy header dict without ast.literal_eval.

    literal_eval on this one small fixed string costs ~14us - more than everything else
    in the decode put together. What numpy writes is rigid, so it is read directly, with
    literal_eval kept for anything unexpected (a structured dtype writes a list for
    'descr', not a string).
    """
    descr_at = header.find("'descr': '")
    order_at = header.find("'fortran_order': ")
    shape_at = header.find("'shape': (")

    if descr_at < 0 or order_at < 0 or shape_at < 0:
        return ast.literal_eval(header)

    descr = header[descr_at + 10:header.index("'", descr_at + 10)]
    fortran_order = header[order_at + 17] == 'T'

    shape_text = header[shape_at + 10:header.index(')', shape_at)].strip().rstrip(',')
    shape = tuple(int(extent) for extent in shape_text.split(',')) if shape_text else ()

    return {'descr': descr, 'fortran_order': fortran_order, 'shape': shape}


def npy_bytes_to_ndarray(payload: bytes) -> np.ndarray:
    """Inverse of ndarray_to_npy_bytes, avoiding np.load.

    np.load is ~38us an array regardless of size - header parsing, version negotiation
    and the allow_pickle machinery - against ~0.9us for reading the header here and
    calling frombuffer. The copy is not optional: frombuffer aliases the payload and
    hands back a read-only array, where np.load returns a writable one.

    Falls back to np.load for anything not a plain C-ordered non-object array, so an
    unrecognised payload is slow rather than wrong.
    """
    try:
        if payload[:8] != _NPY_V1_HEADER:
            raise ValueError('not a version 1.0 .npy payload')

        header_length = int.from_bytes(payload[8:10], 'little')
        header_start = 10
        header = parse_npy_header(payload[header_start:header_start + header_length].decode('latin1'))

        if header['fortran_order']:
            raise ValueError('fortran ordered')

        dtype = np.dtype(header['descr'])
        if dtype.hasobject:
            raise ValueError('object dtype')

        data = payload[header_start + header_length:]
        # .copy(): frombuffer would alias the payload and be read-only, np.load is not.
        return np.frombuffer(data, dtype=dtype).reshape(header['shape']).copy()

    except Exception:
        memfile = BytesIO()
        memfile.write(payload)
        memfile.seek(0)
        return np.load(memfile, allow_pickle=True, fix_imports=False)


def custom_encoder(obj):

    if isinstance(obj, set):
        return {'__set__': True, 'as_list': list(obj)}

    if isinstance(obj, np.ndarray):
        return {'__ndarray__': True, 'as_list': ndarray_to_npy_bytes(obj)}
    
    # Any numpy scalar, not just int64. msgpack serialises np.float64 by itself
    # because it subclasses python float, but np.float32 does not subclass
    # anything msgpack knows, so it used to reach here and fail. Residue data is
    # held as float32 now, so scalars read out of it end up in packed objects.
    if isinstance(obj, np.generic):
        return obj.item()

    """
    if isinstance(obj, ObjectRef):
        memfile = BytesIO()
        cloudpickle.dump(obj, memfile)
        del obj
        return {'__ObjectRef__': True, 'as_list': memfile.getvalue()}
    """
    try:
        serialized_object = obj.__serialize__()
        return {f'__{fullname(obj)}__': True, 'as_list': serialized_object}
    except AttributeError:
        pass

    return obj

# The attribute names a class carries in its serialised form, in order. Walking
# __slots__ against slot_mask for every decoded object - and taking an
# AttributeError for the classes that have no mask - was a quarter of the unpack
# cost, and one packed position row decodes dozens of objects. The layout is a
# property of the class, so it is derived once; deactivate_slot_mask and
# reactivate_slot_mask clear this, as they are what can change it.
_DECODE_SLOTS_CACHE: dict[type, tuple[str, ...]] = {}


def clear_decode_caches() -> None:
    _DECODE_SLOTS_CACHE.clear()


def _decoded_attribute_names(class_init: callable) -> tuple[str, ...]:
    names = _DECODE_SLOTS_CACHE.get(class_init)
    if names is None:
        slots = class_init.__slots__
        slot_mask = getattr(class_init, 'slot_mask', None)
        if slot_mask is None:
            names = tuple(slots)
        else:
            names = tuple(name for name, keep in zip(slots, slot_mask) if keep)
        _DECODE_SLOTS_CACHE[class_init] = names
    return names


def standard_decode(serialized_object: list[any], class_init: callable) -> any:
    rebuild_obj = class_init()
    # Indexed rather than zipped, so a truncated payload still raises IndexError.
    for i, attribute_name in enumerate(_decoded_attribute_names(class_init)):
        setattr(rebuild_obj, attribute_name, serialized_object[i])
    return rebuild_obj


# Class objects keyed by the dunder name the encoder writes.
_DECODE_CLASS_CACHE: dict[str, type] = {}


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

    if '__ndarray__' == dunder_class_name:
        return npy_bytes_to_ndarray(serialized_object)
    
    """
    if '__ObjectRef__' == dunder_class_name:
        memfile = BytesIO()
        memfile.write(serialized_object)
        memfile.seek(0)
        new_obj_ref = cloudpickle.load(memfile)
        return new_obj_ref
    """
        
    if dunder_class_name[:2] == '__':
        # Resolved once per class rather than once per decoded object. Every packed
        # position row carries dozens of them - the interaction profile alone is a
        # fifteen object tree - and each one used to rsplit the name and walk
        # import_module, which is a full importlib lookup even when the module is
        # already in sys.modules.
        class_obj = _DECODE_CLASS_CACHE.get(dunder_class_name)
        if class_obj is None:
            full_class_name: str = dunder_class_name[2:-2]
            module_path, class_name = full_class_name.rsplit('.', 1)
            module = import_module(module_path)
            try:
                class_obj = getattr(module, class_name)
            except AttributeError as e:
                print(f'Error in decoding {full_class_name}: {e}')
                return obj
            _DECODE_CLASS_CACHE[dunder_class_name] = class_obj

        return standard_decode(serialized_object, class_obj)

    return obj