"""RIN files: one per assembly, a zstandard-compressed msgpack map.

    format          FORMAT_VERSION
    engine_version  constants.ENGINE_VERSION of the engine that wrote it
    structure_id    '1A04-2'
    source          {'name', 'size', 'mtime_ns'} of the assembly file it was computed from
    params          constants.THRESHOLDS
    families, classes, kinds   the code tables
    mode, stats     how it was computed (region restriction, centrality mode, timings, counts)
    nodes           columns: chain, res_id, ins_code, res_name, label_asym_id, label_seq_id,
                    entity_id, kind, kept, sas_area
    edges           columns: a, b, family, class_a, class_b, count, score, min_distance
    metrics         {'names': [...], 'values': n_nodes x k float32, NaN = not computed}

Numeric columns are stored as {'dtype', 'shape', 'data'} (raw bytes), string columns as lists.
"""

from __future__ import annotations

import os

import msgpack
import numpy as np
import zstandard

FORMAT_VERSION = 1
SUFFIX = '.rin.zst'
_COMPRESSION_LEVEL = 10


def _pack_array(array: np.ndarray) -> dict:
    array = np.ascontiguousarray(array)
    return {'dtype': array.dtype.str, 'shape': list(array.shape), 'data': array.tobytes()}


def _unpack_array(packed: dict) -> np.ndarray:
    return np.frombuffer(packed['data'], dtype = np.dtype(packed['dtype'])).reshape(packed['shape']).copy()


def _encode_columns(columns: dict) -> dict:
    encoded = {}
    for name, column in columns.items():
        column = np.asarray(column)
        if column.dtype.kind in 'US' or column.dtype == object:
            encoded[name] = {'strings': [str(value) for value in column.tolist()]}
        else:
            encoded[name] = _pack_array(column)
    return encoded


def _decode_columns(encoded: dict) -> dict:
    columns = {}
    for name, column in encoded.items():
        if 'strings' in column:
            columns[name] = np.array(column['strings'], dtype = str) if column['strings'] else np.zeros(0, dtype = str)
        else:
            columns[name] = _unpack_array(column)
    return columns


def rin_file_path(db_path: str, structure_id: str) -> str:
    """<db>/<xx>/<structure_id>.rin.zst, xx = the divided folder of the PDB id (as in the mirror)."""
    pdb_id = structure_id.rsplit('-', 1)[0].lower()
    return os.path.join(db_path, pdb_id[-3:-1], f'{structure_id}{SUFFIX}')


def encode(document: dict) -> bytes:
    payload = dict(document)
    payload['nodes'] = _encode_columns(document['nodes'])
    payload['edges'] = _encode_columns(document['edges'])
    metrics = document.get('metrics') or {}
    payload['metrics'] = {
        'names': list(metrics.get('names', [])),
        'values': _pack_array(np.asarray(metrics.get('values', np.zeros((0, 0), dtype = np.float32)), dtype = np.float32)),
    }
    return zstandard.ZstdCompressor(level = _COMPRESSION_LEVEL).compress(msgpack.packb(payload, use_bin_type = True))


def decode(blob: bytes) -> dict:
    document = msgpack.unpackb(zstandard.ZstdDecompressor().decompress(blob), raw = False, strict_map_key = False)
    document['nodes'] = _decode_columns(document['nodes'])
    document['edges'] = _decode_columns(document['edges'])
    document['metrics'] = {
        'names': document['metrics']['names'],
        'values': _unpack_array(document['metrics']['values']),
    }
    return document


def write_rin(path: str, document: dict) -> int:
    """Writes atomically; returns the size in bytes."""
    os.makedirs(os.path.dirname(path), exist_ok = True)
    blob = encode(document)
    tmp_path = f'{path}.tmp.{os.getpid()}'
    try:
        with open(tmp_path, 'wb') as f:
            f.write(blob)
        os.replace(tmp_path, path)
    except BaseException:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return len(blob)


def read_rin(path: str) -> dict:
    with open(path, 'rb') as f:
        return decode(f.read())
