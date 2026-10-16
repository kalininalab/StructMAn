import time
import sys
import traceback
import random
import math
import itertools
import numpy as np
from numba import njit, int64
from numba.typed import List

from filelock import FileLock

#from structman.base_utils.config_class import Config
from structman.base_utils.base_utils import add_to_times, aggregate_times
from typing import TypeAlias
#from jenkspy import JenksNaturalBreaks
from collections import deque
from sys import getsizeof, stderr

Bin: TypeAlias = tuple[set[int], int, int]

# Taken from https://code.activestate.com/recipes/577504/
def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    def dict_handler(d):
        return itertools.chain.from_iterable(d.items())

    all_handlers = {
        tuple: iter,
        list: iter,
        deque: iter,
        dict: dict_handler,
        set: iter,
        frozenset: iter,
    }
    all_handlers.update(handlers)  # user handlers take precedence
    seen = set()  # track which object id's have already been seen
    default_size = getsizeof(0)  # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:  # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)

def binningSelect(
        keys: set[int] | list[int],
        rows: list[str],
        table: str,
        config: any,
        density: float =0.5,
        locked: bool = False,
        db_lock = None,
        return_run_times = False,
        full_select = False,
        picky_select = False,
        ) -> list[list]:
    
    # Use this for huge selects, the first entry of rows has to be the key by which the entries are selected
    key_name: str = rows[0]

    times = []
    t0 = time.time()

    singletons: list[int]
    bins: list[np.ndarray[int]]

    #singletons, bins = median_focus_binning(keys, density_thresh = density)
    #singletons, bins = jnb_binning(keys, density_thresh = density)

    if config.verbosity >= 7:
        config.logger.info(f'binningSelect keys:\n{keys=}')

    singletons, bins = naive_neighbor_binning(keys, density_thresh = density, full_select=full_select, picky_select = picky_select)

    if config.verbosity >= 7:
        config.logger.info(f'binning results:\n{singletons=}\n{bins=}')

    t0 = add_to_times(times, t0, 'binning')

    if len(singletons) > 0:
        if len(singletons) == 1:
            equals_rows = {key_name: singletons[0]}
            total_results = list(select(config, rows, table, equals_rows=equals_rows, locked=locked, db_lock=db_lock))
        else:
            in_rows = {key_name: singletons}
            total_results = list(select(config, rows, table, in_rows=in_rows, locked=locked, db_lock=db_lock))
    else:
        total_results = []

    t0 = add_to_times(times, t0, 'singleton select')

    total_loop_times = []

    for ids in bins:
        loop_times = []
        tl = time.time()
        min_id = ids[0]
        max_id = ids[-1]
        if config.verbosity >= 6:
            config.logger.info(f'Binning select from {min_id} to {max_id}, {len(ids)=}')

        between_rows = {key_name: (min_id, max_id)}

        results = select(config, rows, table, between_rows=between_rows, locked=locked, db_lock=db_lock)

        tl = add_to_times(loop_times, tl, 'select')

        id_set = set(ids)

        results = [row for row in results if row[0] in id_set]
        
        total_results += results

        tl = add_to_times(loop_times, tl, 'filter results')
        
        total_loop_times = aggregate_times(total_loop_times, loop_times)


    times.append(total_loop_times)
    t0 = add_to_times(times, t0, 'bin selects')

    if return_run_times:
        return total_results, times

    return total_results


def insert(
        table: str,
        columns: list[str],
        values: list | tuple,
        config: any,
        n_trials: int = 3,
        mapping_db: bool = False,
        db_name: str | None = None,
        locked = False,
        db_lock = None,
        buffer_factor = 3.0):
    if config.db_address != '-':
        wildcard_symbol = '%s'
        insert_statement = 'INSERT IGNORE'
        max_number_of_params = math.inf
    else:
        wildcard_symbol = '?'
        insert_statement = 'INSERT OR IGNORE'
        max_number_of_params = 32000

    columns_str = ','.join(columns)

    if len(values) == 0:
        if config.verbosity >= 3:
            config.logger.info(f'Empty insert try {table=}')
        return

    max_package_size = config.max_package_size
    size_of_values = size_estimation(values, buffer_factor=buffer_factor)

    number_of_packages = (size_of_values // max_package_size)
    if not size_of_values % max_package_size == 0:
        number_of_packages += 1

    package_length = (len(values) // number_of_packages)
    if not len(values) % number_of_packages == 0:
        package_length += 1

    if config.verbosity >= 6:
        config.logger.info(f'Insert to {table} with total estimated size {size_of_values / 1024. / 1024.}Mb, since max size is {max_package_size / 1024. / 1024.}Mb this makes {number_of_packages} packages in total (size per package: {size_of_values / 1024. / 1024. / number_of_packages})')

    def iter_packages():
        # Yields the packages one at a time instead of collecting them all first.
        # For inserts of tens of millions of rows, holding every package at once
        # duplicates the entire value set as parameter lists plus one large
        # placeholder string each, on top of an already very large process.
        # The generator is re-created per retry attempt, so retries still work.
        placeholder_cache = {}
        params = []
        value_strs = []
        n = 0
        for value_list in values:
            arity = len(value_list)
            row_placeholder = placeholder_cache.get(arity)
            if row_placeholder is None:
                # This is constant per row width; it used to be rebuilt from a
                # temporary list and a join for every single row.
                row_placeholder = '(%s)' % ','.join([wildcard_symbol] * arity)
                placeholder_cache[arity] = row_placeholder

            params.extend(value_list)
            value_strs.append(row_placeholder)
            n += 1
            if n == package_length or len(params) >= max_number_of_params:
                n = 0
                yield ','.join(value_strs), params
                value_strs = []
                params = []
        if value_strs != []:
            yield ','.join(value_strs), params

    if db_lock is None:
        db_lock_file = f'{config.tmp_folder}/insert_{table}_lock.lock'
    else:
        db_lock_file = f'{config.tmp_folder}/{db_lock}'

    statement = ''
    params = []

    def run_insert():
        nonlocal statement, params
        db, cursor = config.getDB(mapping_db = mapping_db, db_name = db_name)
        # getDB returns (None, None) once its own retries are exhausted; using
        # that as a cursor buries the real cause under an AttributeError.
        if cursor is None:
            if db is not None:
                db.close()
            raise ConnectionError(
                f'No database connection for insert into {table} '
                f'({config.db_address}, {db_name if db_name is not None else config.db_name})'
            )
        try:
            for value_str, params in iter_packages():
                statement = f'{insert_statement} INTO {table} ({columns_str}) VALUES {value_str}'

                if config.verbosity >= 6:
                    config.logger.info(f'Insert with {len(params)} parameters, {db_name=}')

                if config.verbosity >= 7:
                    config.logger.info(f'Better size estimation of the package: {sys.getsizeof(str(params))} + {sys.getsizeof(statement)}')

                cursor.execute(statement, params)
                db.commit()
        finally:
            db.close()

    n = 0
    while n < n_trials:  # Repeat the querry if fails for n_trials times
        try:
            if locked:
                with FileLock(db_lock_file):
                    run_insert()
            else:
                run_insert()
            break
        except:
            n += 1
            if n == 1:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
            if n < n_trials:
                time.sleep(0.01 * (2 ** n))
    if n == n_trials:
        raise NameError(f'Invalid Insert: f{statement[:500]=}\n{len(params)=}\n{e}\n{f}\n{g}')


def remove(config, table, between_rows={}, in_rows={}, equals_rows={}, null_columns=set(), n_trials=3, from_mapping_db=False):
    if config.db_address != '-':
        wildcard_symbol = '%s'
    else:
        wildcard_symbol = '?'

    params = []

    if len(between_rows) == 0 and len(in_rows) == 0 and len(equals_rows) == 0:
        where_str = ''
    else:
        where_parts = []

        for equals_row in equals_rows:
            params.append(equals_rows[equals_row])
            where_parts.append(equals_row + f' = {wildcard_symbol}')

        for null_column in null_columns:
            where_parts.append(null_column + ' IS NULL')

        for in_row in in_rows:
            for param in in_rows[in_row]:
                params.append(param)
            where_parts.append(in_row + ' IN (%s)' % ','.join([wildcard_symbol] * len(in_rows[in_row])))  # There have to be as many %s placeholders in the statement as there are parameters for the IN clasue

        for bet_row in between_rows:
            (low, high) = between_rows[bet_row]
            params.append(low)
            params.append(high)
            where_parts.append(bet_row + f' BETWEEN {wildcard_symbol} AND {wildcard_symbol}')

        where_str = ' WHERE %s' % ' AND '.join(where_parts)

        if len(params) == 0:
            return []

    statement = 'DELETE FROM %s%s' % (table, where_str)

    n = 0
    while n < n_trials:  # Repeat the querry if fails for n_trials times
        db, cursor = config.getDB(mapping_db=from_mapping_db)
        try:
            cursor.execute(statement, params)
            results = cursor.fetchall()
            db.commit()
            break
        except:
            if n == 0:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
            n += 1
        db.close()
    if n == n_trials:
        raise NameError('Invalid Delete: %s\nParam size:%s\n%s\n%s, %s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:50]), from_mapping_db, config.mapping_db, e, str(f), g))
    return results


def set_null(config, table: str, column: str, equals_rows: dict = {}, in_rows: dict = {},
             n_trials: int = 3, locked: bool = False, db_lock = None) -> None:
    """UPDATE table SET column = NULL on the rows that match all the given filters.

    A plain UPDATE rather than update(): that one writes through INSERT ... ON
    DUPLICATE KEY UPDATE, which adds rows instead of changing them on a table without
    a unique key over the filter columns, such as the Alignment tables. Filters as in
    remove(). A call without filter values does nothing, so it can never touch a whole
    table by accident. The lock is insert()'s for the table, so this cannot interleave
    with writes to it.
    """
    if config.db_address != '-':
        wildcard_symbol = '%s'
    else:
        wildcard_symbol = '?'

    params = []
    where_parts = []
    for equals_row in equals_rows:
        params.append(equals_rows[equals_row])
        where_parts.append(f'{equals_row} = {wildcard_symbol}')
    for in_row in in_rows:
        if len(in_rows[in_row]) == 0:
            return
        params.extend(in_rows[in_row])
        where_parts.append(f'{in_row} IN ({",".join([wildcard_symbol] * len(in_rows[in_row]))})')

    if len(params) == 0:
        return

    statement = f'UPDATE {table} SET {column} = NULL WHERE {" AND ".join(where_parts)}'

    if db_lock is None:
        db_lock_file = f'{config.tmp_folder}/insert_{table}_lock.lock'
    else:
        db_lock_file = f'{config.tmp_folder}/{db_lock}'

    def run_update():
        db, cursor = config.getDB()
        if cursor is None:
            if db is not None:
                db.close()
            raise ConnectionError(f'No database connection for update of {table} ({config.db_address}, {config.db_name})')
        try:
            cursor.execute(statement, params)
            db.commit()
        finally:
            db.close()

    n = 0
    while n < n_trials:  # Repeat the querry if fails for n_trials times
        try:
            if locked:
                with FileLock(db_lock_file):
                    run_update()
            else:
                run_update()
            break
        except:
            n += 1
            if n == 1:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
            if n < n_trials:
                time.sleep(0.01 * (2 ** n))
    if n == n_trials:
        raise NameError(f'Invalid Update: {statement[:500]=}\n{len(params)=}\n{e}\n{f}\n{g}')

def size_estimation(values: list[any], buffer_factor: float = 1.5, exact = False):
    # Select 100 random rows and calculate their str size
    n_of_r_values = 100
    if len(values) <= n_of_r_values or exact:
        #size_of_values = sys.getsizeof(str(values))
        size_of_values = total_size(values)
    else:
        #size_of_values = (sys.getsizeof(str(random.sample(values,n_of_r_values)))*len(values))/n_of_r_values
        size_of_values = (total_size(random.sample(values, n_of_r_values)) * len(values)) / n_of_r_values

    # Add X% to the size estimation, just to be sure
    size_of_values = size_of_values * buffer_factor
    return size_of_values


def update(config, table, columns, values, mapping_db = False, locked=False, db_lock=None, buffer_factor = 4.0):
    if config.db_address != '-':
        wildcard_symbol = '%s'
        insert_statement = 'INSERT IGNORE'
        conflict_statement = 'ON DUPLICATE KEY UPDATE'
    else:
        wildcard_symbol = '?'
        insert_statement = 'INSERT OR IGNORE'
        conflict_statement = 'ON CONFLICT DO UPDATE SET'

    # Updating with an insert statement (see: https://stackoverflow.com/questions/20255138/sql-update-multiple-records-in-one-query)
    params = []

    if config.db_address != '-':
        column_str = ','.join(columns)

        update_str_parts = []
        for column in columns:
            update_str_parts.append('%s=VALUES(%s)' % (column, column))

        update_str = ','.join(update_str_parts)
    else:
        where_str = f'{columns[0]}={wildcard_symbol}'

        columns_parts = []
        for column in columns[1:]:
            columns_parts.append(f'{column}={wildcard_symbol}')
        column_str = ','.join(columns_parts)

    parts = []
    n = 0
    value_strs = []

    if len(values) == 0:
        return

    if db_lock is None:
        db_lock_file = f'{config.tmp_folder}/update_{table}_lock.lock'
    else:
        db_lock_file = f'{config.tmp_folder}/{db_lock}'

    if config.db_address != '-':

        max_package_size = config.max_package_size
        size_of_values = size_estimation(values, buffer_factor=buffer_factor)

        number_of_packages = (size_of_values // max_package_size)
        if not size_of_values % max_package_size == 0:
            number_of_packages += 1

        package_length = (len(values) // number_of_packages)
        if not len(values) % number_of_packages == 0:
            package_length += 1

        
        if config.verbosity >= 6:
            config.logger.info(f'Update to {table} with total estimated size {size_of_values / 1024. / 1024.}Mb, since max size is {max_package_size / 1024. / 1024.}Mb this makes {number_of_packages} packages in total (size per package: {size_of_values / 1024. / 1024. / number_of_packages})')

        n = 0
        for value_list in values:
            for value in value_list:
                params.append(value)
            value_str_part = '(%s)' % ','.join([wildcard_symbol] * len(value_list))
            value_strs.append(value_str_part)
            n += 1
            if n == package_length:
                n = 0
                value_str = ','.join(value_strs)
                parts.append((value_str, params))
                value_strs = []
                params = []
        if value_strs != []:
            value_str = ','.join(value_strs)
            parts.append((value_str, params))

        for value_str, params in parts:
            if len(params) == 0:
                continue

            statement = f'{insert_statement} INTO {table} ({column_str}) VALUES {value_str} {conflict_statement} {update_str}'

            if config.verbosity >= 8:
                config.logger.info(f'In update {statement=}')

            def run_update():
                db, cursor = config.getDB(mapping_db = mapping_db)
                # getDB returns (None, None) once its own retries are exhausted;
                # using that as a cursor buries the real cause under an
                # AttributeError.
                if cursor is None:
                    if db is not None:
                        db.close()
                    raise ConnectionError(
                        f'No database connection for update on {table} ({config.db_address})'
                    )
                try:
                    cursor.execute(statement, params)
                    db.commit()
                finally:
                    db.close()

            max_number_of_tries = 3 if locked else 1
            try_count = 0
            while True:
                try:
                    if locked:
                        with FileLock(db_lock_file):
                            run_update()
                    else:
                        run_update()
                    break
                except:
                    try_count += 1
                    # This used to compare try_count against max_number_of_tries
                    # before incrementing, so the raise was unreachable and a
                    # permanently failing update was silently dropped.
                    if try_count >= max_number_of_tries:
                        [e, f, g] = sys.exc_info()
                        g = traceback.format_exc()
                        exact_size_of_values = size_estimation(values, buffer_factor=buffer_factor, exact=True)
                        exact_size_of_params = size_estimation(params, buffer_factor=buffer_factor, exact=True)
                        raise NameError(f'Invalid Update {mapping_db=} {size_of_values=} {exact_size_of_values=} {exact_size_of_params=} {max_package_size=} {len(parts)=}:\n{statement[:300]} ... {conflict_statement} {update_str}\nParam size:{len(params)}\n{e}\n{f}\n{g}')
                    time.sleep(0.1 * (2 ** try_count))
    else:
        if locked:
            with FileLock(db_lock_file):
                db, cursor = config.getDB(mapping_db = mapping_db)
                sql = "PRAGMA synchronous = OFF;"
                cursor.execute(sql)
                db.commit()
                sql = "PRAGMA journal_mode = MEMORY;"
                cursor.execute(sql)
                db.commit()
                for value_tuple in values:
                    params = list(value_tuple[1:])
                    params.append(value_tuple[0])
                    statement = f"UPDATE OR IGNORE {table} SET {column_str} WHERE {where_str}"
                    
                    try:
                        cursor.execute(statement, params)
                        db.commit()
                    except:
                        [e, f, g] = sys.exc_info()
                        g = traceback.format_exc()
                        raise NameError('Invalid Update: %s\nParam size:%s\n%s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:500]), e, f, g))
                db.close()
        else:
            db, cursor = config.getDB(mapping_db = mapping_db)
            sql = "PRAGMA synchronous = OFF;"
            cursor.execute(sql)
            db.commit()
            sql = "PRAGMA journal_mode = MEMORY;"
            cursor.execute(sql)
            db.commit()
            for value_tuple in values:
                params = list(value_tuple[1:])
                params.append(value_tuple[0])
                statement = f"UPDATE OR IGNORE {table} SET {column_str} WHERE {where_str}"
                
                try:
                    cursor.execute(statement, params)
                    db.commit()
                except:
                    [e, f, g] = sys.exc_info()
                    g = traceback.format_exc()
                    raise NameError('Invalid Update: %s\nParam size:%s\n%s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:500]), e, f, g))
            db.close()


def calc_density(sorted_unique_ints: list[int]) -> float:
    l : int = len(sorted_unique_ints)
    if l == 0:
        return 0.
    if l == 1:
        return 1.
    
    density: float = l / (1 + sorted_unique_ints[-1] - sorted_unique_ints[0])
    return density
"""
def jnb_binning(id_set: set[int] | list[int], density_thresh: float =0.5) -> tuple[list[int], list[Bin]]:
    if len(id_set) < 10:
        return list(id_set), []
    
    sorted_ids = sorted(id_set)

    # If the given set is dense enough, just return it as one interval
    density: float = calc_density(sorted_ids)
    if density > density_thresh:
        return [], [(id_set, sorted_ids[0], sorted_ids[-1])]

    print(f'{density=}')

    jnb = JenksNaturalBreaks(3)

    jnb.fit(sorted_ids)

    singletons = []
    bins = []

    for group in jnb.groups_:
        #print(group)
        d = calc_density(group)

        if d < density_thresh:
            if d < density:
                singletons.extend(group)
            else:
                rec_singletons, rec_bins = jnb_binning(group)
                singletons.extend(rec_singletons), bins.extend(rec_bins)
        else:
            bins.append((group, group[0], group[-1]))

    return singletons, bins
"""

@njit
def njit_calc_density(sorted_unique_ints: np.ndarray) -> np.float64:
    id_len : int = len(sorted_unique_ints)
    if id_len == 0:
        return np.float64(0.)
    if id_len == 1:
        return np.float64(1.)
    
    density: np.float64 = id_len / (1 + sorted_unique_ints[-1] - sorted_unique_ints[0])
    return density

@njit
def nnb(
    id_set: np.ndarray[int],
    density_thresh: float,
    minimal_bin_size: int,
    big_set_thresh: int,
    big_set_density_thresh: float,
    ) -> tuple[np.ndarray[int], list[np.ndarray[int]]]:

    bins = List.empty_list(np.int64[:])

    sorted_ids = np.unique(id_set)

    density: np.float64 = njit_calc_density(sorted_ids)

    singletons: np.ndarray[int]

    if density > density_thresh:
        bins.append(sorted_ids)
        return (np.empty(0, dtype=np.int64), bins)

    # A large id set that is merely moderately dense still belongs in one range scan.
    # Position ids are fused as pos_db_id * 100 + db_nr with db_nr in 0..19, so any
    # set of them has a density pinned near 20/101 - just under the main threshold,
    # forever. The run lengths that fall out of that (20, with gaps of 81) are far
    # below any sensible minimal_bin_size, so every id used to end up a singleton and
    # the caller asked for all of them in one IN clause. Reading 1/density times too
    # many rows in one indexed range beats hundreds of millions of point lookups,
    # as long as the over-fetch stays bounded - which is what big_set_density_thresh
    # guarantees, and why this does not fire for a handful of ids scattered over the
    # whole key space.
    if len(sorted_ids) >= big_set_thresh and density > big_set_density_thresh:
        bins.append(sorted_ids)
        return (np.empty(0, dtype=np.int64), bins)

    singletons = [np.int64(x) for x in range(0)]
    #bins: list[np.ndarray[int]] =  []

    current_bin = [np.int64(x) for x in range(0)]
    binning_active: bool = False

    pos: int
    val: int

    for pos_minus, val in enumerate(sorted_ids[1:]):
        pos = pos_minus + 1
        prev_val: int = sorted_ids[pos - 1]

        d: int = val-prev_val
        if d <= 2:
            current_bin.append(prev_val)
            binning_active = True
        elif binning_active:
            current_bin.append(prev_val)
            binning_active = False
            if len(current_bin) >= minimal_bin_size:
                bins.append(np.array(current_bin, dtype=np.int64))
            else:
                for singleton in current_bin:
                    singletons.append(singleton)
            current_bin = [np.int64(x) for x in range(0)] 
        else:
            singletons.append(prev_val)
        

    if d <= 2:
        current_bin.append(val)
        bins.append(np.array(current_bin))
    else:
        singletons.append(val)

    singletons = np.array(singletons, dtype=np.int64)

    #if len(bins) == 0:
    #    bins = [np.array([])]

    return singletons, bins

# A bin has to be worth a round trip of its own: one BETWEEN select instead of this
# many parameters in the shared IN clause. This used to be len(id_set) // 50, which
# reads as "at least 2% of the input" but scales away to nothing - at 165 million ids
# it demanded runs of 3.3 million consecutive values, so no run ever qualified and
# every id fell through to the singleton list. Capped, so it keeps its meaning for
# small inputs and stops growing past the point where a bin is obviously worthwhile.
MINIMAL_BIN_SIZE_CAP = 1_000

# Above this many ids, a single range scan with bounded over-fetch is preferred over
# per-id lookups. See the comment at the second density check in nnb.
BIG_ID_SET_THRESH = 1_000_000
BIG_ID_SET_DENSITY_THRESH = 0.1


def naive_neighbor_binning(id_set: set[int] | list[int], density_thresh: float =0.5, full_select=False, picky_select=False) -> tuple[np.ndarray[int], list[np.ndarray[int]]]:
    # small sets are returned as a list of singletons
    if len(id_set) < 10:
        return list(id_set), []

    try:
        id_set = np.array(id_set, dtype=int)
    except TypeError:
        id_set = np.array([x for x in id_set if x is not None], dtype=int)

    if full_select:
        if len(id_set) == 0:
            return [], []
        if len(id_set) == 1:
            return id_set, []
        return [], [np.unique(id_set)]

    if picky_select:
        sorted_ids = np.unique(id_set)
        density: np.float64 = njit_calc_density(sorted_ids)

        if density > density_thresh:
            return ([], [sorted_ids])
        return sorted_ids, []

    return nnb(
        id_set, density_thresh,
        min(len(id_set) // 50, MINIMAL_BIN_SIZE_CAP),
        BIG_ID_SET_THRESH, BIG_ID_SET_DENSITY_THRESH,
    )


def median_focus_binning(id_set: set[int] | list[int], density_thresh: float =0.5) -> tuple[list[int], list[Bin]]:
    # small sets are returned as a list of singletons
    if len(id_set) < 10:
        return list(id_set), []

    sorted_ids = sorted(id_set)

    # If the given set is dense enough, just return it as one interval
    density: float = len(sorted_ids) / (1 + sorted_ids[-1] - sorted_ids[0])
    if density > density_thresh:
        return [], [(id_set, sorted_ids[0], sorted_ids[-1])]

    l_quartile_pos: int = len(id_set) // 4
    r_quartile_pos: int = 3 * l_quartile_pos

    avg_quartile_dist: int = 1 + (2 * (1 + sorted_ids[r_quartile_pos] - sorted_ids[l_quartile_pos]) // len(sorted_ids))

    median_pos: int = len(id_set) // 2
    median: int = sorted_ids[median_pos]

    id_set: set[int] = set([median])
    singletons: list[int] = []

    l_minus_1_value: int = median
    r_plus_1_value: int = median  # needed for giving the boundaries of the set if only median is in the set
    l_value: int
    for i in range(median_pos):
        l_value = sorted_ids[median_pos - i]
        l_minus_1_value = sorted_ids[median_pos - i - 1]
        if (l_value - l_minus_1_value) < avg_quartile_dist:
            id_set.add(l_minus_1_value)
        else:
            l_minus_1_value = l_value  # needed for giving the boundaries of the set
            break

    for j in range(median_pos - i):
        singletons.append(sorted_ids[j])

    r_value: int
    r_plus_1_value: int
    for i in range(median_pos, len(sorted_ids) - 1):
        r_value = sorted_ids[i]
        r_plus_1_value = sorted_ids[i + 1]
        if (r_plus_1_value - r_value) < avg_quartile_dist:
            id_set.add(r_plus_1_value)
        else:
            r_plus_1_value = r_value
            break

    for j in range(i + 1, len(sorted_ids)):
        singletons.append(sorted_ids[j])

    if l_minus_1_value == r_plus_1_value:
        singletons.append(l_minus_1_value)
        return singletons, []

    return singletons, [(id_set, l_minus_1_value, r_plus_1_value)]

# Even where max_package_size would allow more, an IN list this long stops being a
# good idea: MySQLdb builds the whole escaped query in memory before sending it, and
# the server's optimizer degrades on huge IN sets. sqlite has a hard parameter limit
# (SQLITE_MAX_VARIABLE_NUMBER, 32766 on current builds) and gets its own cap.
MAX_IN_CLAUSE_PARAMS = 100_000
MAX_IN_CLAUSE_PARAMS_SQLITE = 30_000


def select_chunk_length(config: any, chunked_values: list[any], n_other_params: int) -> int:
    """How many IN values one statement may carry."""
    if len(chunked_values) == 0:
        return 1

    if config.db_address != '-':
        cap = MAX_IN_CLAUSE_PARAMS
    else:
        cap = MAX_IN_CLAUSE_PARAMS_SQLITE

    # Leave room for the parameters the other clauses already contribute; they are
    # repeated in every chunk.
    cap = max(1, cap - n_other_params)

    max_package_size = getattr(config, 'max_package_size', None)
    if max_package_size:
        # Per value: its escaped literal plus the ',%s' placeholder it needs in the
        # statement. size_estimation samples the list rather than walking all of it,
        # which matters when 'all of it' is a hundred million entries.
        bytes_per_value = (size_estimation(chunked_values, buffer_factor=1.5) / len(chunked_values)) + 3
        by_package_size = int((max_package_size * 0.9) / max(1.0, bytes_per_value))
        cap = min(cap, max(1, by_package_size))

    return cap


def select(
        config: any,
        rows: list[str],
        table: str,
        between_rows: dict[str, tuple[int, int]] = {},
        in_rows: dict[str, list] = {},
        equals_rows: dict[str, any] = {},
        null_columns: list[str] = [],
        n_trials: int = 3,
        from_mapping_db: bool = False,
        db_name: str | None = None,
        locked: bool = False,
        db_lock = None) -> list[list]:
    
    if len(rows) == 0:
        row_str = '*'
    else:
        row_str = ','.join(rows)

    if config.db_address != '-':
        wildcard_symbol = '%s'
    else:
        wildcard_symbol = '?'

    params: list[any] = []

    # The one IN list that gets split across statements, if any. insert() and
    # update() have always sized their packages against max_package_size; select()
    # did not, so a big IN list went to the server as a single statement. A whole-PDB
    # output run reached 165 million parameters in one SELECT, which MySQLdb mogrifies
    # into one bytes query - hundreds of GB in the driver - before the server drops
    # the connection with 'Server has gone away'. The longest list is the one worth
    # splitting; the others are short (multi-column in_rows only comes from
    # suggest.py) and are repeated verbatim in every chunk, which is sound because
    # `A IN (chunk) AND B IN (...)` over disjoint chunks of A unions back to exactly
    # the rows of `A IN (all) AND B IN (...)`.
    chunked_row: str | None = None
    chunked_values: list[any] = []
    if len(in_rows) > 0:
        chunked_row = max(in_rows, key=lambda in_row: len(in_rows[in_row]))
        chunked_values = in_rows[chunked_row]
        if not isinstance(chunked_values, (list, tuple)):
            # dict_keys, sets and numpy arrays all arrive here; slicing needs a
            # sequence.
            chunked_values = list(chunked_values)

    # where_parts and part_params stay index-aligned: a clause's parameters sit at
    # its own slot, so the flattened parameter order always matches the placeholder
    # order in the statement. The chunked column's slot is None in both, and the
    # chunk is spliced in at exactly that position - appending it at the end instead
    # binds the wrong values whenever another clause follows it.
    where_parts: list[str | None] = []
    part_params: list[list[any] | None] = []

    if len(between_rows) == 0 and len(in_rows) == 0 and len(equals_rows) == 0:
        where_str = ''
    else:
        for equals_row in equals_rows:
            part_params.append([equals_rows[equals_row]])
            where_parts.append(equals_row + f' = {wildcard_symbol}')

        for null_column in null_columns:
            part_params.append([])
            where_parts.append(null_column + ' IS NULL')

        for in_row in in_rows:
            if in_row == chunked_row:
                # Placeholder; the real clause is built per chunk below.
                part_params.append(None)
                where_parts.append(None)
                continue
            part_params.append(list(in_rows[in_row]))
            where_parts.append(in_row + ' IN (%s)' % ','.join([wildcard_symbol] * len(in_rows[in_row])))  # There have to be as many %s placeholders in the statement as there are parameters for the IN clasue

        for bet_row in between_rows:
            (low, high) = between_rows[bet_row]
            part_params.append([low, high])
            where_parts.append(bet_row + f' BETWEEN {wildcard_symbol} AND {wildcard_symbol}')

        params = [p for part in part_params if part is not None for p in part]

        if len(params) == 0 and len(chunked_values) == 0:
            return []

        where_str = ' WHERE %s' % ' AND '.join(
            [f'{chunked_row} IN (...)' if part is None else part for part in where_parts])

    # Only ever executed when nothing is chunked; otherwise it is the shape reported
    # if the select fails before the first chunk is built.
    statement = f'SELECT {row_str} FROM {table}{where_str}'

    chunk_length = select_chunk_length(config, chunked_values, len(params))

    if len(chunked_values) > chunk_length:
        # IN is set membership, so dropping repeats changes nothing about the result,
        # and it is what makes the chunks disjoint - without it a value appearing in
        # two chunks would return its rows twice. Only worth the extra container when
        # the select is actually about to be split; a single statement lets the
        # server collapse the repeats itself.
        chunked_values = list(dict.fromkeys(chunked_values))

    def iter_statements():
        # Yields (statement, params) one chunk at a time rather than materialising
        # every chunk first, so only one chunk's parameters are alive at once.
        if chunked_row is None:
            yield statement, params
            return

        for start in range(0, len(chunked_values), chunk_length):
            chunk = chunked_values[start:start + chunk_length]
            chunk_clause = chunked_row + ' IN (%s)' % ','.join([wildcard_symbol] * len(chunk))
            chunk_where = ' WHERE %s' % ' AND '.join(
                [chunk_clause if part is None else part for part in where_parts])
            chunk_params = [p
                            for part in part_params
                            for p in (chunk if part is None else part)]
            yield f'SELECT {row_str} FROM {table}{chunk_where}', chunk_params

    if len(chunked_values) > chunk_length and config.verbosity >= 4:
        n_chunks = 1 + (len(chunked_values) - 1) // chunk_length
        config.logger.info(
            f'Select from {table} on {chunked_row}: {len(chunked_values)} ids '
            f'split into {n_chunks} statements of up to {chunk_length}'
        )

    if db_lock is None:
        db_lock_file = f'{config.tmp_folder}/select_lock.lock'
  
    else:
        db_lock_file = f'{config.tmp_folder}/{db_lock}'


    # Reported in the failure message, so they have to survive run_select.
    failed_statement = statement
    failed_params: list[any] = params

    def run_select():
        # getDB gives back (None, None) once it has exhausted its own retries;
        # using that as a cursor raised an AttributeError that hid the actual cause
        # (the database refusing the connection) behind a bogus 'Invalid Select'.
        nonlocal failed_statement, failed_params
        db, cursor = config.getDB(mapping_db=from_mapping_db, db_name = db_name)
        if cursor is None:
            if db is not None:
                db.close()
            raise ConnectionError(
                f'No database connection for select on {table} '
                f'({config.db_address}, {db_name if db_name is not None else config.db_name})'
            )
        try:
            # One connection for every chunk, the way insert() does it. The
            # generator is re-created per retry attempt, so a retry redoes the
            # whole select rather than resuming half a result set.
            fetched = []
            for failed_statement, failed_params in iter_statements():
                cursor.execute(failed_statement, failed_params)
                fetched.extend(cursor.fetchall())
            db.commit()
            return fetched
        finally:
            db.close()

    n = 0
    while n < n_trials:  # Repeat the querry if fails for n_trials times
        broken = False
        try:
            if locked:
                with FileLock(db_lock_file):
                    results = run_select()
            else:
                results = run_select()
            broken = True
        except:
            if n == 0:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
            n += 1
            time.sleep(0.01 * (2 ** n))
        if broken:
            break
    if n == n_trials:
        raise NameError(f'Invalid Select: {failed_statement[:500]=} {len(failed_params)=}\n{failed_params[:50]=}\n{from_mapping_db=}\n{e}\n{f}\n{g}')
    return results