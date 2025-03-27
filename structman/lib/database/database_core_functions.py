import time
import sys
import traceback
import random
import math
from structman.lib.sdsc.sdsc_utils import total_size
#from structman.base_utils.config_class import Config
from typing import TypeAlias
#from jenkspy import JenksNaturalBreaks

Bin: TypeAlias = tuple[set[int], int, int]

def binningSelect(
        keys: set[int] | list[int],
        rows: list[str],
        table: str,
        config: any,
        density: float =0.5) -> list[list]:
    
    # Use this for huge selects, the first entry of rows has to be the key by which the entries are selected
    key_name: str = rows[0]

    t0 = time.time()

    singletons: list[int]
    bins: list[Bin]

    #singletons, bins = median_focus_binning(keys, density_thresh = density)
    #singletons, bins = jnb_binning(keys, density_thresh = density)
    singletons, bins = naive_neighbor_binning(keys, density_thresh = density)


    if config.verbosity >= 6:
        print('\nbinningSelect keys:\n', keys, 'and binning results:\n', singletons, '\n', bins)

    t1 = time.time()
    if config.verbosity >= 3:
        print(f'Time for binning in binningSelect {table}: {t1 - t0} {len(keys)=}')

    if len(singletons) > 0:
        if len(singletons) == 1:
            equals_rows = {key_name: singletons[0]}
            total_results = list(select(config, rows, table, equals_rows=equals_rows))
        else:
            in_rows = {key_name: singletons}
            total_results = list(select(config, rows, table, in_rows=in_rows))
    else:
        total_results = []

    t2 = time.time()
    if config.verbosity >= 3:
        print('Time for singleton select in binningSelect:', t2 - t1, 'Amount of singletons:', len(singletons))

    t3 = 0.
    t4 = 0.
    t5 = 0.
    for ids, min_id, max_id in bins:
        if config.verbosity >= 3:
            print(f'Binning select from {min_id} to {max_id}, {len(ids)=}')
        t3 += time.time()
        between_rows = {key_name: (min_id, max_id)}

        results = select(config, rows, table, between_rows=between_rows)

        t4 += time.time()
        if isinstance(ids, list):
            ids = set(ids)
        for row in results:
            if not row[0] in ids:
                continue
            total_results.append(row)
        t5 += time.time()

    if config.verbosity >= 3:
        print('Time for between selects in binningSelect:', t4 - t3, 'Amount of bins:', len(bins))
        print('Time for id checks in binningSelect:', t5 - t4)

    return total_results


def insert(table: str, columns: list[str], values: list | tuple, config: any, n_trials: int = 3, mapping_db: bool = False, db_name: str | None = None):
    if config.db_address != '-':
        wildcard_symbol = '%s'
        insert_statement = 'INSERT IGNORE'
        max_number_of_params = math.inf
    else:
        wildcard_symbol = '?'
        insert_statement = 'INSERT OR IGNORE'
        max_number_of_params = 32000

    params = []

    columns_str = ','.join(columns)

    parts = []
    value_strs = []

    if len(values) == 0:
        if config.verbosity >= 3:
            print('Empty insert try', table)
        return

    max_package_size = config.max_package_size
    size_of_values = size_estimation(values)

    number_of_packages = (size_of_values // max_package_size)
    if not size_of_values % max_package_size == 0:
        number_of_packages += 1

    package_length = (len(values) // number_of_packages)
    if not len(values) % number_of_packages == 0:
        package_length += 1

    if config.verbosity >= 3:
        print('Insert to', table, 'with total estimated size', size_of_values / 1024. / 1024., 'Mb,since max size is',
              max_package_size / 1024. / 1024., 'Mb this makes', number_of_packages,
              'packages in total (size per package:', size_of_values / 1024. / 1024. / number_of_packages, ')')

    n = 0
    for value_list in values:
        for value in value_list:
            params.append(value)
        value_str_part = '(%s)' % ','.join([wildcard_symbol] * len(value_list))
        value_strs.append(value_str_part)
        n += 1
        if n == package_length or len(params) >= max_number_of_params:
            n = 0
            value_str = ','.join(value_strs)
            parts.append((value_str, params))
            value_strs = []
            params = []
    if value_strs != []:
        value_str = ','.join(value_strs)
        parts.append((value_str, params))

    for value_str, params in parts:
        statement = f'{insert_statement} INTO {table} ({columns_str}) VALUES {value_str}'

        if config.verbosity >= 2:
            print(f'Insert with {len(params)} parameters, {db_name=}')

        if config.verbosity >= 5:
            print(f'Better size estimation of the package: {sys.getsizeof(str(params))} + {sys.getsizeof(statement)}')

        n = 0
        while n < n_trials:  # Repeat the querry if fails for n_trials times
            db, cursor = config.getDB(mapping_db = mapping_db, db_name = db_name)
            try:
                cursor.execute(statement, params)
                db.commit()
                db.close()
                break
            except:
                db.close()
                n += 1
                if n == 1:
                    [e, f, g] = sys.exc_info()
                    g = traceback.format_exc()
        if n == n_trials:
            raise NameError('Invalid Insert: %s\nParam size:%s\n%s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:500]), e, str(f), g))


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


def update(config, table, columns, values, mapping_db = False):
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

    if config.db_address != '-':

        max_package_size = config.max_package_size
        size_of_values = size_estimation(values, buffer_factor=2)

        number_of_packages = (size_of_values // max_package_size)
        if not size_of_values % max_package_size == 0:
            number_of_packages += 1

        package_length = (len(values) // number_of_packages)
        if not len(values) % number_of_packages == 0:
            package_length += 1

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
            db, cursor = config.getDB(mapping_db = mapping_db)
            try:
                cursor.execute(statement, params)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                exact_size_of_values = size_estimation(values, buffer_factor=2, exact=True)
                exact_size_of_params = size_estimation(params, buffer_factor=2, exact=True)
                raise NameError(f'Invalid Update {mapping_db=} {size_of_values=} {exact_size_of_values=} {exact_size_of_params=} {max_package_size=} {len(parts)=}:\n{statement[:300]} ... {conflict_statement} {update_str}\nParam size:{len(params)}\n{e}\n{f}\n{g}')
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

def naive_neighbor_binning(id_set: set[int] | list[int], density_thresh: float =0.5) -> tuple[list[int], list[Bin]]:
    # small sets are returned as a list of singletons
    if len(id_set) < 10:
        return list(id_set), []

    sorted_ids = sorted(id_set)
    density: float = calc_density(sorted_ids)
    if density > density_thresh:
        return [], [(id_set, sorted_ids[0], sorted_ids[-1])]

    singletons = []
    bins = []

    current_bin = []
    binning_active = False
    minimal_bin_size = len(id_set) // 50

    for pos, val in enumerate(sorted_ids[1:], start=1):
        prev_val = sorted_ids[pos - 1]

        d = val-prev_val
        if d <= 2:
            current_bin.append(prev_val)
            binning_active = True
        elif binning_active:
            current_bin.append(prev_val)
            binning_active = False
            if len(current_bin) >= minimal_bin_size:
                bins.append((current_bin, current_bin[0], current_bin[-1]))
            else:
                singletons.extend(current_bin)
            current_bin = []
        else:
            singletons.append(prev_val)
        

    if d <= 2:
        current_bin.append(val)
        bins.append((current_bin, current_bin[0], current_bin[-1]))
    else:
        singletons.append(val)

    return singletons, bins


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

def select(
        config: any,
        rows: list[str],
        table: str,
        between_rows: dict[str, tuple[int, int]] = {},
        in_rows: dict[str, list] = {},
        equals_rows: dict[str, any] = {},
        null_columns: list[str] = [],
        n_trials: int = 3,
        from_mapping_db: bool = False) -> list[list]:
    
    if len(rows) == 0:
        row_str = '*'
    else:
        row_str = ','.join(rows)

    if config.db_address != '-':
        wildcard_symbol = '%s'
    else:
        wildcard_symbol = '?'

    params: list[any] = []

    if len(between_rows) == 0 and len(in_rows) == 0 and len(equals_rows) == 0:
        where_str = ''
    else:
        where_parts: list[str] = []

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

    statement = f'SELECT {row_str} FROM {table}{where_str}'

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
        raise NameError('Invalid Select: %s\nParam size:%s\n%s\n%s, %s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:50]), from_mapping_db, config.mapping_db, e, str(f), g))
    return results