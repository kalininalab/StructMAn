import time
import sys
import traceback
import random
import math
from structman.lib.sdsc.sdsc_utils import total_size

def binningSelect(keys, rows, table, config, binning_function='median_focus', density=0.5):
    # Use this for huge selects, the first entry of rows has to be the key by which the entries are selected
    key_name = rows[0]
    t0 = time.time()
    if binning_function == 'split_fusion':
        singletons, bins = split_fusion_binning(keys, density=density, fusion=True)
    elif binning_function == 'split':
        singletons, bins = split_fusion_binning(keys, density=density, fusion=False)
    elif binning_function == 'median_focus':
        singletons, bins = median_focus_binning(keys)
    else:
        config.errorlog.add_error('Unknown binning function:', binning_function)
        return []

    if config.verbosity >= 6:
        print('\nbinningSelect keys:\n', keys, 'and binning results:\n', singletons, '\n', bins)

    t1 = time.time()
    if config.verbosity >= 3:
        print(f'Time for binning in binningSelect {table}: {t1 - t0}')

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
        print('Time for between select in binningSelect:', t4 - t3, 'Amount of bins:', len(bins))
        print('Time for id check in binningSelect:', t5 - t4)

    return total_results


def insert(table, columns, values, config, n_trials=3, mapping_db = False):
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
    size_of_values = size_estimation(config, values)

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
            print('Insert with', len(params), 'parameters')

        if config.verbosity >= 3:
            print('Better size estimation of the package:', sys.getsizeof(str(params)) + sys.getsizeof(statement))

        n = 0
        while n < n_trials:  # Repeat the querry if fails for n_trials times
            db, cursor = config.getDB(mapping_db = mapping_db)
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


def size_estimation(config, values):
    # Select 100 random rows and calculate their str size
    n_of_r_values = 100
    if len(values) <= n_of_r_values:
        #size_of_values = sys.getsizeof(str(values))
        size_of_values = total_size(values)
    else:
        #size_of_values = (sys.getsizeof(str(random.sample(values,n_of_r_values)))*len(values))/n_of_r_values
        size_of_values = (total_size(random.sample(values, n_of_r_values)) * len(values)) / n_of_r_values

    # Add X% too the size estimation, just to be sure
    size_of_values = size_of_values * 2.5
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
        size_of_values = size_estimation(config, values)

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


def id_set_split(id_set, min_id, max_id):
    if len(id_set) == 0:
        return []
    if len(id_set) == 1:
        return [(id_set, min_id, max_id)]
    if len(id_set) == 2:
        return [(set([min_id]), min_id, min_id), (set([max_id]), max_id, max_id)]
    split_id = (min_id + max_id) // 2
    l = set()
    min_l = min_id
    max_l = min_id
    r = set()
    min_r = max_id
    max_r = max_id
    for i in id_set:
        if i < split_id:
            l.add(i)
            if i > max_l:
                max_l = i
        else:
            r.add(i)
            if i < min_r:
                min_r = i
    return [(l, min_l, max_l), (r, min_r, max_r)]


def split_fusion_binning(id_set, tablesize=None, binsize=500000, density=0.5, fusion=True):
    if tablesize == 0:
        return set(), []

    if tablesize is not None:
        if tablesize > 0:
            factor = 1
            density = float(len(id_set)) / (float(tablesize) * factor)
        else:
            density = 0.
    else:
        factor = 1

    bins = {}
    for i in id_set:
        bin_number = i // binsize
        if bin_number not in bins:
            bins[bin_number] = set()
        bins[bin_number].add(i)
    sorted_bins = []
    for bin_number in sorted(bins.keys()):
        id_set = bins[bin_number]
        sorted_bins.append((id_set, min(id_set), max(id_set)))

    dense_enough = False
    while not dense_enough:
        dense_enough = True
        split_bins = []
        for (id_set, min_id, max_id) in sorted_bins:
            min_amount = ((1 + max_id - min_id) * density)
            # split set if smaller than min_amount
            if len(id_set) < min_amount:
                for (iset, mi, ma) in id_set_split(id_set, min_id, max_id):
                    if len(iset) > 0:
                        split_bins.append((iset, mi, ma))
                    dense_enough = False
            else:
                split_bins.append((id_set, min_id, max_id))

        sorted_bins = split_bins

    if fusion:

        fusion_done = False
        while not fusion_done:
            fusion_done = True
            fused_bins = []
            last_bin_fused = False
            for pos, (id_set, min_id, max_id) in enumerate(sorted_bins):
                if pos == 0 or last_bin_fused:
                    last_bin_fused = False
                    if pos == (len(sorted_bins) - 1):
                        fused_bins.append((id_set, min_id, max_id))
                    continue
                pre_id_set, pre_min_id, pre_max_id = sorted_bins[pos - 1]
                if ((len(id_set) + len(pre_id_set))) > ((max_id - pre_min_id) * density * 2 * factor):

                    fused_bins.append(((pre_id_set | id_set), pre_min_id, max_id))
                    last_bin_fused = True
                    fusion_done = False
                else:
                    fused_bins.append((pre_id_set, pre_min_id, pre_max_id))
                    last_bin_fused = False
                    if pos == (len(sorted_bins) - 1):
                        fused_bins.append((id_set, min_id, max_id))

            sorted_bins = fused_bins

    singletons = []
    non_singleton_bins = []
    for (id_set, min_id, max_id) in sorted_bins:
        if min_id == max_id:
            singletons.append(min_id)
        else:
            non_singleton_bins.append((id_set, min_id, max_id))
    return singletons, non_singleton_bins


def median_focus_binning(id_set, density_thresh=0.5):
    # small sets are returned as a list of singletons
    if len(id_set) < 10:
        return list(id_set), []

    sorted_ids = sorted(id_set)

    # If the given set is dense enough, just return it as one interval
    density = len(sorted_ids) / (1 + sorted_ids[-1] - sorted_ids[0])
    if density > density_thresh:
        return [], [(id_set, sorted_ids[0], sorted_ids[-1])]

    l_quartile_pos = len(id_set) // 4
    r_quartile_pos = 3 * l_quartile_pos

    avg_quartile_dist = 1 + (2 * (1 + sorted_ids[r_quartile_pos] - sorted_ids[l_quartile_pos]) / len(sorted_ids))

    median_pos = len(id_set) // 2
    median = sorted_ids[median_pos]

    id_set = set([median])
    singletons = []

    l_minus_1_value = median
    r_plus_1_value = median  # needed for giving the boundaries of the set if only median is in the set
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

def select(config, rows, table, between_rows={}, in_rows={}, equals_rows={}, null_columns=set(), n_trials=3, from_mapping_db=False):
    if len(rows) == 0:
        row_str = '*'
    else:
        row_str = ','.join(rows)

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

    statement = 'SELECT %s FROM %s%s' % (row_str, table, where_str)

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