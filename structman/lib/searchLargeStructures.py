import os

import psutil
import ray

from structman import settings
from structman.lib import pdbParser, serializedPipeline
from structman.base_utils import ray_utils


@ray.remote(num_cpus=1)
def parse(parse_dump, sub_folder):
    ray_utils.ray_hack()

    divide_folder, in_set, chain_limit, config = parse_dump
    pdb_path = config.pdb_path

    files = os.listdir("%s%s" % (divide_folder, sub_folder))

    outlines = []
    bio_entries = set()
    for filename in files:
        if not filename.count('.pdb1.gz') > 0:
            continue
        pdb_id = filename[:4].upper()

        if pdb_id not in in_set:
            continue

        bio_entries.add(pdb_id)

        parse_out = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path)

        if parse_out is None:
            print(pdb_id, 'pdbParser failed')

        (template_page, interaction_partners, chain_type_map, oligo, rare_residues) = parse_out

        N = 0
        for chain in chain_type_map:
            if chain_type_map[chain] != 'Protein':
                continue
            example_chain = chain
            N += 1
        if N > chain_limit:
            outlines.append('%s:%s\n' % (pdb_id, example_chain))

    return outlines, bio_entries


@ray.remote(num_cpus=1)
def parseAU(sub_folder, AU_dump):
    # hack proposed by the devs of ray to prevent too many processes being spawned
    # resources = ray.ray.get_resource_ids()
    # cpus = [v[0] for v in resources['CPU']]
    # psutil.Process().cpu_affinity(cpus)

    AU_folder, bio_entries, in_set, chain_limit, config = AU_dump
    pdb_path = config.pdb_path

    files = os.listdir("%s%s" % (AU_folder, sub_folder))
    outlines = []

    for filename in files:
        if not filename.count('.ent.gz') > 0:
            continue
        pdb_id = filename[3:7].upper()
        if pdb_id in bio_entries:
            continue

        if pdb_id not in in_set:
            continue

        pdb_id = '%s_AU' % pdb_id

        parse_out = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path)

        if parse_out is None:
            print(pdb_id, 'pdbParser failed')

        (template_page, interaction_partners, chain_type_map, oligo, rare_residus) = parse_out

        N = 0
        for chain in chain_type_map:
            if chain_type_map[chain] != 'Protein':
                continue
            example_chain = chain
            N += 1
        if N > chain_limit:
            outlines.append('%s:%s\n' % (pdb_id, example_chain))
    return outlines


def search(config, chain_limit, infile=None):

    complexes = set()
    if infile is not None:
        f = open(infile, 'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            words = line.replace('\t', ' ').split()
            pdb_tuple = words[0]
            if pdb_tuple.count(':') != 1:
                continue
            pdb_id = pdb_tuple.split(':')[0]
            complexes.add(pdb_id)

    outlines = []
    pdb_path = config.pdb_path

    ray_utils.ray_init(config)

    divide_folder = "%s/data/biounit/PDB/divided/" % pdb_path

    sub_folders = os.listdir(divide_folder)

    parse_results = []
    parse_dump = ray.put((divide_folder, complexes, chain_limit, config))

    for sub_folder in sub_folders:
        parse_results.append(parse.remote(parse_dump, sub_folder))

    bio_entries = set()
    parse_out = ray.get(parse_results)
    for outlines_part, bio_entries_part in parse_out:
        outlines += outlines_part
        bio_entries = bio_entries | bio_entries_part

    AU_folder = "%s/data/structures/divided/pdb/" % pdb_path

    AU_sub_folders = os.listdir(AU_folder)

    sub_folders = sorted(AU_sub_folders)

    AU_dump = ray.put((AU_folder, bio_entries, complexes, chain_limit, config))
    parse_results = []

    for sub_folder in sub_folders:
        parse_results.append(parseAU.remote(sub_folder, AU_dump))

    parse_out = ray.get(parse_results)
    for outlines_part in parse_out:
        outlines += outlines_part

    ray.shutdown()

    f = open('large_structures_more_than_%s_chains.smlf' % str(chain_limit), 'w')
    f.write(''.join(outlines))
    f.close()

    return
