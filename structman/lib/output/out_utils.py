import time
from ast import literal_eval

from structman.lib.output import out_generator
from structman.lib.output import plots

type_factor = {
                'substitution' : 0.5,
                'deletion'     : 1.0
              }

function_factor = {
                    'Protein interaction' : 1.0,
                    'DNA interaction' : 1.0,
                    'RNA interaction' : 1.0,
                    'ligand interaction'  : 0.5,
                    'Peptide interaction' : 0.5,
                    'metal interaction'   : 0.2,
                    'ion interaction'     : 0.2,
                    'Surface'             : 0.0,
                    'Buried'              : 0.4,
                    'Core'                : 0.8,
                    '-'                   : 0.0,
                    None                  : 0.0,
                  }

def calculatePositionScore(mutation_type, function_annotation):
    pos_score = type_factor[mutation_type] * function_factor[function_annotation]
    return pos_score

def generate_ppi_filename(outfolder, session_name):
    prot_prot_file = f'{outfolder}/{session_name}.protein_protein_interactions.tsv'
    return prot_prot_file

def generate_gis_filename(outfolder, session_name):
    gi_score_file = f'{outfolder}/{session_name}_gene_isoform_scores.tsv'
    return gi_score_file

def parse_indel_analysis_file(filename, query_tags = []):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    indel_map = {'sub':{}, 'del':{}, 'ins':{}}

    indel_tag_map = {'sub':{}, 'del':{}, 'ins':{}}

    for indel_type in indel_tag_map:
        for tag in query_tags:
            indel_tag_map[indel_type][tag] = {}

    for pos, column_name in enumerate(lines[0].split('\t')):
        if column_name == 'Tags':
            tag_pos = pos
        if column_name == 'Indel':
            indel_pos = pos
        if column_name == 'Wildtype RIN-based simple classification':
            wt_class_pos = pos
        if column_name == 'Mutant RIN-based simple classification':
            mut_class_pos = pos

    for line in lines[1:]:
        words = line.replace('\n', '').split('\t')

        indel_notation = words[indel_pos]
        if indel_notation.count('delins') == 1:
            indel_type = 'sub'
        elif indel_notation.count('del') == 1:
            indel_type = 'del'
        elif indel_notation.count('ins') == 1:
            indel_type = 'ins'
        else:
            print(f'Unknown indel type for: {indel_notation}')
            continue
        tags = words[tag_pos].split(',')
        wt_class = words[wt_class_pos]
        mut_class = words[mut_class_pos]

        

        if indel_type == 'ins':
            important_class = mut_class
        else:
            important_class = wt_class

        if not important_class in indel_map[indel_type]:
            indel_map[indel_type][important_class] = 0
        indel_map[indel_type][important_class] += 1

        for tag in tags:
            if not tag in indel_tag_map[indel_type]:
                continue
            if not important_class in indel_tag_map[indel_type][tag]:
                indel_tag_map[indel_type][tag][important_class] = 0
            indel_tag_map[indel_type][tag][important_class] += 1
    return indel_map, indel_tag_map

def parse_annotation_file(annotationfile, by_conf=False, rin_classes=False):
    f = open(annotationfile, 'r')
    lines = f.readlines()
    f.close()

    tag_map = {}
    simple_tag_map = {}
    tag_sizes = {}

    class_map = {}
    simple_class_map = {}
    size = 0

    simple_high_confidence_map = {}
    hc_size = 0
    hc_threshs = [0.0, 0.2, 0.4, 0.6, 0.8]

    violins = {}
    comp_violins = {}

    tag_pos = None
    class_pos = None
    simple_class_pos = None
    confidence_pos = None
    rin_class_pos = None
    rin_simple_class_pos = None

    for pos, column_name in enumerate(lines[0].split('\t')):
        if column_name == 'Tags':
            tag_pos = pos
        if column_name == 'Class':
            class_pos = pos
        if column_name == 'Simple Class':
            simple_class_pos = pos
        if column_name == 'Confidence Value':
            confidence_pos = pos
        if column_name == 'RIN Class':
            rin_class_pos = pos
        if column_name == 'RIN Simple Class':
            rin_simple_class_pos = pos

    for line in lines[1:]:
        words = line.replace('\n', '').split('\t')

        tags = words[tag_pos].split(',')
        if not rin_classes:
            classification = words[class_pos]
            simple_classification = words[simple_class_pos]
        else:
            classification = words[rin_class_pos]
            simple_classification = words[rin_simple_class_pos]
        if words[confidence_pos] != out_generator.OutputGenerator.null_symbol:
            confidence = float(words[confidence_pos])
        else:
            confidence = 0.

        size += 1

        if classification not in class_map:
            class_map[classification] = 1
        else:
            class_map[classification] += 1

        if simple_classification not in simple_class_map:
            simple_class_map[simple_classification] = 1
        else:
            simple_class_map[simple_classification] += 1



        for tag in tags:
            if tag == '':
                continue
            if tag[0] == '#':
                if tag.count(':') == 1:
                    violin_tag, violin_value = tag[1:].split(':')
                else:
                    violin_tag, violin_value = tag[1:].split('=')
                try:
                    violin_value = float(violin_value)
                except:
                    continue

                if violin_tag not in violins:
                    violins[violin_tag] = {}

                if simple_classification not in violins[violin_tag]:
                    violins[violin_tag][simple_classification] = []
                violins[violin_tag][simple_classification].append(violin_value)

                if violin_tag not in comp_violins:
                    comp_violins[violin_tag] = {}

                if classification not in comp_violins[violin_tag]:
                    comp_violins[violin_tag][classification] = []
                comp_violins[violin_tag][classification].append(violin_value)

                continue

            if tag not in tag_map:
                tag_map[tag] = {}
                simple_tag_map[tag] = {}
                tag_sizes[tag] = 0
            if classification not in tag_map[tag]:
                tag_map[tag][classification] = 1
            else:
                tag_map[tag][classification] += 1
            if simple_classification not in simple_tag_map[tag]:
                simple_tag_map[tag][simple_classification] = 1
            else:
                simple_tag_map[tag][simple_classification] += 1
            tag_sizes[tag] += 1

        if by_conf:  # TODO
            for hc_thresh in hc_threshs:
                if confidence > hc_thresh:
                    if simple_classification not in simple_high_confidence_map:
                        simple_high_confidence_map[simple_classification] = 1
                    else:
                        simple_high_confidence_map[simple_classification] += 1
                    hc_size += 1

    return tag_map, simple_tag_map, tag_sizes, class_map, simple_class_map, simple_high_confidence_map, violins, comp_violins, size

def classDistributionFromFile(annotationfile, outfolder, session_name, config, by_conf=False, rin_classes=False):
    #"Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"
    outfile = '%s/%s' % (outfolder, session_name)

    print(outfile)

    t0 = time.time()

    tag_map, simple_tag_map, tag_sizes, class_map, simple_class_map, simple_high_confidence_map, violins, comp_violins, size = parse_annotation_file(annotationfile, by_conf=by_conf, rin_classes=rin_classes)

    t1 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part1: %s' % str(t1 - t0)))

    plots.makeViolins(violins, outfile, session_name)
    plots.makeViolins(comp_violins, outfile, session_name, add='_complex_classes')

    t2 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part2: %s' % str(t2 - t1)))

    classes = list(class_map.keys())
    outlines = ['Tag\t%s' % '\t'.join(classes)]

    words = ['total']
    for classification in classes:
        r = float(class_map[classification]) / float(size)
        words.append(str(r))
    outlines.append('\t'.join(words))

    if len(tag_map) > 1:
        for tag in tag_map:
            words = [tag]
            for classification in classes:
                if classification not in tag_map[tag]:
                    r = 0.0
                else:
                    r = float(tag_map[tag][classification]) / float(tag_sizes[tag])
                words.append(str(r))
            outlines.append('\t'.join(words))

    simple_classes = list(simple_class_map.keys())
    simple_outlines = ['Tag\t%s' % '\t'.join(simple_classes)]

    words = ['total']
    for classification in simple_classes:
        r = float(simple_class_map[classification]) / float(size)
        words.append(str(r))
    simple_outlines.append('\t'.join(words))

    if len(simple_tag_map) > 1:
        for tag in simple_tag_map:
            words = [tag]
            for classification in simple_classes:
                if classification not in simple_tag_map[tag]:
                    r = 0.0
                else:
                    r = float(simple_tag_map[tag][classification]) / float(tag_sizes[tag])
                words.append(str(r))
            simple_outlines.append('\t'.join(words))

    if by_conf:
        simple_high_confidence_outlines = ['\t'.join(list(simple_high_confidence_map.keys()))]
        words = []
        for simple_classification in simple_high_confidence_map:
            r = float(simple_high_confidence_map[simple_classification]) / float(hc_size)
            words.append(str(r))
        simple_high_confidence_outlines.append('\t'.join(words))

    if rin_classes:
        file_name_tag = 'rin_'
    else:
        file_name_tag = ''

    f = open('%s.%sclass_distribution.tsv' % (outfile, file_name_tag), 'w')
    f.write('\n'.join(outlines))
    f.close()

    f = open('%s.%ssimple_class_distribution.tsv' % (outfile, file_name_tag), 'w')
    f.write('\n'.join(simple_outlines))
    f.close()

    if by_conf:
        f = open('%s.simple_high_confidence_class_distribution.tsv' % outfile, 'w')
        f.write('\n'.join(simple_high_confidence_outlines))
        f.close()

    t3 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part3: %s' % str(t3 - t2)))

def InteractionScoreAveragesFromFile(InteractionProfilesfile, outfile, session_name, by_tag=False):
    outfile = "%s/%s" % (outfile, session_name)
    f = open(InteractionProfilesfile, 'r')
    lines = f.readlines()
    f.close()

    min_degree = None
    max_degree = None
    min_score = None
    max_score = None

    if len(lines) == 1:
        return

    for line in lines[1:]:
        row = line.replace('\n', '').split('\t')
        Ligand_Interaction_Degree = float(row[4])
        Ligand_Interaction_Score = float(row[5])
        Chain_Interaction_Degree = float(row[6])
        Chain_Interaction_Score = float(row[7])
        Short_Interaction_Degree = float(row[8])
        Short_Interaction_Score = float(row[9])
        Medium_Interaction_Degree = float(row[10])
        Medium_Interaction_Score = float(row[11])
        Long_Interaction_Degree = float(row[12])
        Long_Interaction_Score = float(row[13])

        degrees = [Ligand_Interaction_Degree, Chain_Interaction_Degree, Short_Interaction_Degree, Medium_Interaction_Degree, Long_Interaction_Degree]
        scores = [Ligand_Interaction_Score, Chain_Interaction_Score, Short_Interaction_Score, Medium_Interaction_Score, Long_Interaction_Score]
        for degree in degrees:
            if min_degree is None or degree < min_degree:
                min_degree = degree
            if max_degree is None or degree > max_degree:
                max_degree = degree

        for score in scores:
            if min_score is None or score < min_score:
                min_score = score
            if max_score is None or score > max_score:
                max_score = score

    bins = 50
    degree_bin_size = (max_degree - min_degree) / float(bins)
    score_bin_size = (max_score - min_score) / float(bins)

    degree_tag_map = {}
    score_tag_map = {}
    degree_tag_sizes = {}
    score_tag_sizes = {}

    scatter_plots = {}

    for line in lines[1:]:
        row = line.replace('\n', '').split('\t')
        Ligand_Interaction_Degree = float(row[4])
        Ligand_Interaction_Score = float(row[5])
        Chain_Interaction_Degree = float(row[6])
        Chain_Interaction_Score = float(row[7])
        Short_Interaction_Degree = float(row[8])
        Short_Interaction_Score = float(row[9])
        Medium_Interaction_Degree = float(row[10])
        Medium_Interaction_Score = float(row[11])
        Long_Interaction_Degree = float(row[12])
        Long_Interaction_Score = float(row[13])
        Class = row[14]
        comp_class = row[15]

        degrees = {'LI degree': Ligand_Interaction_Degree, 'CI degree': Chain_Interaction_Degree, 'SI degree': Short_Interaction_Degree, 'MI degree': Medium_Interaction_Degree, 'LoI degree': Long_Interaction_Degree}
        scores = {'LI score': Ligand_Interaction_Score, 'CI score': Chain_Interaction_Score, 'SI score': Short_Interaction_Score, 'MI score': Medium_Interaction_Score, 'LoI score': Long_Interaction_Score}

        if by_tag:
            tags = row[3]
            tags = '%s,All' % tags
            for tag in tags.split(','):
                if tag[0] == '#':
                    scatter_tag, scatter_value = tag[1:].split(':')
                    scatter_value = float(scatter_value)

                    if scatter_tag not in scatter_plots:
                        scatter_plots[scatter_tag] = {}

                        for score in scores:
                            scatter_plots[scatter_tag][score] = {}

                        for degree in degrees:
                            scatter_plots[scatter_tag][degree] = {}

                    for score in scores:
                        if Class not in scatter_plots[scatter_tag][score]:
                            scatter_plots[scatter_tag][score][Class] = [[], []]
                    for score in scores:
                        sc = scores[score]
                        scatter_plots[scatter_tag][score][Class][0].append(sc)
                        scatter_plots[scatter_tag][score][Class][1].append(scatter_value)

                    for degree in degrees:
                        if Class not in scatter_plots[scatter_tag][degree]:
                            scatter_plots[scatter_tag][degree][Class] = [[], []]
                    for degree in degrees:
                        deg = degrees[degree]
                        scatter_plots[scatter_tag][degree][Class][0].append(deg)
                        scatter_plots[scatter_tag][degree][Class][1].append(scatter_value)

                    continue
                for degree_name in degrees:
                    combined_tag = "%s %s" % (tag, degree_name)
                    if combined_tag not in degree_tag_map:
                        degree_tag_map[combined_tag] = [0] * bins
                    degree = degrees[degree_name]
                    for i in range(1, bins + 1):
                        if degree < min_degree + i * degree_bin_size:
                            degree_tag_map[combined_tag][i - 1] += 1
                            break
                    if combined_tag not in degree_tag_sizes:
                        degree_tag_sizes[combined_tag] = 1
                    else:
                        degree_tag_sizes[combined_tag] += 1
                for score_name in scores:
                    combined_tag = "%s %s" % (tag, score_name)
                    if combined_tag not in score_tag_map:
                        score_tag_map[combined_tag] = [0] * bins
                    score = scores[score_name]
                    for i in range(1, bins + 1):
                        if score < min_score + i * score_bin_size:
                            score_tag_map[combined_tag][i - 1] += 1
                            break
                    if combined_tag not in score_tag_sizes:
                        score_tag_sizes[combined_tag] = 1
                    else:
                        score_tag_sizes[combined_tag] += 1
        else:
            pass  # TODO

    plotScatter(scatter_plots, outfile, session_name)

    startline_words = ['']
    for i in range(0, bins):
        startline_words.append('>%s' % str(min_degree + i * degree_bin_size))

    outlines = ['\t'.join(startline_words)]
    for tag in degree_tag_map:
        words = [tag]
        for amount in degree_tag_map[tag]:
            r = amount / float(degree_tag_sizes[tag])
            words.append(str(r))
        outlines.append('\t'.join(words))

    startline_words = ['']
    for i in range(0, bins):
        startline_words.append('>%s' % str(min_score + i * score_bin_size))

    outlines.append('\t'.join(startline_words))
    for tag in score_tag_map:
        words = [tag]
        for amount in score_tag_map[tag]:
            r = amount / float(score_tag_sizes[tag])
            words.append(str(r))
        outlines.append('\t'.join(words))

    f = open('%s.interaction_profile_means.tsv' % outfile, 'w')
    f.write('\n'.join(outlines))
    f.close()
