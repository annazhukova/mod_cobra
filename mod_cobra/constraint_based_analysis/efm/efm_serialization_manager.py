from collections import defaultdict
import logging

import libsbml
import openpyxl
from mod_cobra.constraint_based_analysis.efm.control_effective_flux_calculator import get_fm_yield

from mod_sbml.utils.misc import invert_map
from mod_cobra.constraint_based_analysis.efm.EFM import EFM
from mod_sbml.sbml.sbml_manager import get_r_comps
from mod_sbml.annotation.kegg.kegg_annotator import get_kegg_r_id
from mod_sbml.sbml.submodel_manager import submodel
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula
from mod_sbml.serialization.xlsx_helper import BASIC_STYLE, save_data

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'

basic_r_style = lambda r_id: BASIC_STYLE

SORT_BY_PATTERN_LENGTH = 1
SORT_BY_EFM_NUMBER = 0
SORT_BY_WEIGHTED_PRODUCT_OF_EFM_NUMBER_AND_PATTERN_LENGTH = 2


def get_pattern_sorter(id2pattern, p_id2efm_ids, sort=SORT_BY_EFM_NUMBER):
    if not id2pattern:
        return SIMPLE_PATTERN_SORTER
    if SORT_BY_PATTERN_LENGTH == sort:
        return lambda p_id: -len(id2pattern[p_id])
    elif SORT_BY_EFM_NUMBER == sort:
        return lambda p_id: -len(p_id2efm_ids[p_id])
    elif SORT_BY_WEIGHTED_PRODUCT_OF_EFM_NUMBER_AND_PATTERN_LENGTH == sort:
        max_pattern_len = max(len(pattern) for pattern in id2pattern.itervalues())
        min_pattern_len = min(len(pattern) for pattern in id2pattern.itervalues())
        max_efm_number = max(len(efm_ids) for efm_ids in p_id2efm_ids.itervalues())
        min_efm_num = min(len(efm_ids) for efm_ids in p_id2efm_ids.itervalues())
        return lambda p_id: \
            -(1.0 * (len(id2pattern[p_id]) - min_pattern_len) / (max_pattern_len - min_pattern_len)
              if max_pattern_len != min_pattern_len else 1) * \
            ((1.0 * len(p_id2efm_ids[p_id]) - min_efm_num) / (max_efm_number - min_efm_num)
             if max_efm_number != min_efm_num else 1)
    else:
        logging.error('Unknown pattern sorting function %s. Will sort by pattern id.' % sort)
        return SIMPLE_PATTERN_SORTER


def serialize_efms_txt(id2efm, path, efm_id2efficiency, all_efm_intersection):
    with open(path, 'w+') as f:
        f.write('Found %d EFMs:\n--------------------------------------------------------------\n\n' % len(id2efm))
        if all_efm_intersection and len(all_efm_intersection):
            f.write(('Intersection of all EFMs has length %d:\n\t%s.\n'
                     % (len(all_efm_intersection), all_efm_intersection.to_string(binary=True))) +
                    '--------------------------------------------------------------\n\n')
        for efm_id in sorted(id2efm.iterkeys()):
            efm = id2efm[efm_id]
            f.write('EFM %d of length %d of efficiency %g:\n\t%s\n\n'
                    % (efm_id, len(efm), efm_id2efficiency[efm_id], efm.to_string(subpattern=all_efm_intersection)))


def serialize_fms_txt(id2fm, id2efm, fm_id2key, key2efm_ids, path, efm_id2efficiency, fm_id2efficiency,
                      all_efm_intersection, model, r_id, r_rev, in_r_id, in_r_rev):
    with open(path, 'w+') as f:
        f.write('Found %d EFMs, grouped into %d pathways:\n==============================================================\n\n'
                % (len(id2efm), len(id2fm)))
        if all_efm_intersection and len(all_efm_intersection):
            f.write(('All EFMs contain the following %d reactions:\n\t%s.\n'
                     % (len(all_efm_intersection), all_efm_intersection.to_string(binary=True))) +
                    '==============================================================\n\n')
        for fm_id in sorted(id2fm.iterkeys()):
            fm = id2fm[fm_id]
            r_id2st, p_id2st = fm_id2key[fm_id]
            efm_ids = key2efm_ids[r_id2st, p_id2st]
            yield_str = '' if not in_r_id else ' of yield %g' % get_fm_yield(fm, r_id, in_r_id, r_rev, in_r_rev)
            if len(efm_ids) == 1:
                efm_id = next(iter(efm_ids))
                efm = id2efm[efm_id]

                f.write('EFM %s (a.k.a. pathway %s) of length %d of efficiency %g%s.\n\n'
                        % (efm_id, fm_id, len(efm), efm_id2efficiency[efm_id], yield_str))
            else:
                f.write('Pathway %s of length %d of efficiency %g%s.\n\n'
                        % (fm_id, len(fm), fm_id2efficiency[fm_id], yield_str))
            f.write('Inputs: %s;\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                for (m_id, st) in r_id2st))
            f.write('Outputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                 for (m_id, st) in p_id2st))
            f.write('Structure: %s\n\n' % fm.to_string(subpattern=all_efm_intersection))
            if len(efm_ids) > 1:
                f.write('Contains %d EFM%s:\n\n' % (len(efm_ids), 's' if len(efm_ids) != 1 else ''))
                for efm_id in efm_ids:
                    efm = id2efm[efm_id]
                    yield_str = '' if not in_r_id else ' of yield %g' % get_fm_yield(efm, r_id, in_r_id, r_rev, in_r_rev)
                    f.write('\tEFM %s of length %d of efficiency %g%s:\n\n\t%s\n\n\n'
                            % (efm_id, len(efm), efm_id2efficiency[efm_id], yield_str,
                               efm.to_string(subpattern=all_efm_intersection)))
            f.write('-------------------------------------------------------------\n\n')


def serialize_communities(id2cluster, id2intersection, id2imp_rns, total_len, all_fm_intersection, path):
    with open(path, 'w+') as f:
        f.write('Analysed %d pathways\n-------------------------------------------------------------------\n\n'
                % total_len)
        if all_fm_intersection and len(all_fm_intersection):
            f.write(('All pathways contain following %d reactions:\n\n\t%s.\n\n'
                     % (len(all_fm_intersection), all_fm_intersection.to_string(binary=True))) +
                    '===================================================================\n\n')
        f.write('Found %d communities\n--------------------------------------------------------------\n\n' %
                len(id2cluster))
        all_r_ids = set(all_fm_intersection.to_r_id2coeff().iterkeys())
        for clu_id in sorted(id2cluster.iterkeys()):
            cluster = id2cluster[clu_id]
            intersection = id2intersection[clu_id]
            imp_rns = id2imp_rns[clu_id]
            f.write('Community %d contains %d following pathways:\n\n\t%s\n\n'
                    % (clu_id, len(cluster), ', '.join(sorted(cluster))))
            f.write('all of which contain following reactions:\n\n\t%s\n\n'
                    % (intersection.to_string(binary=True, subpattern=all_fm_intersection)))
            f.write(('which form a reaction community, which also contains following reactions:\n\n\t %s\n\n' +
                     '--------------------------------------------------------------\n\n') %
                    imp_rns.to_string(binary=True, subpattern=intersection, show_subpattern=False,
                                     key=lambda r_id: (0 if r_id in all_r_ids else 1, r_id)))


def read_efms(output_efm_file, r_ids, rev_r_ids):
    id2efm = {}
    with open(output_efm_file, 'r') as f:
        for line in f:
            line = line.replace('\n', '').strip()
            # header line
            if line.find('Found ') == 0:
                continue
            efm = line.split('\t')
            if not efm:
                continue
            efm_id = int(efm[0])
            efm = efm[1:]
            id2efm[efm_id] = EFM(r_id2coeff={r_id: float(coeff) for (coeff, r_id) in (it.split(' ') for it in efm)},
                                 r_ids=r_ids, rev_r_ids=rev_r_ids)
    return id2efm


def r_ids2sbml(r_ids, sbml, out_sbml, suffix='', r_updater=lambda r: True):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    submodel(r_ids, model)
    model.setId('%s_%s' % (model.getId(), suffix))
    model.setName('%s_%s' % (model.getName(), suffix))
    for r_id in r_ids:
        r = model.getReaction(r_id)
        if r:
            r_updater(r)
    libsbml.SBMLWriter().writeSBMLToFile(doc, out_sbml)


def serialize_important_reactions(r_id2efm_ids, model, output_file, imp_rn_threshold=0):
    len2r_ids = defaultdict(list)
    for r_id, efm_ids in r_id2efm_ids.iteritems():
        len2r_ids[len(efm_ids)].append(r_id)

    with open(output_file, 'w+') as f:
        f.write('Found %d reactions that participate in at least %d EFMs\n-----------------------------\n\n'
                % (len(r_id2efm_ids), imp_rn_threshold))
        for num, r_ids in sorted(len2r_ids.iteritems(), key=lambda (n, _): -n):
            for r_id in sorted(r_ids):
                f.write('in %d EFMs\t%s:%s\n' %
                        (num, r_id, get_sbml_r_formula(model, model.getReaction(r_id[1:] if '-' == r_id[0] else r_id))))
            f.write('\n')


def serialize_pattern(p_id, pattern, efm_ids, model, output_file):
    """
    Serializes a pattern to a file.

    :param p_id: pattern id

    :param pattern: EFM.

    :param output_file: path to the file where the pattern should be saved

    :param efm_ids: a collection of ids of EFMs that contain this pattern
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = pattern.to_r_id2coeff(binary=True)
        f.write('Pattern %d of length %d found in %d EFMs:\n\t%s\n------------------\n\n'
                % (p_id, len(r_id2coeff), len(efm_ids), ', '.join(sorted(efm_ids))))
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
                                                                           show_metabolite_ids=True)))
            f.write('\n')


def serialize_clique(cl_id, clique, efm_ids, model, output_file):
    """
    Serializes a clique to a file.

    :param cl_id: clique id

    :param clique: EFM.

    :param output_file: path to the file where the clique should be saved
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = clique.to_r_id2coeff(binary=True)
        f.write('Reaction group %d of length %d found in %d EFMs:\n\t%s\n------------------\n\n'
                % (cl_id, len(r_id2coeff), len(efm_ids), ', '.join(sorted(efm_ids))))
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
                                                                           show_metabolite_ids=True)))
            f.write('\n')


def serialize_community(c_id, community_intersection, imp_rns, efm_ids, model, output_file):
    """
    Serializes a community intersection to a file.

    :param c_id: community id

    :param community_intersection: community intersection.

    :param output_file: path to the file where the community intersection should be saved
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = community_intersection.to_r_id2coeff(binary=True)
        imp_r_id2coeff = {r_id: coeff for (r_id, coeff) in imp_rns.to_r_id2coeff(binary=True).iteritems()
                          if r_id not in r_id2coeff}
        f.write('Community %d of %d following EFMs:\n\n\t%s\n\n------------------\n\n' % (c_id, len(efm_ids),
                                                                                          ', '.join(sorted(efm_ids))))
        f.write('Reactions that participate in all EFMs of this community (%d):\n\n------------------\n\n'
                % len(r_id2coeff))
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id),
                                                                           show_compartments=False,
                                                                           show_metabolite_ids=True)))
            f.write('\n')
        if imp_r_id2coeff:
            f.write('------------------\n\n')
            f.write('Other reactions that participate in the same reaction community (%d):\n\n------------------\n\n'
                    % len(imp_r_id2coeff))
            coeff2r_id = invert_map(imp_r_id2coeff)
            for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
                for r_id in sorted(r_ids):
                    f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id),
                                                                               show_compartments=False,
                                                                               show_metabolite_ids=True)))
                f.write('\n')


def serialize_fm(fm_id, fm, model, output_file, efficiency, r_id, in_r_id, r_rev, in_r_rev):
    """
    Serializes a FM to a file.

    :param fm_id: FM id

    :param fm: FM.

    :param output_file: path to the file where the FM should be saved
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = fm.to_r_id2coeff()
        yield_str = '' if not in_r_id else ' of yield %g' % get_fm_yield(fm, r_id, in_r_id, r_rev, in_r_rev)
        f.write('Pathway %s of length %d of efficiency %g%s\n------------------\n\n'
                % (fm_id, len(r_id2coeff), efficiency, yield_str))
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' %
                        (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
                                                         show_metabolite_ids=True)))
            f.write('\n')


def serialize_patterns(p_id2efm_ids, id2pattern, output_file, min_pattern_len, all_efm_intersection=None,
                       min_efm_num=2, sorter=SIMPLE_PATTERN_SORTER):
    """
    Serializes patterns to a file, one pattern per line. Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Patterns are sorted according to the sorter function.

    :param p_id2efm_ids: dict, {pattern id: ids of EFMs containing this pattern}

    :param id2pattern: dict, {pattern id: pattern}.

    :param output_file: path to the file where the patterns should be saved

    :param sorter: a function that given a pattern id return a value that will be used to sort a collection of patterns
    """

    def log_pattern(p_id, f):
        pattern = id2pattern[p_id]
        p_length = len(pattern)
        efm_len = len(p_id2efm_ids[p_id])
        p_string = pattern.to_string(binary=True, subpattern=all_efm_intersection)
        f.write("Pattern %d of length %d:\n\t%s,\nfound in %d EFMs:\n\t%s\n\n"
                % (p_id, p_length, p_string, efm_len, ', '.join(sorted(p_id2efm_ids[p_id]))))

    with open(output_file, 'w+') as f:
        f.write(('Found %d patterns of length >= %d, contained in at least %d EFMs.\n'
                 % (len(id2pattern), min_pattern_len, min_efm_num)) +
                '--------------------------------------------------------------\n\n')
        if all_efm_intersection and len(all_efm_intersection):
            f.write(('Intersection of all EFMs has length %d: %s.\n'
                     % (len(all_efm_intersection), all_efm_intersection.to_string(binary=True))) +
                    '--------------------------------------------------------------\n\n')

        for p_id in sorted(p_id2efm_ids.iterkeys(), key=sorter):
            log_pattern(p_id, f)


def serialize_cliques(model, id2clique, cl_id2efm_ids, id2key, key2cl_ids, output_file):
    """
    Serializes cliques to a file, one clique per line. Cliques are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.

    :param id2clique: dict, {clique id: clique}.

    :param output_file: path to the file where the cliques should be saved """

    def log_clique(cl_id, f):
        clique = id2clique[cl_id]
        cl_length = len(clique)
        cl_string = clique.to_string()
        efm_ids = cl_id2efm_ids[cl_id]
        f.write("--------------------------\n\nReaction group %d of length %d:\n\n\t%s\n\nfound in %d EFMs:\n\n\t%s\n\n"
                % (cl_id, cl_length, cl_string, len(efm_ids), ', '.join(sorted(efm_ids))))

    with open(output_file, 'w+') as f:
        f.write('Found %d coupled reaction groups of %d types (based on inputs and outputs).\n\n'
                % (len(id2clique), len(key2cl_ids)))
        f.write('======================================================================\n\n')

        for key_id in sorted(id2key.iterkeys()):
            r_id2st, p_id2st = id2key[key_id]
            cl_ids = key2cl_ids[(r_id2st, p_id2st)]
            f.write('Type %d, contains %d reaction group%s.\n\n'
                    % (key_id, len(cl_ids), 's' if len(cl_ids) != 1 else ''))
            f.write('Inputs: %s;\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                for (m_id, st) in r_id2st))
            f.write('Outputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                   for (m_id, st) in p_id2st))

            for cl_id in sorted(cl_ids):
                log_clique(cl_id, f)

            f.write('======================================================================\n\n')
