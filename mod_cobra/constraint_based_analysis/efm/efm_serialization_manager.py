import libsbml

from mod_cobra.constraint_based_analysis.efm.control_effective_flux_calculator import get_fm_yield
from mod_sbml.utils.misc import invert_map
from mod_cobra.constraint_based_analysis.efm.EFM import EFM
from mod_sbml.sbml.submodel_manager import submodel
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula

THICK_DELIMITER = '==============================================================\n\n'
THIN_DELIMITER = '------------------------------------------\n\n'

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'


def serialize_fms_txt(id2fm, id2efm, fm_id2key, key2folded_efm_ids, path, efm_id2efficiency,
                      id2folded_efm, folded_efm_id2efm_ids, new_r_id2cl_id, id2clique,
                      all_efm_intersection, model, r_id, r_rev, in_r_id, in_r_rev):
    with open(path, 'w+') as f:
        f.write('Found %d EFMs, folded into %d EFMs, grouped into %d pathways:\n\n'
                % (len(id2efm), len(id2folded_efm), len(id2fm)))
        f.write(THICK_DELIMITER)
        if all_efm_intersection and len(all_efm_intersection):
            f.write('All EFMs contain following %d reactions:\n\t%s.\n\n'
                    % (len(all_efm_intersection), all_efm_intersection.to_string(binary=True)))
            f.write(THICK_DELIMITER)
        for fm_id in sorted(id2fm.iterkeys()):
            fm = id2fm[fm_id]
            r_id2st, p_id2st = fm_id2key[fm_id]
            folded_efm_ids = key2folded_efm_ids[r_id2st, p_id2st]
            yield_str = '' if not in_r_id else ' of yield %g' % get_fm_yield(fm, r_id, in_r_id, r_rev, in_r_rev)
            if len(folded_efm_ids) == 1:
                folded_efm_id = next(iter(folded_efm_ids))
                folded_efm = id2efm[folded_efm_id]

                f.write('Folded EFM %s (a.k.a. pathway %s) of length %d%s.\n\n'
                        % (folded_efm_id, fm_id, len(folded_efm), yield_str))
                caption = 'Structure:\n\n\t'
            else:
                f.write('Pathway %s of length %d%s.\n\n' % (fm_id, len(fm), yield_str))
                caption = None

            f.write('Inputs: %s;\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                for (m_id, st) in r_id2st))
            f.write('Outputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                   for (m_id, st) in p_id2st))

            if len(folded_efm_ids) > 1:
                f.write('Structure:\n\n\t%s\n\n' % fm.to_string(subpattern=all_efm_intersection))
                f.write('Contains %d folded FMs:\n\n' % len(folded_efm_ids))

            started = False
            for folded_efm_id in folded_efm_ids:
                if started:
                    f.write(THIN_DELIMITER)
                else:
                    started = True
                folded_efm = id2efm[folded_efm_id]
                f.write(caption if caption else 'Folded EFM %s of length %d:\n\n\t' % (folded_efm_id, len(folded_efm)))
                f.write('%s\n\n' % folded_efm.to_string(subpattern=all_efm_intersection))

                f.write('Unfolds into %d EFM%s:\n\n' % (len(folded_efm_id2efm_ids[folded_efm_id]),
                                                        's' if len(folded_efm_id2efm_ids[folded_efm_id]) != 1 else ''))
                for efm_id in folded_efm_id2efm_ids[folded_efm_id]:
                    efm = id2efm[efm_id]
                    f.write('EFM %s of length %d of efficiency %g:\n\n\t%s\n\n'
                            % (efm_id, len(efm), efm_id2efficiency[efm_id],
                               efm.to_string(subpattern=all_efm_intersection)))
                    f.write(THIN_DELIMITER)
            f.write(THICK_DELIMITER)


def serialize_communities(id2cluster, id2intersection, id2imp_rns, total_len, all_fm_intersection, path):
    with open(path, 'w+') as f:
        f.write('Analysed %d pathways\n\n' % total_len)
        f.write(THIN_DELIMITER)
        if all_fm_intersection and len(all_fm_intersection):
            f.write('All pathways contain following %d reactions:\n\n\t%s.\n\n'
                    % (len(all_fm_intersection), all_fm_intersection.to_string(binary=True)))
            f.write(THICK_DELIMITER)
        f.write('Found %d communities\n\n' % len(id2cluster))
        all_r_ids = set(all_fm_intersection.to_r_id2coeff().iterkeys())
        for clu_id in sorted(id2cluster.iterkeys()):
            f.write(THIN_DELIMITER)
            cluster = id2cluster[clu_id]
            intersection = id2intersection[clu_id]
            imp_rns = id2imp_rns[clu_id]
            f.write('Community %d contains %d following pathways:\n\n\t%s\n\n'
                    % (clu_id, len(cluster), ', '.join(sorted(cluster))))
            f.write('all of which contain following reactions:\n\n\t%s\n\n'
                    % (intersection.to_string(binary=True, subpattern=all_fm_intersection)))
            f.write('which form a reaction community, which also contains following reactions:\n\n\t %s\n\n'
                    % imp_rns.to_string(binary=True, subpattern=intersection, show_subpattern=False,
                                        key=lambda r_id: (0 if r_id in all_r_ids else 1, r_id)))


def read_efms(output_efm_file):
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
            id2efm[efm_id] = EFM(r_id2coeff={r_id: float(coeff) for (coeff, r_id) in (it.split(' ') for it in efm)})
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


def serialize_clique(cl_id, clique, efm_ids, model, output_file):
    """
    Serializes a clique to a file.

    :param cl_id: clique id

    :param clique: EFM.

    :param output_file: path to the file where the clique should be saved
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = clique.to_r_id2coeff(binary=True)
        f.write('Reaction group %d of length %d found in %d EFMs:\n\t%s\n\n'
                % (cl_id, len(r_id2coeff), len(efm_ids), ', '.join(sorted(efm_ids))))
        f.write(THIN_DELIMITER)
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n'
                        % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
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
        f.write('Community %d of %d following EFMs:\n\n\t%s\n\n' % (c_id, len(efm_ids), ', '.join(sorted(efm_ids))))
        f.write(THIN_DELIMITER)
        f.write('Reactions that participate in all EFMs of this community (%d):\n\n' % len(r_id2coeff))
        f.write(THIN_DELIMITER)
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id),
                                                                           show_compartments=False,
                                                                           show_metabolite_ids=True)))
            f.write('\n')
        if imp_r_id2coeff:
            f.write(THIN_DELIMITER)
            f.write('Other reactions that participate in the same reaction community (%d):\n\n' % len(imp_r_id2coeff))
            f.write(THIN_DELIMITER)
            coeff2r_id = invert_map(imp_r_id2coeff)
            for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
                for r_id in sorted(r_ids):
                    f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id),
                                                                               show_compartments=False,
                                                                               show_metabolite_ids=True)))
                f.write('\n')


def serialize_fm(fm_id, fm, model, output_file, r_id, in_r_id, r_rev, in_r_rev):
    """
    Serializes a FM to a file.

    :param fm_id: FM id

    :param fm: FM.

    :param output_file: path to the file where the FM should be saved
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = fm.to_r_id2coeff()
        yield_str = '' if not in_r_id else ' of yield %g' % get_fm_yield(fm, r_id, in_r_id, r_rev, in_r_rev)
        f.write('Pathway %s of length %d%s\n\n' % (fm_id, len(r_id2coeff), yield_str))
        f.write(THIN_DELIMITER)
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' %
                        (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
                                                         show_metabolite_ids=True)))
            f.write('\n')


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
        f.write(THIN_DELIMITER)
        f.write("Reaction group %d of length %d:\n\n\t%s\n\nfound in %d EFMs:\n\n\t%s\n\n"
                % (cl_id, cl_length, cl_string, len(efm_ids), ', '.join(sorted(efm_ids))))

    with open(output_file, 'w+') as f:
        f.write('Found %d coupled reaction groups of %d types (based on inputs and outputs).\n\n'
                % (len(id2clique), len(key2cl_ids)))
        f.write(THICK_DELIMITER)

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

            f.write(THICK_DELIMITER)
