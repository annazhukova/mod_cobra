import libsbml

from mod_cobra.constraint_based_analysis.efm.control_effective_flux_calculator import get_fm_yield
from mod_sbml.utils.misc import invert_map
from mod_cobra.constraint_based_analysis.efm.EFM import EFM, TYPE_EFM, TYPE_FOLDED_EFM
from mod_sbml.sbml.submodel_manager import submodel
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula

THICK_DELIMITER = '==============================================================\n\n'
THIN_DELIMITER = '------------------------------------------\n\n'
TINY_DELIMITER = '-------------\n\n'

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'


def serialize_fms_txt(id2fm, id2efm, fm_id2key, key2folded_efm_ids, path, efm_id2efficiency,
                      id2folded_efm, folded_efm_id2efm_ids, r_id2new_r_id, r_id2cl_id,
                      all_fm_intersection, all_fm_intersection_folded, model, in_m_id, out_m_id):
    with open(path, 'w+') as f:
        f.write('Found %d EFMs, folded into %d EFMs, grouped into %d pathways:\n\n'
                % (len(id2efm), len(id2folded_efm), len(id2fm)))
        f.write(THICK_DELIMITER)

        def get_key(r_id, c):
            key = [0] if r_id in all_fm_intersection or r_id in all_fm_intersection_folded else []
            if (r_id, c) in r_id2new_r_id:
                key.append(r_id2new_r_id[(r_id, c)])
            if (r_id, c) in r_id2cl_id:
                key.append(r_id2cl_id[(r_id, c)])
            key.append(r_id)
            return tuple(key)

        get_group = lambda r_id, c: r_id2cl_id[(r_id, c)] if (r_id, c) in r_id2cl_id else None

        if all_fm_intersection and len(all_fm_intersection):
            f.write('All EFMs contain following %d reactions:\n\t%s,\n\n'
                    % (len(all_fm_intersection),
                       all_fm_intersection.to_string(binary=True, get_key=get_key, get_group=get_group)))
            if len(all_fm_intersection) > len(all_fm_intersection_folded):
                f.write('or (with reaction groups folded):\n\t%s.\n\n'
                        % all_fm_intersection_folded.to_string(binary=True, get_key=get_key))
            f.write(THICK_DELIMITER)

        fm_id2yield = {fm_id: get_fm_yield(fm_id2key[fm_id], in_m_id, out_m_id) for fm_id in id2fm.iterkeys()}
        for fm_id, fm in sorted(id2fm.iteritems(),
                                key=lambda (fm_id, fm): (-fm_id2yield[fm_id] if fm_id2yield[fm_id] else 0, fm_id)):
            write_pathway(fm_id, fm, fm_id2key, model, efm_id2efficiency,
                          folded_efm_id2efm_ids, id2efm, key2folded_efm_ids, id2folded_efm, fm_id2yield[fm_id],
                          get_key, get_group, f)
            f.write(THICK_DELIMITER)


def write_pathway(fm_id, fm, fm_id2key, model, efm_id2efficiency, folded_efm_id2efm_ids,
                  id2efm, key2folded_efm_ids, id2folded_efm, fm_yield, get_key, get_group, f):
    r_id2st, p_id2st = fm_id2key[fm_id]
    f.write('Inputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                        for (m_id, st) in r_id2st))
    f.write('Outputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                           for (m_id, st) in p_id2st))
    if fm_yield is not None:
        f.write('Yield: %g;\n\n' % fm_yield)

    if TYPE_FOLDED_EFM == fm.type:
        write_folded_efm(fm_id, fm, efm_id2efficiency, folded_efm_id2efm_ids, id2efm, get_key, get_group, f)
        return

    if TYPE_EFM == fm.type:
        write_efm(fm_id, fm, efm_id2efficiency, get_key, get_group, f)
        return

    f.write('%s of length %d:\n\n\t%s\n\n' % (fm_id, len(fm), fm.to_string(get_key=get_key, get_group=get_group)))
    folded_efm_ids = key2folded_efm_ids[r_id2st, p_id2st]
    f.write('\tContains %d elements:\n\n' % len(folded_efm_ids))
    started = False
    for f_efm_id in folded_efm_ids:
        if started:
            f.write('\t\t%s' % THIN_DELIMITER)
        else:
            started = True
        write_folded_efm(f_efm_id, id2folded_efm[f_efm_id], efm_id2efficiency, folded_efm_id2efm_ids, id2efm,
                         get_key, get_group, f, tab='\t\t')


def write_folded_efm(f_efm_id, f_efm, efm_id2efficiency, folded_efm_id2efm_ids, id2efm, get_key, get_group, f, tab=''):
    if TYPE_EFM == f_efm.type:
        write_efm(f_efm_id, f_efm, efm_id2efficiency, get_key, get_group, f, tab)
        return

    f.write('%s%s of length %d:\n\n%s\t%s\n\n'
            % (tab, f_efm_id, len(f_efm), tab, f_efm.to_string(get_key=get_key, get_group=get_group)))

    efm_ids = folded_efm_id2efm_ids[f_efm_id]
    f.write('%s\tUnfolds into %s' % (tab, ('the following %d EFMs:\n\n' % len(efm_ids)) if len(efm_ids) != 1 else ''))

    started = False
    for efm_id in efm_ids:
        if started:
            f.write('%s\t\t%s' % (tab, TINY_DELIMITER))
        else:
            started = True
        write_efm(efm_id, id2efm[efm_id], efm_id2efficiency, get_key, get_group, f,
                  tab='%s\t\t' % tab, no_first_tab=len(efm_ids) == 1)


def write_efm(efm_id, efm, efm_id2efficiency, get_key, get_group, f,  tab='', no_first_tab=False):
    f.write('%s%s of length %d of efficiency %g:\n\n%s\t%s\n\n'
            % ('' if no_first_tab else tab, efm_id, len(efm), efm_id2efficiency[efm_id], tab,
               efm.to_string(get_key=get_key, get_group=get_group)))


def serialize_communities(id2cluster, id2intersection, id2imp_rns, total_len, all_fm_intersection, path):
    with open(path, 'w+') as f:
        f.write('Analysed %d pathways. ' % total_len)
        if all_fm_intersection and len(all_fm_intersection):
            f.write('All pathways contain following %d reactions:\n\n\t%s.\n\n'
                    % (len(all_fm_intersection), all_fm_intersection.to_string(binary=True)))
        else:
            f.write('\n\n')
        f.write(THICK_DELIMITER)

        get_key = lambda r_id, c: (0, r_id) if all_fm_intersection and r_id in all_fm_intersection else (1, r_id)
        get_group = lambda r_id, c: 1 if r_id in all_fm_intersection else None

        f.write('Found %d communities\n\n' % len(id2cluster))
        for clu_id in sorted(id2cluster.iterkeys()):
            f.write(THIN_DELIMITER)
            cluster = id2cluster[clu_id]
            intersection = id2intersection[clu_id]
            hidden_r_ids = set(intersection.to_r_id2coeff().iterkeys())
            imp_rns = id2imp_rns[clu_id]
            f.write('Community %d contains %d following pathways:\n\n\t%s\n\n'
                    % (clu_id, len(cluster), ', '.join(sorted(cluster))))
            f.write('all of which contain following reactions:\n\n\t%s\n\n'
                    % (intersection.to_string(binary=True, get_key=get_key, get_group=get_group)))
            f.write('which form a reaction community, which also contains following reactions:\n\n\t %s\n\n'
                    % imp_rns.to_string(binary=True, hidden_r_ids=hidden_r_ids))


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
            id2efm[efm_id] = EFM(r_id2coeff={r_id: float(coeff) for (coeff, r_id) in (it.split(' ') for it in efm)},
                                 type=TYPE_EFM)
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
        f.write('Reaction group %d of length %d found in %d EFMs:\n\n\t%s\n\n'
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


def serialize_fm(fm_id, fm, (r_id2st, p_id2st), in_m_id, out_m_id, model, output_file):
    """
    Serializes a FM to a file.

    :param fm_id: FM id

    :param fm: FM.

    :param output_file: path to the file where the FM should be saved
    """
    with open(output_file, 'w+') as f:
        r_id2coeff = fm.to_r_id2coeff()
        yield_str = '' if not in_m_id or not out_m_id else \
            ' of yield %g' % get_fm_yield((r_id2st, p_id2st), in_m_id, out_m_id)
        f.write('%s of length %d%s\n\n' % (fm_id, len(r_id2coeff), yield_str))
        f.write(THIN_DELIMITER)
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' %
                        (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), show_compartments=False,
                                                         show_metabolite_ids=True)))
            f.write('\n')


def serialize_cliques(model, id2clique, cl_id2efm_ids, id2key, key2cl_ids, key2r_ids, output_file):
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
        f.write('\t%s' % THIN_DELIMITER)
        f.write("\tReaction group %d of length %d:\n\n\t\t%s\n\n\tfound in %d EFMs:\n\n\t\t%s\n\n"
                % (cl_id, cl_length, cl_string, len(efm_ids), ', '.join(sorted(efm_ids))))

    with open(output_file, 'w+') as f:
        f.write('Found %d coupled reaction groups of %d types (Types are based on inputs and outputs).\n\n'
                % (len(id2clique), len(key2cl_ids)))
        f.write(THICK_DELIMITER)

        for key_id in sorted(id2key.iterkeys()):
            r_id2st, p_id2st = id2key[key_id]
            cl_ids = key2cl_ids[(r_id2st, p_id2st)]
            f.write('Type %d, contains %d reaction group%s.\n\n' % (key_id, len(cl_ids), 's' if len(cl_ids) != 1 else ''))
            f.write('Inputs: %s;\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                for (m_id, st) in r_id2st))
            f.write('Outputs: %s;\n\n' % ', '.join('%g %s (%s)' % (st, model.getSpecies(m_id).getName(), m_id)
                                                   for (m_id, st) in p_id2st))

            for cl_id in sorted(cl_ids):
                log_clique(cl_id, f)

            if (r_id2st, p_id2st) in key2r_ids:
                f.write('\t%s' % THIN_DELIMITER)
                r_ids = key2r_ids[(r_id2st, p_id2st)]
                f.write('\tReaction%s %s also ha%s the same structure.\n\n'
                        % ('s' if len(r_ids) != 1 else '',
                           ', '.join((('-%s' % r_id) if c < 0 else r_id for (r_id, c) in r_ids)),
                           's' if len(r_ids) == 1 else 've'))

            f.write(THICK_DELIMITER)
