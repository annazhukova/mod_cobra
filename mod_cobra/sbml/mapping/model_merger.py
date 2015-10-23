from collections import defaultdict
from numpy import hstack, vstack, zeros

import libsbml
from scipy.linalg import block_diag

from mod_cobra.efm.System import System
from mod_sbml.sbml.sbml_manager import create_reaction, create_species, create_compartment

__author__ = 'anna'


def simple_merge_models(S, model_id2c_id2group, model_id2dfs, out_sbml):
    doc = libsbml.SBMLDocument(2, 4)
    model = doc.createModel()

    model.setId('merged_model')

    model_id2id2id = defaultdict(dict)
    common_ids = set()
    c_group2id = {}

    new_m_id2i, new_r_id2i, new_efm_id2i, new_boundary_m_ids = {}, {}, {}, []

    for model_id, [_, _, df] in model_id2dfs.iteritems():
        for index, row in df.iterrows():
            c_id, name = row['Id'], row['Name']
            if model_id in model_id2c_id2group and c_id in model_id2c_id2group[model_id]:
                group = model_id2c_id2group[model_id][c_id]
                if group in c_group2id:
                    new_id = c_group2id[group]
                else:
                    new_id = create_compartment(model, name=name, id_='merged_%s_%s' % (model_id, c_id)).getId()
                    c_group2id[group] = new_id
                    common_ids.add(new_id)
            else:
                new_id = create_compartment(model, name=name, id_='%s_%s' % (model_id, c_id)).getId()
            model_id2id2id[model_id][c_id] = new_id

    id2id = {}

    m_id_group_ids = set(S.m_id2gr_id.itervalues())
    for (model_id, m_id), i in ((it, i) for (it, i) in S.m_id2i.iteritems() if it not in m_id_group_ids):
        c_id = model_id2dfs[model_id][0].at[m_id, 'Compartment']
        c_id = model_id2id2id[model_id][c_id]
        name = model_id2dfs[model_id][0].at[m_id, 'Name']
        is_boundary = (model_id, m_id) in S.boundary_m_ids
        new_id = create_species(model, compartment_id=c_id, name=name, bound=is_boundary,
                                id_='%s_%s' % (model_id, m_id)).getId()
        model_id2id2id[model_id][m_id] = new_id
        id2id[(model_id, m_id)] = new_id
        new_m_id2i[new_id] = i
        if is_boundary:
            new_boundary_m_ids.append(new_id)

    for it in m_id_group_ids:
        model_id, m_ids = next(iter(it))
        m_id = next(iter(m_ids))
        is_boundary = it in S.boundary_m_ids
        new_id = \
            create_species(model,
                           compartment_id=model_id2id2id[model_id][model_id2dfs[model_id][0].at[m_id, 'Compartment']],
                           name=model_id2dfs[model_id][0].at[m_id, 'Name'], bound=is_boundary,
                           id_='merged_%s_%s' % (model_id, m_id)).getId()
        for model_id, m_ids in it:
            model_id2id2id[model_id].update({m_id: new_id for m_id in m_ids})
        id2id[it] = new_id
        new_m_id2i[new_id] = S.m_id2i[it]
        if is_boundary:
            new_boundary_m_ids.append(new_id)
        common_ids.add(new_id)

    for ((model_id, r_id), i) in ((it, i) for (it, i) in S.r_id2i.iteritems() if it not in S.gr_id2r_id2c.keys()):
        r_id2st, p_id2st = S.st_matrix.get_inputs_outputs((model_id, r_id))
        new_id = create_reaction(model, {id2id[m_id]: st for (m_id, st) in r_id2st.iteritems()},
                                 {id2id[m_id]: st for (m_id, st) in p_id2st.iteritems()},
                                 model_id2dfs[model_id][1].at[r_id, 'Name'], reversible=True,
                                 id_='%s_%s' % (model_id, r_id)).getId()
        model_id2id2id[model_id][r_id] = new_id
        new_r_id2i[new_id] = i

    for gr, it2c in S.gr_id2r_id2c.iteritems():
        model_id, r_id = next(it2c.iterkeys())
        r_id2st, p_id2st = S.st_matrix.get_inputs_outputs(gr)
        new_id = \
            create_reaction(model, {id2id[m_id]: st for (m_id, st) in r_id2st.iteritems()},
                            {id2id[m_id]: st for (m_id, st) in p_id2st.iteritems()},
                            model_id2dfs[model_id][1].at[r_id, 'Name'], reversible=True,
                            id_='merged_%s_%s' % (model_id, r_id)).getId()
        for model_id, r_id in it2c.iterkeys():
            model_id2id2id[model_id][r_id] = new_id
        new_r_id2i[new_id] = S.r_id2i[gr]
        common_ids.add(new_id)

    for ((model_id, efm_id), i) in ((it, i) for (it, i) in S.efm_id2i.iteritems() if it not in S.gr_id2efm_ids.keys()):
        new_id = '%s_%s' % (model_id, efm_id)
        new_efm_id2i[new_id] = i

    for gr, efm_ids in S.gr_id2efm_ids.iteritems():
        model_id, efm_id = next(efm_ids)
        new_id = 'merged_%s_%s' % (model_id, efm_id)
        new_efm_id2i[new_id] = S.r_id2i[gr]

    libsbml.writeSBMLToFile(doc, out_sbml)

    return model_id2id2id, common_ids, System(m_id2i=new_m_id2i, r_id2i=new_r_id2i, efm_id2i=new_efm_id2i, N=S.N, V=S.V,
                                              boundary_m_ids=new_boundary_m_ids)


def join(model_id2m_ids_groups, model_id2S):
    common_m_ids = set()
    for model_id2m_ids in model_id2m_ids_groups:
        for model_id, m_ids in model_id2m_ids.iteritems():
            common_m_ids |= {(model_id, m_id) for m_id in m_ids}

    N = block_diag(*[S.N[sorted(S.m_id2i[m_id] for m_id in S.m_ids if (model_id, m_id) not in common_m_ids), :]
                     [:, sorted(S.r_id2i[r_id] for r_id in S.r_ids)]
                     for (model_id, S) in model_id2S.iteritems()])
    m_id2i, r_id2i, efm_id2i = {}, {}, {}
    r_shift, m_shift, efm_shift = 0, 0, 0
    boundary_m_ids = []
    for model_id, S in model_id2S.iteritems():
        r_id2i.update({(model_id, r_id): i + r_shift
                       for (i, r_id) in enumerate(sorted(S.r_ids, key=lambda r_id: S.r_id2i[r_id]))})
        r_shift += len(S.r_ids)
        efm_id2i.update({(model_id, efm_id): i + efm_shift
                         for (i, efm_id) in enumerate(sorted(S.efm_ids, key=lambda efm_id: S.efm_id2i[efm_id]))})
        efm_shift += len(S.efm_ids)
        m_ids = {m_id for m_id in S.m_ids if (model_id, m_id) not in common_m_ids}
        m_id2i.update({(model_id, m_id): i + m_shift
                       for (i, m_id) in enumerate(sorted(m_ids, key=lambda m_id: S.m_id2i[m_id]))})
        boundary_m_ids += [(model_id, m_id) for m_id in set(S.boundary_m_ids) & m_ids]
        m_shift += len(m_ids)
    m_id2new_m_id = dict()
    for model_id2m_ids in model_id2m_ids_groups:
        model_id2m_row = {model_id: zeros(shape=(1, len(S.r_ids))) for model_id, S in model_id2S.iteritems()}
        new_m_id = tuple(((model_id, tuple(m_ids)) for (model_id, m_ids) in model_id2m_ids.iteritems()))
        m_id2i[new_m_id] = m_shift
        if next((True for (model_id, m_ids) in model_id2m_ids.iteritems()
                 if m_ids & set(model_id2S[model_id].boundary_m_ids)), False):
            boundary_m_ids.append(new_m_id)
        for model_id, m_ids in model_id2m_ids.iteritems():
            S = model_id2S[model_id]
            r_is = sorted(S.r_id2i[r_id] for r_id in S.r_ids)
            for m_id in m_ids:
                m_id2new_m_id[(model_id, m_id)] = new_m_id
                model_id2m_row[model_id] += S.N[S.m_id2i[m_id], r_is]
        row = hstack(tuple(model_id2m_row[model_id] for model_id in model_id2S.iterkeys()))
        N = vstack((N, row))
        m_shift += 1
    V = block_diag(*[S.V[sorted(S.r_id2i[r_id] for r_id in S.r_ids), :]
                     [:, sorted(S.efm_id2i[efm_id] for efm_id in S.efm_ids)]
                     for (model_id, S) in model_id2S.iteritems()])
    return System(N=N, m_id2i=m_id2i, r_id2i=r_id2i, V=V, efm_id2i=efm_id2i, m_id2gr_id=m_id2new_m_id,
                  boundary_m_ids=boundary_m_ids)


def merge(S, ignored_ids=None):
    S.remove_reaction_duplicates(ignored_ids)
    S.m_ids -= ignored_ids
