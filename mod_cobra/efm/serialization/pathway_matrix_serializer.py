import io
import os
from mod_cobra.html import describe
from mod_cobra import round_value
from mod_sbml.serialization import format_m_name

__author__ = 'anna'


def serialize(S, model, in_m_id, out_m_id, model_name, main_dir):
    def format_r(r_id):
        r2st, p2st = S.st_matrix.get_inputs_outputs(r_id)
        formatter = lambda (m_id, st): '%s %s' % ("" if st == 1 else ("%g" % st),
                                                  format_m_name(model.getSpecies(m_id), model, show_id=False))
        rs, ps = sorted(formatter(it) for it in r2st.iteritems()), sorted(formatter(it) for it in p2st.iteritems())
        r = model.getReaction(r_id)
        reversible = r.getReversible() if r else True

        return " + ".join(rs) + (" &#8596; " if reversible else " &#8594; ") + " + ".join(ps)

    def format_p(p_id):
        fm_yield = S.get_yield(p_id, in_m_id, out_m_id)
        if fm_yield is not None:
            return 'Yield: %g' % round_value(fm_yield)
        return ''

    r_id2tooltip = {r_id: format_r(r_id) for r_id in S.r_id2i.iterkeys()}
    p_id2tooltip = {p_id: format_p(p_id) for p_id in S.efm_id2i.iterkeys()}

    efm_id2tr = {}

    def to_tr(efm_id, header_class='header', color_class='pw', super_class='', hidden=False):
        if efm_id in efm_id2tr:
            return
        if efm_id in S.gr_id2efm_ids:
            for sub_efm_id in S.gr_id2efm_ids[efm_id]:
                to_tr(sub_efm_id, header_class='sub_%s' % header_class, super_class='sub_%s' % efm_id,
                      color_class='sub_%s' % color_class, hidden=True)
        efm_id2tr[efm_id] = describe('foldable_rows.html', S=S, efm_id=efm_id, header_class=header_class,
                                     super_class=super_class, efm_id2tr=efm_id2tr, color_class=color_class,
                                     p_id2tooltip=p_id2tooltip, hidden=hidden)

    for efm_id in S.efm_ids:
        to_tr(efm_id)

    content = describe('V.html', S=S, efm_id2tr=efm_id2tr, r_id2tooltip=r_id2tooltip, model_name=model_name)
    html_file_name = '%s_pathways.html' % model_name
    with io.open(os.path.join(main_dir, html_file_name), 'w+', encoding='utf-8') as f:
        f.write(content)
    return html_file_name
