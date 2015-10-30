from collections import defaultdict, Counter
import re
from mod_sbml.annotation.kegg.kegg_annotator import get_pathway2r_ids
from mod_sbml.sbml.sbml_manager import get_subsystem2r_ids

__author__ = 'anna'


def get_pathways(model, pts, name2pw, root_ids):
    pw2rs = defaultdict(set)
    pw2r_ids_map = get_pathway2r_ids(model=model)[0]
    if not pw2r_ids_map:
        pw2r_ids_map = get_subsystem2r_ids(model=model)[0]
    else:
        pw2r_ids_map = {re.sub('^[A-Za-z]{3}', 'map', pw): r_ids for (pw, r_ids) in pw2r_ids_map.iteritems()}
    for pw, r_ids in pw2r_ids_map.iteritems():
            pw = pw.lower().strip()
            # term = None
            term = pts.get_term(pw, check_only_ids=False)
            if not term and name2pw and pw in name2pw:
                term = pts.get_term(name2pw[pw], check_only_ids=False)
            if not term:
                keys = [pw]
            else:
                if not term.get_id() in root_ids:
                    keys = {t.get_name() for t in pts.get_ancestors(term, direct=False) if t.get_id() in root_ids}
                else:
                    keys = [term.get_name()]
            for key in keys:
                pw2rs[key] |= r_ids
    return pw2rs


def r_ids2pws(r_ids, pw2rs):
    return {pw for pw, rs in pw2rs.iteritems() if len(r_ids & rs) >= max(2, min(len(r_ids), len(rs)) * 0.25)}


def count_pathways(efm_ids, efm_id2pws):
    result = Counter()
    for efm_id in efm_ids:
        result.update({pw: 1 for pw in efm_id2pws[efm_id]})
    return result

