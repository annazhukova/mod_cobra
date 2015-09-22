__author__ = 'anna'


def get_fm_efficiency(fm, r_id, r_rev=False):
    """
    We follow the definition of EFM efficiency by Stelling et al. Nature. 2002;420(6912):190-3.
    that relates the system's output to the sum of all absolute fluxes:
    For an EFM $e_i$ its efficiency for a reaction $\mu$ is defined as
    $$\varepsilon_i(\mu)=\frac{e^{\mu}_i}{\sum_{l}|e^l_i|}
    \mbox{, where } e^{j}_{i} \mbox{ is a rate of the reaction } j.$$

    :param fm: EFM ($e_i$)
    :param r_id: id of the reaction of interest (reaction $\mu$)
    :param r_rev: wether the reaction of interest (reaction $\mu$) should be considered in reversed direction
    :return: int, efficiency
    """
    r_id2coeff = fm.to_r_id2coeff()
    if not r_id2coeff or r_id not in r_id2coeff:
        return 0
    return (-1.0 if r_rev else 1.0) * r_id2coeff[r_id] / sum(abs(coeff) for coeff in r_id2coeff.itervalues())


def get_fm_yield(fm, r_id, in_r_id, r_rev=False, in_r_rev=False):
    r_id2coeff = fm.to_r_id2coeff()
    if not r_id2coeff or r_id not in r_id2coeff or not in_r_id in r_id2coeff:
        return 0
    return (-1.0 if r_rev else 1.0) * r_id2coeff[r_id] / ((-1.0 if in_r_rev else 1.0) * r_id2coeff[in_r_id])


def classify_efm_by_efficiency(id2efm, r_id, r_rev=False):
    return {efm_id: get_fm_efficiency(efm, r_id, r_rev) for (efm_id, efm) in id2efm.iteritems()}