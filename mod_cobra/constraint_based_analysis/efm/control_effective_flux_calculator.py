__author__ = 'anna'


def get_efm_efficiency(efm, r_id):
    """
    We follow the definition of EFM efficiency by Stelling et al. Nature. 2002;420(6912):190-3.
    that relates the system's output to the sum of all absolute fluxes:
    For an EFM $e_i$ its efficiency for a reaction $\mu$ is defined as
    $$\varepsilon_i(\mu)=\frac{e^{\mu}_i}{\sum_{l}|e^l_i|}
    \mbox{, where } e^{j}_{i} \mbox{ is a rate of the reaction } j.$$

    :param efm: EFM ($e_i$)
    :param r_id: id of the reaction of interest (reaction $\mu$)
    :return: int, efficiency
    """
    r_id2coeff = efm.to_r_id2coeff()
    if not r_id2coeff or r_id not in r_id2coeff:
        return 0
    return r_id2coeff[r_id] / sum(abs(coeff) for coeff in r_id2coeff.itervalues())


def classify_efm_by_efficiency(id2efm, r_id):
    return {efm_id: get_efm_efficiency(efm, r_id) for (efm_id, efm) in id2efm.iteritems()}