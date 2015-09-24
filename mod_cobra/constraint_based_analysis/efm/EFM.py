from collections import Counter
import math
import sys
from mod_cobra.constraint_based_analysis import ZERO_THRESHOLD

__author__ = 'anna'


def get_int_size():
    """
    Calculates the maximal number of bits in an int:
    math.log(sys.maxint) / math.log(2).

    :return: int, the maximal number of bits in an int.
    """
    return math.log(sys.maxint) / math.log(2)


def coeff_to_binary(coeff):
    return 1 if coeff > 0 else -1


class EFM(object):
    def __init__(self, r_id2coeff=None):
        """
        Creates an EFM representation from {r_id: coefficient} map.

        :param r_id2coeff: a EFM represented as a dictionary {r_id: coefficient}.
        """
        self.r_id2coeff = r_id2coeff
        self.hash = hash(tuple(sorted(self.r_id2coeff.iteritems())))

    def translate(self, id2id):
        return EFM(r_id2coeff={id2id[r_id]: coeff for (r_id, coeff) in self.to_r_id2coeff().iteritems()})

    def join(self, efms, zero_threshold=ZERO_THRESHOLD):
        if not efms:
            return self
        r_id2coeff = Counter(self.to_r_id2coeff())
        for efm in efms:
            r_id2coeff.update(efm.to_r_id2coeff())
        return EFM(r_id2coeff={r_id: coeff for (r_id, coeff) in r_id2coeff.iteritems() if abs(coeff) > zero_threshold})

    def __str__(self, binary=False):
        return self.to_string()

    def to_string(self, binary=False, subpattern=None, key=None, show_subpattern=True):
        r_id2coefficient = self.to_r_id2coeff()
        subkeys = set()
        if subpattern:
            subkeys = set(r_id2coefficient.iterkeys()) & set(subpattern.to_r_id2coeff().iterkeys())
        keys = sorted(set(r_id2coefficient.iterkeys()) - subkeys, key=key)
        if binary:
            result = ''
            if subkeys and show_subpattern:
                result += '(%s)\t' % '\t'.join('%s%s' % ('-' if r_id2coefficient[r_id] < 0 else '', r_id)
                                             for r_id in sorted(subkeys, key=key))
            return result + '\t'.join('%s%s' % ('-' if r_id2coefficient[r_id] < 0 else '', r_id) for r_id in keys)
        result = ''
        if subkeys and show_subpattern:
            result += '(%s)\t' % '\t'.join('%g %s' % (r_id2coefficient[r_id], r_id) for r_id in sorted(subkeys, key=key))
        return result + '\t'.join('%g %s' % (r_id2coefficient[r_id], r_id) for r_id in keys)

    def to_r_id2coeff(self, binary=False):
        """
        Returns a representation of this EFM as a dictionary
        that maps ids of active reactions to their coefficients: {r_id: coefficient}.
        If binary is set to True, the coefficient values are 1 (for any active reaction in the standard direction)
        or -1 (for reactions that are active in the reversed direction).

        :param binary: boolean, if is set to True, the coefficient values in the result will be
        1 (for any active reaction in the standard direction)
        or -1 (for reactions that are active in the reversed direction).
        Otherwise, any float coefficient values will be possible.

        :return: dict, {r_id: coefficient}.
        """
        if binary:
            return {r_id: coeff_to_binary(coeff) for (r_id, coeff) in self.r_id2coeff.iteritems()}
        return self.r_id2coeff

    def get_sample_r_id_coefficient_pair(self, binary=False):
        r_id, coeff = next(self.r_id2coeff.iteritems())
        return r_id, coeff_to_binary(coeff) if binary else coeff

    def intersection(self, other):
        if not other or not isinstance(other, EFM):
            raise AttributeError('Other should be of type EFM')
        r_id2coeff = \
            dict(set(self.to_r_id2coeff(binary=True).iteritems()) & set(other.to_r_id2coeff(binary=True).iteritems()))
        return EFM(r_id2coeff=r_id2coeff)

    def fold_cliques(self, id2clique, cl_id2new_r_id):
        new_r_id2coeff = Counter()
        replaced_r_ids = set()
        for cl_id, clique in id2clique.iteritems():
            r_id, coeff = clique.get_sample_r_id_coefficient_pair()
            if r_id in self.r_id2coeff:
                ratio = 1.0 * self.r_id2coeff[r_id] / coeff
                if ratio > 0:
                    replaced_r_ids |= set(clique.r_id2coeff.iterkeys())
                    new_r_id2coeff.update({cl_id2new_r_id[cl_id]: ratio})
        if not replaced_r_ids:
            return self
        new_r_id2coeff.update({r_id: coeff for (r_id, coeff) in self.r_id2coeff.iteritems()
                               if r_id not in replaced_r_ids})
        return EFM(r_id2coeff=new_r_id2coeff)

    def intersection_len(self, other):
        return len(set(self.to_r_id2coeff(binary=True).iteritems()) & set(other.to_r_id2coeff(binary=True).iteritems()))

    def __len__(self):
        """
        Returns the length of this EFM, i.e. the number of active reactions.
        """
        return len(self.r_id2coeff)

    def __eq__(self, other):
        if not other or not isinstance(other, EFM):
            return False
        return self.hash == other.hash

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self.hash




