from collections import defaultdict
import logging

__author__ = 'anna'


def classify_efms(id2efm, min_pattern_len, min_efm_num=2, max_pattern_num=None):
    """
    Classifies EFMs to find common patterns.

    :param id2efm: dict {efm_id: efm}, EFMs are instances of class EFM.

    :param min_pattern_len: int, minimal length for a pattern to be considered.

    :param min_efm_num: int, (optional) at least how many EFMs should include a pattern for it to be considered.
    (The default value is 2.)

    :param max_pattern_num: int, (optional) at most how many patterns should be returned.
    If not set, all the patterns will be returned.

    :return: 2 dictionaries: p_id2efm_ids and id2pattern.
    p_id2efm_ids maps a pattern_id to ids of the EFMs containing this pattern,
    id2pattern maps a pattern_id to the pattern (represented as an instance of class EFM).
    """

    pattern2efm_ids = defaultdict(set)
    pattern2efm_ids.update({efm: {efm_id} for (efm_id, efm) in id2efm.iteritems()})

    def process_pattern_intersection(pattern1, pattern2, detected_patterns):
        intersection = pattern1.intersection(pattern2)
        if len(intersection) >= min_pattern_len:
            if intersection not in pattern2efm_ids:
                detected_patterns.add(intersection)
            pattern2efm_ids[intersection] |= pattern2efm_ids[pattern1] | pattern2efm_ids[pattern2]
            return intersection
        return None

    efm_ids = set(id2efm.iterkeys())
    patterns = set(pattern2efm_ids.iterkeys())

    level = 2
    while patterns:
        logging.info("Calculating patterns contained in at least %d EFMs..." % level)
        new_patterns = set()
        patterns_to_remove = set()
        for pattern in patterns:
            p_efm_ids = pattern2efm_ids[pattern]
            for efm_id in efm_ids - p_efm_ids:
                efm = id2efm[efm_id]
                new_pattern = process_pattern_intersection(pattern, efm, new_patterns)
                # if we found a subpattern that is long enough and is common for more EFMs
                if new_pattern and new_pattern != pattern:
                    patterns_to_remove.add(pattern)
        patterns = new_patterns
        for p in patterns_to_remove:
            del pattern2efm_ids[p]
        logging.info("... found %d patterns." % len(patterns))
        level += 1
    max_pattern_num = min(len(pattern2efm_ids), max_pattern_num if max_pattern_num else len(pattern2efm_ids))
    patterns = sorted((p for p in pattern2efm_ids.iterkeys() if len(pattern2efm_ids[p]) >= min_efm_num),
                      key=lambda p: -len(pattern2efm_ids[p]))[0: max_pattern_num]
    id2pattern = dict(zip(xrange(1, len(patterns) + 1), patterns))
    p_id2efm_ids = {p_id: pattern2efm_ids[p] for (p_id, p) in id2pattern.iteritems()}
    return p_id2efm_ids, id2pattern




