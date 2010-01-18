"""Sample estimated branch lengths using three different amounts of information.

The first method observes all of the changes along the branch (infinite sites model).
For the second method multiple changes are observed as a single change (infinite alleles model).
For the third method multiple changes are observed as a single change unless they are reverted,
in which case they are observed as no change (Jukes-Cantor model).
"""

import StringIO
import math
import random

import numpy.random

from SnippetUtil import HandlingError
import RUtil
import Form

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('branch_length', 'branch length', 1.0, low_exclusive=0.0),
            Form.Integer('sequence_length', 'sequence length', 1000, low=1)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the branch length and the sequence length
    branch_length = fs.branch_length
    sequence_length = fs.sequence_length
    # sample sequence changes at three levels of informativeness
    sequence_changes = sample_sequence_changes(branch_length, sequence_length)
    # get a distance estimate for each level of informativeness
    first_estimate, second_estimate, third_estimate = sample_distance(*sequence_changes)
    # begin the response
    out = StringIO.StringIO()
    print >> out, 'distance estimates:'
    print >> out, 'using all change information:', first_estimate
    print >> out, 'without multiple change information:', second_estimate
    print >> out, 'without reverted change information:', third_estimate
    # return the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

def sample_distance(mean_changes, p_observed, p_nonreverted):
    """
    Get a distance esimate for each of three data sources.
    @param mean_changes: the
    @return: a triple of distance estimates
    """
    # when each change is observed, the estimate is the mean number of changes
    first_estimate = mean_changes
    # when multiple changes are hidden, the distance can still be estimated
    if p_observed == 1.0:
        second_estimate = float('inf')
    else:
        second_estimate = -math.log(1 - p_observed)
    # when changes can be reverted, the distance can still be estimated
    if p_nonreverted >= 0.75:
        third_estimate = float('inf')
    else:
        third_estimate = - .75 * math.log(1 - (4.0/3.0) * p_nonreverted)
    return first_estimate, second_estimate, third_estimate

def sample_sequence_changes(branch_length, nsites):
    """
    @param branch_length: the expected number of changes along the branch at each site
    @param nsites: the number of sites in the sequence
    @return: a triple of (mean changes, change frequency, non-reverted change frequency)
    """
    samples = [sample_site_changes(branch_length) for i in range(nsites)]
    change_counts, observed_changes, nonreverted_changes = zip(*samples)
    mean_changes = sum(change_counts) / float(nsites)
    p_observed = sum(observed_changes) / float(nsites)
    p_nonreverted = sum(nonreverted_changes) / float(nsites)
    return mean_changes, p_observed, p_nonreverted

def sample_site_changes(branch_length):
    """
    @param branch_length: the expected number of changes along the branch
    @return: a triple of (change count, 1 if any changes, 1 if non-reverted changes)
    """
    # initialize the random choice table
    choice_from = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    # the number of changes is poisson distributed
    change_count = numpy.random.poisson(branch_length)
    # note whether or not a (possibly reverted) change was observed
    observed_change = min(change_count, 1)
    # simulate the changes
    state = 0
    for i in range(change_count):
        state = random.choice(choice_from[state])
    # note whether a nonreverted change was observed
    nonreverted_change = min(state, 1)
    return change_count, observed_change, nonreverted_change

def hard_coded_analysis():
    branch_length = 5.0
    sequence_length = 1000
    nsequences = 1000
    estimate_triple_list = []
    column_headers = ('most.info', 'less.info', 'least.info')
    for i in range(nsequences):
        # sample sequence changes at three levels of informativeness
        sequence_changes = sample_sequence_changes(branch_length, sequence_length)
        # get a distance estimate for each level of informativeness
        estimate_triple = sample_distance(*sequence_changes)
        estimate_triple_list.append(estimate_triple)
    print RUtil.get_table_string(estimate_triple_list, column_headers)

def main():
    hard_coded_analysis()

if __name__ == '__main__':
    main()

