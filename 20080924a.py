"""Sample estimated branch lengths using three different information sources.

The first method observes all of the changes
along the branch (infinite sites model).
For the second method multiple changes
are observed as a single change (infinite alleles model).
For the third method multiple changes
are observed as a single change unless they are reverted,
in which case they are observed as no change (Jukes-Cantor model).
"""

from StringIO import StringIO
import math
import random

import numpy as np

from SnippetUtil import HandlingError
import RUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('branch_length', 'branch length',
                1.0, low_exclusive=0.0),
            Form.Integer('sequence_length', 'sequence length',
                1000, low=1)]
    return form_objects

def get_form_out():
    return FormOut.Report()

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
    estimate_a, estimate_b, estimate_c = sample_distance(*sequence_changes)
    # begin the response
    out = StringIO()
    print >> out, 'distance estimates:'
    print >> out, 'using all change information:', estimate_a
    print >> out, 'without multiple change information:', estimate_b
    print >> out, 'without reverted change information:', estimate_c
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
    Return three values.
    First, the mean number of changes.
    Second, the change frequence,
    Third, the non-reverted change frequency
    The branch length is the expected number of changes
    along the branch at each site
    @param branch_length: the expected number of changes along the branch
    @param nsites: the number of sites in the sequence
    @return: the triple of results
    """
    samples = [sample_site_changes(branch_length) for i in range(nsites)]
    change_counts, observed_changes, nonreverted_changes = zip(*samples)
    mean_changes = sum(change_counts) / float(nsites)
    p_observed = sum(observed_changes) / float(nsites)
    p_nonreverted = sum(nonreverted_changes) / float(nsites)
    return mean_changes, p_observed, p_nonreverted

def sample_site_changes(branch_length):
    """
    Return three values.
    First, the change count.
    Second, 1 if any changes.
    Third, 1 if non-reverted changes.
    @param branch_length: the expected number of changes along the branch
    @return: the triple of results
    """
    # initialize the random choice table
    choice_from = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    # the number of changes is poisson distributed
    change_count = np.random.poisson(branch_length)
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
        sequence_changes = sample_sequence_changes(
                branch_length, sequence_length)
        # get a distance estimate for each level of informativeness
        estimate_triple = sample_distance(*sequence_changes)
        estimate_triple_list.append(estimate_triple)
    print RUtil.get_table_string(estimate_triple_list, column_headers)

def main():
    hard_coded_analysis()

if __name__ == '__main__':
    main()

