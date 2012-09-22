"""
Check confounding of h and s using diffusion theory for large populations.

The confounding is checked with respect to fixation probability.
This script assumes a large background population size
with a single mutant allele and no subsequent mutation.
It uses the notation in the paper
"Effects of dominance on the probability of fixation of a mutant allele."
Christina T. L. Chen and Quo-Shin Chi and Stanley A. Sawyer.
2008, J. Math. Biol.
Fitnesses are
AA: 1 + sigma
aA: 1 + h * sigma
aa: 1
The allele 'a' is the background allele.
The scaled selection coefficient is
s = 2*N*sigma
where N is the effective population size.
"""

from StringIO import StringIO
import cgi

import Form
import FormOut
import kimura

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size',
                10000, low=1),
            Form.Sequence('s_values', 'selection values',
                ('1', '2', '4')),
            Form.Sequence('h_values', 'dominance values',
                ('0.5', '0.25')),
            ]

def get_form_out():
    return FormOut.Html()

def get_html_table_row(arr):
    out = StringIO()
    print >> out, '<tr>'
    for value in arr:
        print >> out, '<td>%s</td>' % cgi.escape(str(value))
    print >> out, '</tr>'
    return out.getvalue().rstrip()

def get_response_content(fs):
    out = StringIO()
    #
    nmutants = 1
    N_hap = 2 * fs.N_diploid
    if nmutants >= N_hap:
        raise Exception('too many initial mutant alleles')
    p0 = nmutants / float(N_hap)
    h_values = [float(h) for h in fs.h_values]
    s_values = [float(s) for s in fs.s_values]
    #
    print >> out, '<html><body><table border="1" cellpadding="10">'
    #
    headers = (
            'h', 's',
            'fAA', 'faA', 'faa',
            'pfix',
            )
    print >> out, get_html_table_row(headers)
    #
    for h in h_values:
        for s in s_values:
            sigma = s / float(fs.N_diploid)
            #
            fAA = 1.0 + sigma
            faA = 1.0 + h * sigma
            faa = 1.0
            #
            pfix = kimura.get_fixation_probability_chen(p0, s, h)
            values = (
                    h, s,
                    fAA, faA, faa,
                    pfix,
                    )
            print >> out, get_html_table_row(values)
    print >> out, '</table></body></html>'
    return out.getvalue()

