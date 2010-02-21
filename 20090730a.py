"""Explore some remedial statistics concepts. [UNFINISHED]
"""

from StringIO import StringIO
import random

import numpy as np

from SnippetUtil import HandlingError
import Form

def get_form():
    """
    @return: the body of a form
    """
    return []

def sample_stuff():
    """
    Make some sample phenotypes and two sample SNPs.
    Each of the three returned values is a list of length nstrains.
    @return: phenotypes, snpa, snpb
    """
    mu = 1.0
    alpha = 1.0
    beta = 2.0
    gamma = 3.0
    nstrains = 40
    snpa = [random.randrange(2) for i in range(nstrains)]
    snpb = [random.randrange(2) for i in range(nstrains)]
    phen = [mu + alpha*ma + beta*mb + random.gauss(0, 0.1) for ma, mb in zip(snpa, snpb)]
    return phen, snpa, snpb

def get_estimates_a(phen, snpa, snpb):
    """
    Do the least squares all in one step.
    @param phen: the phenotypes
    @param snpa: zero or one for genotypes in a snp
    @param snpb: zero or one for genotypes in a snp
    @return: [mu_hat, alpha_hat, beta_hat]
    """
    nstrains = len(phen)
    X = np.array([np.ones(nstrains), snpa, snpb]).T
    Y = np.array(phen)
    x, residues, rank, s = np.linalg.lstsq(X, Y)
    return x.tolist()

def get_estimates_b(phen, snpa, snpb):
    """
    Remove the additive effects one by one.
    @param phen: the phenotypes
    @param snpa: zero or one for genotypes in a snp
    @param snpb: zero or one for genotypes in a snp
    @return: [mu_hat, alpha_hat, beta_hat]
    """
    nstrains = len(phen)
    # copy the phenotype list to prepare to remove effects
    ph = phen[:]
    # estimate the first additive effect and remove it from the phenotypes
    alpha_hat = np.mean([x for x, a in zip(ph, snpa) if a]) - np.mean([x for x, a in zip(ph, snpa) if not a])
    ph = [x - alpha_hat*a for x, a in zip(ph, snpa)]
    # estimate the second additive effect and remove it from the phenotype
    beta_hat = np.mean([x for x, b in zip(ph, snpb) if b]) - np.mean([x for x, b in zip(ph, snpb) if not b])
    ph = [x - beta_hat*b for x, b in zip(ph, snpb)]
    # estimate the mean and remove it from the phenotypes
    mu_hat = np.mean(ph)
    ph = [x - mu_hat for x in ph]
    # return the estimates
    return [mu_hat, alpha_hat, beta_hat]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # create the response
    out = StringIO()
    phen, snpa, snpb = sample_stuff()
    print >> out, get_estimates_a(phen, snpa, snpb)
    print >> out, get_estimates_b(phen, snpa, snpb)
    # return the response
    response_text = out.getvalue().strip()
    return [('Content-Type', 'text/plain')], response_text

