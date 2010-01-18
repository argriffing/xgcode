"""This module has some stuff related to the MAPP inference procedure by Stone and Sidow.
"""

import math

import numpy as np
import scipy

g_property_names = (
        'hydropathy',
        'polarity',
        'charge',
        'volume',
        'free energy alpha',
        'free energy beta')

# each row is an amino acid and each column is a property
g_property_array = np.array((
        (0.7667, 0.607, -0.0529, -1.2721, -0.91777, 0.67374),
        (1.001, 0.6887, -0.0529, -0.79141, 0.34486, -0.80639),
        ( -1.0077, -1.6003, -2.1704, -0.7286, 0.32539, 1.4671), 
        ( -1.0077, -1.3959, -2.1704, -0.069091, -1.6297, 0.69506),
        (1.1015, 1.0362, -0.0529, 1.175, -0.13071, -1.0111), 
        (0.0301, 0.4844, -0.0529, -1.9606, 1.1375, 0.77611), 
        ( -0.9073, -0.3331, 1.0058, 0.28844, 0.13628, -0.068461), 
        (1.6706, 0.9136, -0.0529, 0.61457, -0.52007, -1.4462),
        ( -1.1416, -1.5185, 2.0646, 0.66047, -0.3838, 0.4434), 
        (1.4362, 0.8522, -0.0529, 0.61457, -0.72031, -0.56752), 
        (0.8001, 0.9749, -0.0529, 0.52277, -0.85937, -0.46942), 
        ( -1.0077, -0.701, -0.0529, -0.65612, 0.42551, 1.1728), 
        ( -0.3716, 0.2391, -0.0529, -0.68995, 3.3151, 2.0899), 
        ( -1.0077, -0.5579, -0.0529, 0.061361, -0.4422, 0.54577), 
        ( -1.3425, -2.2338, 2.0646, 0.77643, -0.69528, 0.23012), 
        ( -0.1038, 0.4026, -0.0529, -1.2625, 0.4867, 0.27704), 
        ( -0.0703, 0.5252, -0.0529, -0.60781, 0.36433, -0.58885),
        (1.5702, 0.8114, -0.0529, -0.030439, -0.13071, -1.604), 
        ( -0.1373, 0.6683, -0.0529, 2.0906, -0.17521, -0.71255), 
        ( -0.2712, 0.1369, -0.0529, 1.2644, 0.069528, -1.0964)))

def get_mean(values):
    """
    @param values: some numbers
    @return: the mean of the numbers
    """
    n = float(len(values))
    return sum(values) / n

def get_variance(values):
    """
    @param values: some numbers
    @return: the variance of the numbers
    """
    n = float(len(values))
    mean = get_mean(values)
    variance = sum((x - mean)**2 for x in values) / (n-1)
    return variance

def get_standardized_property_array(property_array):
    """
    This is step 4 in the method section of the MAPP paper.
    @param property_array: each row is an amino acid and each column is a property
    @return: a standardized array with the same dimensions
    """
    # each row of P is a property
    P = property_array.T
    # get the standardized rows
    standardized_rows = []
    for row in P:
        mean = get_mean(row)
        variance = get_variance(row)
        stddev = math.sqrt(variance)
        standardized_row = [(x - mean) / stddev for x in row]
        standardized_rows.append(standardized_row)
    # return the standardized property array
    return np.array(standardized_rows).T

def estimate_aa_distribution(taxon_weights, taxon_amino_acid_indices):
    """
    This is a modification of part of step 3 in the methods section of the MAPP paper.
    @param taxon_weights: a sequence of weights
    @param taxon_amino_acid_indices: a sequence of amino acid indices
    @return: a sequence giving a probability for each amino acid
    """
    # validate the arguments
    if len(taxon_weights) != len(taxon_amino_acid_indices):
        sa = 'the number of taxon weights'
        sb = 'the number of taxon amino acid indices'
        raise ValueError('%s should be the same as %s' % (sa, sb))
    for i in taxon_amino_acid_indices:
        if i < 0 or i > 19:
            raise ValueError('each amino acid index should be between 0 and 19')
    epsilon = .0000000001
    if abs(sum(taxon_weights) - 1) > epsilon:
        raise ValueError('the taxon weights should sum to 1')
    for weight in taxon_weights:
        if weight < 0:
            raise ValueError('each taxon weight should be nonnegative')
    # define the weight to be given to the uniform distribution
    uniform_distribution_prior = .01
    # define the uniform distribution
    uniform_distribution = [.05 for i in range(20)]
    # define the maximum likelihood distribution
    ntaxa = len(taxon_weights)
    ml_distribution = [0]*20
    for i in taxon_amino_acid_indices:
        ml_distribution[i] += 1.0 / ntaxa
    # define the posterior distribution
    posterior_distribution = []
    for i in range(20):
        p = 0
        p += uniform_distribution_prior * uniform_distribution[i]
        p += (1 - uniform_distribution_prior) * ml_distribution[i]
        posterior_distribution.append(p)
    return posterior_distribution

def estimate_property_means(property_array, aa_distribution):
    """
    This is step 5a in the methods section of the MAPP paper.
    @param property_array: each row is an amino acid and each column is a (possibly standardized) property
    @param aa_distribution: a sequence of amino acid weights
    @return: a sequence of property means
    """
    assert np.shape(property_array) == (20, 6)
    assert len(aa_distribution) == 20
    property_means = []
    for property_values in zip(*property_array):
        mean = sum(weight * value for weight, value in zip(aa_distribution, property_values))
        property_means.append(mean)
    return property_means

def estimate_property_variances(property_array, aa_distribution):
    """
    This is step 5b in the methods section of the MAPP paper.
    @param property_array: each row is an amino acid and each column is a (possibly standardized) property
    @param aa_distribution: a sequence of amino acid weights
    @return: a sequence of property means
    """
    # validation of the arguments is done in the function that estimates the means
    property_means = estimate_property_means(property_array, aa_distribution)
    # estimate the variances
    property_variances = []
    property_array_transpose = zip(*property_array)
    for property_values, property_mean in zip(property_array_transpose, property_means):
        variance = sum(weight * (value - property_mean)**2 for weight, value in zip(aa_distribution, property_values))
        property_variances.append(variance)
    return property_variances

def get_deviations(estimated_property_means, estimated_property_variances, property_array):
    """
    This is step 6 in the methods section of the MAPP paper.
    @param estimated_property_means: a sequence of estimated property means
    @param estimated_property_variances: a sequence of estimated property variances
    @param property_array: a 20x6 array that gives (possibly standardized) properties for each amino acid
    @return: a 20x6 array that gives the deviation of each amino acid from the estimated distribution of each property
    """
    assert np.shape(property_array) == (20, 6)
    deviation_array = []
    for properties in property_array:
        amino_acid_deviations = []
        for estimated_mean, estimated_variance, value in zip(estimated_property_means, estimated_property_variances, properties):
            deviation = (value - estimated_mean) / estimated_variance
            amino_acid_deviations.append(deviation)
        deviation_array.append(amino_acid_deviations)
    return np.array(deviation_array)

def get_property_correlation_matrix(standardized_property_array):
    """
    This is step 7a in the methods section of the MAPP paper.
    @param standardized_property_array: a 20x6 array that gives standardized properties for each amino acid
    @return: a 6x6 symmetric correlation matrix
    """
    S = np.array(standardized_property_array)
    R = np.dot(S.T, S) / 19
    return R

def get_impact_scores(property_correlation_matrix, deviations):
    """
    This is step 7b in the methods section of the MAPP paper.
    @param property_correlation_matrix: a 6x6 matrix of property correlations
    @param deviations: a 20x6 array that gives the deviation of each amino acid from the estimated distribution of each property
    @return: a sequence of twenty impact scores
    """
    assert np.shape(property_correlation_matrix) == (6, 6)
    assert np.shape(deviations) == (20, 6)
    # get the inverse of the correlation matrix
    R_inv = np.linalg.inv(property_correlation_matrix)
    # get the impact scores
    impact_scores = []
    for deviation_row in deviations:
        squared_impact = 0
        for i, a in enumerate(deviation_row):
            for j, b in enumerate(deviation_row):
                squared_impact += a*b*R_inv[i][j]
        impact_scores.append(math.sqrt(squared_impact))
    return impact_scores

def get_p_value(impact_score, ntaxa, nproperties=6):
    """
    This is explained in the methods section on the null distribution of impact scores in the MAPP paper.
    The result is not available unless the number of taxa is greater than the number of properties.
    @param impact_score: the statistic that predicts the impact of the mutation
    @param ntaxa: the number of taxa represented in the column of interest
    @param nproperties: the number of amino acid properties used to calculate the impact score
    """
    assert ntaxa > nproperties
    fstat_numerator = (impact_score**2) * (ntaxa - nproperties)
    fstat_denominator = nproperties * (ntaxa + 1)
    fstat = fstat_numerator / fstat_denominator
    pvalue = 1.0 - scipy.stats.distributions.f.cdf(fstat, nproperties, ntaxa - nproperties)
    return pvalue


