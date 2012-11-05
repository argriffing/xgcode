
import numpy
import scipy.special
import mpmath

def absolute_error(expected, observed):
    return abs(expected - observed)

def relative_error(expected, observed):
    return absolute_error(expected, observed) / abs(expected)

def float_mp_hyp1f1(a, b, x):
    try:
        return float(mpmath.hyp1f1(a, b, x))
    except TypeError as e:
        return numpy.nan

def float_mp_hyp2f0(a1, a2, x):
    try:
        return float(mpmath.hyp2f0(a1, a2, x))
    except TypeError as e:
        return numpy.nan

def sp_hyp2f0_type1(a1, a2, x):
    convergence_type = 1
    y, err = scipy.special.hyp2f0(a1, a2, x, convergence_type)
    return y

def sp_hyp2f0_type2(a1, a2, x):
    convergence_type = 2
    y, err = scipy.special.hyp2f0(a1, a2, x, convergence_type)
    return y

def sp_hyp2f0_type1_combo(a1, a2, x):
    if x < 0:
        return sp_hyp2f0_type1(a1, a2, x)
    else:
        return (-1/x)**a1 * scipy.special.hyperu(a1, 1+a1-a2, (-1/x))

def show_error(a, expected, observed):
    print len(a), len(expected), len(observed)
    for ai, ei, oi in zip(a, expected, observed):
        print numpy.array([ai, ei, oi]), relative_error(ei, oi)

def main():
    mp_hyp1f1 = numpy.vectorize(float_mp_hyp1f1)
    mp_hyp2f0 = numpy.vectorize(float_mp_hyp2f0)
    #a_pos = numpy.power(10, numpy.linspace(-10, 10, 21))
    a_pos = numpy.power(10, numpy.linspace(-5, 5, 21))
    a_neg = -a_pos
    a = numpy.hstack([a_neg, a_pos])
    print len(a_neg)
    print len(a_pos)
    print len(a)
    #
    print a
    #
    print absolute_error(
            mp_hyp1f1(0.5, 1.5, a),
            scipy.special.hyp1f1(0.5, 1.5, a))
    print relative_error(
            mp_hyp1f1(0.5, 1.5, a),
            scipy.special.hyp1f1(0.5, 1.5, a))
    #
    print absolute_error(
            mp_hyp2f0(1.0, 0.5, a),
            sp_hyp2f0_type1(1.0, 0.5, a))
    print relative_error(
            mp_hyp2f0(1.0, 0.5, a),
            sp_hyp2f0_type1(1.0, 0.5, a))
    #
    print absolute_error(
            mp_hyp2f0(1.0, 0.5, a),
            sp_hyp2f0_type2(1.0, 0.5, a))
    print relative_error(
            mp_hyp2f0(1.0, 0.5, a),
            sp_hyp2f0_type2(1.0, 0.5, a))
    print
    print 'hyp1f1:'
    expected = mp_hyp1f1(0.5, 1.5, a)
    observed = scipy.special.hyp1f1(0.5, 1.5, a)
    show_error(a, expected, observed)
    print
    print 'hyp2f0 type 1 convergence:'
    expected = mp_hyp2f0(1.0, 0.5, a)
    observed = sp_hyp2f0_type1(1.0, 0.5, a)
    show_error(a, expected, observed)
    print
    print 'hyp2f0 type 2 convergence:'
    expected = mp_hyp2f0(1.0, 0.5, a)
    observed = sp_hyp2f0_type2(1.0, 0.5, a)
    show_error(a, expected, observed)
    print


if __name__ == '__main__':
    main()

