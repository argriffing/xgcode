"""
Web interface stuff.
"""

from StringIO import StringIO

import smallutil


class HandlingError(Exception): pass

def get_state_value_pair(line):
    """
    This function is for reading data from a web form.
    @param line: a line of text that looks something like 'A : 0.245'
    @return: a state string and a value string
    """
    if not line.strip():
        raise HandlingError('found a blank line when colon separated values were expected')
    if line.count(':') != 1:
        raise HandlingError('expected each line to have a single ":" separator')
    state, value = line.split(':')
    return state.strip(), value.strip()

def get_generic_dictionary(lines, key_name, value_name, valid_values):
    """
    This high level function returns a mapping from keys to values.
    Each key is expected to have exactly one value.
    @param lines: lines of text that each look something like 'Human : G'
    @param key_name: something like 'species name'
    @param value_name: something like 'nucleotide'
    @param valid_values: a set of valid values, or None if any state is valid
    @return: a dictionary mapping each state to a value
    """
    # convert the line sequence to a list of non-empty lines
    lines = list(line for line in lines if line.strip())
    # make sure that at least one line is non-empty
    if not lines:
        raise HandlingError('no %s to %s mapping was specified' % (key_name, value_name))
    lines = list(lines)
    d = {}
    for line in lines:
        key, value = get_state_value_pair(line)
        if valid_values is not None:
            if value not in valid_values:
                raise HandlingError('invalid %s: %s' % (value_name, value))
        if key in d:
            raise HandlingError('duplicate %s: %s' % (key_name, key))
        d[key] = value
    return d

def get_dictionary(dictionary_string, state_name, value_name, required_states):
    """
    This high level function returns a mapping from states to floating point values.
    Each required state is expected to have exactly one floating point value.
    @param dictionary_string: a string of lines of text that each look something like 'A : 0.245'
    @param state_name: something like 'amino acid'
    @param value_name: something like 'energy'
    @param required_states: a set of required states
    @return: a dictionary mapping each state to a value
    """
    # convert the multi-line string to a list of non-empty lines
    lines = smallutil.get_stripped_lines(StringIO(dictionary_string))
    # make sure that at least one line is non-empty
    if not lines:
        raise HandlingError('no %s to %s mapping was specified' % (key_name, value_name))
    d = {}
    for line in lines:
        state, value_string = get_state_value_pair(line)
        if state not in required_states:
            raise HandlingError('invalid %s: %s' % (state_name, state))
        try:
            value = float(value_string)
        except ValueError:
            raise HandlingError('this %s could not be interpreted as a number: %s' % (value_name, value_string))
        d[state] = value
    # assert that each state was assigned a value
    if len(d) < len(required_states):
        raise HandlingError('one or more %s was not assigned a %s' % (state_name, value_name))
    return d

def get_distribution(distribution_string, state_name, valid_states):
    """
    This high level function returns a state distribution.
    Each valid state is expected to have exactly one weight.
    The weights in the returned dictionary are normalized to sum to one.
    @param distribution_string: a string of lines of text that each look something like 'A : 0.245'
    @param state_name: something like 'amino acid'
    @param valid_states: a set of valid states
    @return: a dictionary mapping each state to a probability
    """
    if not distribution_string:
        raise HandlingError('no %s distribution was specified' % state_name)
    state_to_weight = {}
    for line in smallutil.stripped_lines(StringIO(distribution_string)):
        state, weight = get_weight_pair(line, state_name, valid_states)
        if state in state_to_weight:
            raise HandlingError('duplicate %s: %s' % (state_name, state))
        state_to_weight[state] = weight
    if len(state_to_weight) < len(valid_states):
        raise HandlingError('one or more %s was not assigned a weight' % state_name)
    total_weight = float(sum(state_to_weight.values()))
    if not total_weight:
        raise HandlingError('each %s weight is zero' % state_name)
    for state in state_to_weight:
        state_to_weight[state] /= total_weight
    return state_to_weight

def get_weight_pair(line, state_name, valid_states):
    """
    This function is for reading data from a web form.
    @param line: a stripped line of text that looks something like 'A : 0.245'
    @param state_name: something like 'amino acid'
    @param valid_states: a set of valid states
    @return: a state and a real valued weight
    """
    state, weight_string = get_state_value_pair(line)
    if state not in valid_states:
        raise HandlingError('invalid %s: %s' % (state_name, state))
    try:
        weight = float(weight_string)
    except ValueError:
        raise HandlingError('this weight could not be interpreted as a number: %s' % weight_string)
    if weight < 0:
        raise HandlingError('found a negative weight: %s' % weight_string)
    return state, weight

def get_sample_html():
    """
    @return: the contents of a placeholder html file
    """
    target = 'http://www.youtube.com/watch?v=oHg5SJYRHA0'
    sio = StringIO()
    print >> sio, '<html>'
    print >> sio, '<head>'
    print >> sio, '<META HTTP-EQUIV=REFRESH CONTENT="1; URL=%s"/>' % target
    print >> sio, '</head>'
    print >> sio, '<body>'
    print >> sio, '</body>'
    print >> sio, '</html>'
    return sio.getvalue()

def docstring_to_title(docstring):
    """
    @param docstring: something like __doc__
    @return: the first line of the docstring as a title, or None
    """
    # get lines of text without whitespace between lines
    lines = smallutil.get_stripped_lines(StringIO(docstring))
    if lines:
        return lines[0]
    else:
        return None

def docstring_to_html(docstring):
    """
    Convert the docstring to an html header.
    @param docstring: something like __doc__
    """
    # get lines of text without whitespace between lines
    lines = smallutil.get_stripped_lines(StringIO(docstring))
    # no docstring
    if not lines:
        return ''
    # a single line docstring
    if len(lines) == 1:
        return lines[0]
    # multiple lines
    # the first line will be separated by an empty line
    arr = []
    arr.append(lines[0])
    arr.append('<br/><br/>')
    arr.extend(lines[1:])
    return '\n'.join(arr)
