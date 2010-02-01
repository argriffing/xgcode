"""
Read and write linearly.
A sequence of objects is read from or written to a file one by one,
either backwards or forwards.
This is especially useful for the forward-backward algorithm which
computes position specific posterior hidden state distributions
given some hidden Markov model parameters and a sequence of observations.
"""

import unittest
import StringIO

import Util

CLOSED = 0
READING = 1
WRITING = 2

class ConversionError(Exception): pass

class Converter:
    def line_to_value(self, line):
        try:
            return self._line_to_value(line)
        except ValueError:
            raise ConversionError(line)
    def value_to_line(self, value):
        try:
            return self._value_to_line(value)
        except ValueError:
            raise ConversionError(line)
    def _line_to_value(self, line):
        raise NotImplementedError()
    def _value_to_line(self, value):
        raise NotImplementedError()

class IntConverter(Converter):
    def _line_to_value(self, line):
        return int(line)
    def _value_to_line(self, value):
        return repr(value)

class FloatConverter(Converter):
    def _line_to_value(self, line):
        return float.fromhex(line)
    def _value_to_line(self, value):
        return value.hex()

class IntTupleConverter(Converter):
    def _line_to_value(self, line):
        row = line.split()
        return tuple(int(x) for x in row)
    def _value_to_line(self, value):
        return '\t'.join(repr(x) for x in value)

class FloatTupleConverter(Converter):
    def _line_to_value(self, line):
        row = line.split()
        return tuple(float.fromhex(x) for x in row)
    def _value_to_line(self, value):
        return '\t'.join(x.hex() for x in value)


class SequentialIO:
    """
    Enable writing forward and reading forward and backward.
    Reading and writing is by value at the API level,
    but is by line internally.
    Writing is only possible in the forward direction,
    but items may be read in the forward or backward direction.
    """
    def __init__(self):
        raise NotImplementedError()
    def read_forward(self):
        if self.state != READING:
            raise IOError('invalid action in the current state')
        for line in self.obj:
            line = line.strip()
            if line:
                yield self.converter.line_to_value(line)
    def read_backward(self):
        if self.state != READING:
            raise IOError('invalid action in the current state')
        for line in Util.read_backwards(self.obj):
            line = line.strip()
            if line:
                yield self.converter.line_to_value(line)
    def write(self, value):
        if self.state != WRITING:
            raise IOError('invalid action in the current state')
        line = self.converter.value_to_line(value)
        self.obj.write(line + '\n')

class SequentialDiskIO(SequentialIO):
    def __init__(self, converter, filename):
        self.converter = converter
        self.filename = filename
        self.state = CLOSED
    def close(self):
        if self.state == CLOSED:
            raise IOError('invalid action in the current state')
        self.obj.close()
        self.state = CLOSED
    def open_read(self, mode='r'):
        if self.state != CLOSED:
            raise IOError('invalid action in the current state')
        self.obj = open(self.filename)
        self.state = READING
    def open_write(self):
        if self.state != CLOSED:
            raise IOError('invalid action in the current state')
        self.obj = open(self.filename, 'w')
        self.state = WRITING

class SequentialStringIO(SequentialIO):
    def __init__(self, converter):
        self.converter = converter
        self.file_contents = ''
        self.state = CLOSED
    def close(self):
        if self.state == CLOSED:
            raise IOError('invalid action in the current state')
        if self.state == WRITING:
            self.file_contents = self.obj.getvalue()
        self.obj.close()
        del self.obj
        self.state = CLOSED
    def open_read(self):
        if self.state != CLOSED:
            raise IOError('invalid action in the current state')
        self.obj = StringIO.StringIO(self.file_contents)
        self.state = READING
    def open_write(self):
        if self.state != CLOSED:
            raise IOError('invalid action in the current state')
        self.obj = StringIO.StringIO()
        self.state = WRITING


class TestLinearIO(unittest.TestCase):

    def test_string_stream(self):
        """
        Test StringIO storage.
        """
        expected = [1, 1, 2, 3, 5]
        # initialize the stream
        o_converter = IntConverter()
        o_stream = SequentialStringIO(o_converter)
        # write the expected list to the stream
        o_stream.open_write()
        for x in expected:
            o_stream.write(x)
        o_stream.close()
        # attempt to get the expected list back from the stream
        o_stream.open_read()
        observed = list(o_stream.read_forward())
        o_stream.close()
        self.assertEqual(expected, observed)
        # attempt to get the reverse list back from the stream
        o_stream.open_read()
        observed = list(o_stream.read_backward())
        o_stream.close()
        self.assertEqual(list(reversed(expected)), observed)


if __name__ == '__main__':
    unittest.main()
