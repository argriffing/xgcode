"""
Read and write linearly.
A sequence of objects is read from or written to a file one by one,
either backwards or forwards.
This is especially useful for the forward-backward algorithm which
computes position specific posterior hidden state distributions
given some hidden Markov model parameters and a sequence of observations.
"""

import unittest
from StringIO import StringIO

import iterutils

CLOSED = 0
READING = 1
WRITING = 2

class ConversionError(Exception): pass

class SequenceIOError(IOError): pass

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
    """
    def __init__(self):
        raise NotImplementedError()
    def read_forward(self):
        raise NotImplementedError()
    def read_backward(self):
        raise NotImplementedError()
    def write(self):
        raise NotImplementedError()
    def open_read(self):
        raise NotImplementedError()
    def open_write(self):
        raise NotImplementedError()
    def close(self):
        raise NotImplementedError()

class SequentialMemoryIO(SequentialIO):
    """
    Enable writing forward and reading forward and backward.
    Reading and writing is by value at the API level,
    and is also by value internally.
    Writing is only possible in the forward direction,
    but items may be read in the forward or backward direction.
    """
    def __init__(self):
        self.state = CLOSED
        self.arr = []
    def read_forward(self):
        if self.state != READING:
            raise SequenceIOError()
        for value in self.arr:
            yield value
    def read_backward(self):
        if self.state != READING:
            raise SequenceIOError()
        for value in reversed(self.arr):
            yield value
    def write(self, value):
        if self.state != WRITING:
            raise SequenceIOError()
        self.arr.append(value)
    def close(self):
        if self.state == CLOSED:
            raise SequenceIOError()
        self.state = CLOSED
    def open_read(self):
        if self.state != CLOSED:
            raise SequenceIOError()
        self.state = READING
    def open_write(self):
        if self.state != CLOSED:
            raise SequenceIOError()
        self.arr = []
        self.state = WRITING

class SequentialFileObjectIO(SequentialIO):
    """
    Enable writing forward and reading forward and backward.
    Reading and writing is by value at the API level,
    but is by line internally.
    Writing is only possible in the forward direction,
    but items may be read in the forward or backward direction.
    """
    def read_forward(self):
        if self.state != READING:
            raise SequenceIOError()
        for line in self.obj:
            line = line.strip()
            if line:
                yield self.converter.line_to_value(line)
    def read_backward(self):
        if self.state != READING:
            raise SequenceIOError()
        for line in iterutils.read_backwards(self.obj):
            line = line.strip()
            if line:
                yield self.converter.line_to_value(line)
    def write(self, value):
        if self.state != WRITING:
            raise SequenceIOError()
        line = self.converter.value_to_line(value)
        self.obj.write(line + '\n')

class SequentialDiskIO(SequentialFileObjectIO):
    def __init__(self, converter, filename):
        self.converter = converter
        self.filename = filename
        self.state = CLOSED
    def close(self):
        if self.state == CLOSED:
            raise SequenceIOError()
        self.obj.close()
        self.state = CLOSED
    def open_read(self):
        if self.state != CLOSED:
            raise SequenceIOError()
        self.obj = open(self.filename)
        self.state = READING
    def open_write(self):
        if self.state != CLOSED:
            raise SequenceIOError()
        self.obj = open(self.filename, 'w')
        self.state = WRITING

class SequentialStringIO(SequentialFileObjectIO):
    def __init__(self, converter, file_contents=''):
        self.converter = converter
        self.file_contents = file_contents
        self.state = CLOSED
    def close(self):
        if self.state == CLOSED:
            raise SequenceIOError()
        if self.state == WRITING:
            self.file_contents = self.obj.getvalue()
        self.obj.close()
        del self.obj
        self.state = CLOSED
    def open_read(self):
        if self.state != CLOSED:
            raise SequenceIOError()
        self.obj = StringIO(self.file_contents)
        self.state = READING
    def open_write(self):
        if self.state != CLOSED:
            raise SequenceIOError()
        self.obj = StringIO()
        self.state = WRITING


class TestLinearIO(unittest.TestCase):

    def stream_testing_helper(self, stream):
        """
        Test a linear stream.
        @param stream: a SequentialIO stream
        """
        expected = [1, 1, 2, 3, 5]
        # write the expected list to the stream
        stream.open_write()
        for x in expected:
            stream.write(x)
        stream.close()
        # attempt to get the expected list back from the stream
        stream.open_read()
        observed = list(stream.read_forward())
        stream.close()
        self.assertEqual(expected, observed)
        # attempt to get the reverse list back from the stream
        stream.open_read()
        observed = list(stream.read_backward())
        stream.close()
        self.assertEqual(list(reversed(expected)), observed)

    def test_string_stream(self):
        """
        Test StringIO storage for the linear stream.
        """
        converter = IntConverter()
        string_stream = SequentialStringIO(converter)
        self.stream_testing_helper(string_stream)

    def test_memory_stream(self):
        """
        Test memory storage for the linear stream.
        """
        memory_stream = SequentialMemoryIO()
        self.stream_testing_helper(memory_stream)


if __name__ == '__main__':
    unittest.main()
