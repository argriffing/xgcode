"""
Utility functions for interfacing with R.
"""

import os
import tempfile
import unittest
import subprocess

import Util

#TODO replace assert with exceptions
#TODO move RTable into this module, with or without row headers

#FIXME
# This is a horrible hack.
# It is a result of the intersection of two bad things.
# First, there is no system-wide configuration file for the python scripts.
# Second, the servbio hpc does not have R on the path
# which is passed to the python script from the mobyle caller.
# Note that I used to try to use
# /usr/local/apps/R/xeon/2.9.0/bin/R
# but this version connected to an invalid readline library.
g_rlocations = (
        'R',
        '/usr/local/apps/R/em64t/R-2.11.1/bin/R')


class RError(Exception): pass

def mk_call_str(name, *args, **kwargs):
    args_v = [str(v) for v in args]
    kwargs_v = ['%s=%s' % kv for kv in kwargs.items()]
    arr = args_v + kwargs_v
    return '%s(%s)' % (name, ', '.join(arr))

def run(pathname):
    """
    Run the R script.
    Redirect .Rout to stderr.
    The return code is 0 for success and 1 for failure.
    The returned stdout and stderr are strings not files.
    @param pathname: name of the R script
    @return: (returncode, r_stdout, r_stderr)
    """
    for rlocation in g_rlocations:
        cmd = [
                rlocation, 'CMD', 'BATCH',
                '--vanilla', '--silent', '--slave',
                pathname, '/dev/stderr']
        try:
            proc = subprocess.Popen(cmd,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except OSError as e:
            continue
        break
    else:
        raise ValueError('could not find R')
    proc_stdout, proc_stderr = proc.communicate()
    return proc.returncode, proc_stdout, proc_stderr

def run_with_table(table, user_data, callback):
    """
    @param table: the table string
    @param user_data: typically a fieldstorage-like object
    @param callback: this callback is like f(user_data, table_filename)
    @return: the R output as a string
    """
    retcode, r_out, r_err = run_with_table_verbose(
            table, user_data, callback)
    if retcode:
        raise RError(r_err)
    return r_err

def run_with_table_verbose(table, user_data, callback):
    """
    @param table: the table string
    @param user_data: typically a fieldstorage-like object
    @param callback: this callback is like f(user_data, table_filename)
    @return: returncode, r_stdout, r_stderr
    """
    # Create a temporary data table file for R.
    f_temp_table = tempfile.NamedTemporaryFile(delete=False)
    f_temp_table.write(table)
    f_temp_table.close()
    # Create a temporary R script file.
    script_content = callback(user_data, f_temp_table.name)
    f_temp_script = tempfile.NamedTemporaryFile(delete=False)
    f_temp_script.write(script_content)
    f_temp_script.close()
    # Call R.
    retcode, r_out, r_err = run(f_temp_script.name)
    # To facilitate debugging, only delete temporary files if R was successful.
    if not retcode:
        # Delete the temporary data table file.
        os.unlink(f_temp_table.name)
        # Delete the temporary script file.
        os.unlink(f_temp_script.name)
    # Return the R results.
    return retcode, r_out, r_err

def run_plotter(table, user_script_content, device_name):
    """
    @param table: the table string
    @param user_script_content: script without header or footer
    @param device_name: an R device function name
    @return: returncode, r_stdout, r_stderr, image_data
    """
    temp_table_name = Util.create_tmp_file(table)
    temp_plot_name = Util.get_tmp_filename()
    script_content = '\n'.join([
        'my.table <- read.table("%s")' % temp_table_name,
        '%s("%s")' % (device_name, temp_plot_name),
        user_script_content,
        'dev.off()']) + '\n'
    temp_script_name = Util.create_tmp_file(script_content)
    retcode, r_out, r_err = run(temp_script_name)
    if retcode:
        image_data = None
    else:
        os.unlink(temp_table_name)
        os.unlink(temp_script_name)
        try:
            with open(temp_plot_name, 'rb') as fin:
                image_data = fin.read()
        except IOError as e:
            msg_a = 'could not open the plot image file'
            msg_b = ' that R was supposed to write'
            raise RError(msg_a + msg_b)
        os.unlink(temp_plot_name)
    return retcode, r_out, r_err, image_data

def get_table_string(M, column_headers):
    """
    Convert a row major rate matrix to a string representing an R table.
    @param M: a row major matrix
    @param column_headers: the labels of the data columns
    """
    # assert that the matrix is rectangular
    assert len(set(len(row) for row in M)) == 1
    # Assert that the number of columns
    # is the same as the number of column headers.
    assert len(M[0]) == len(column_headers)
    # assert that the headers are valid
    for header in column_headers:
        if '_' in header:
            msg_a = 'the header "%s" is invalid '
            msg_b = 'because it uses an underscore' % header
            raise ValueError(msg_a + msg_b)
    # define each line in the output string
    lines = []
    lines.append('\t'.join([''] + list(column_headers)))
    for i, row in enumerate(M):
        R_row = [float_to_R(value) for value in row]
        lines.append('\t'.join([str(i+1)] + R_row))
    return '\n'.join(lines)

def float_to_R(value):
    """
    Convert a python floating point value to a string usable by R.
    @param value: a floating point number
    @return: the R string representing the floating point number
    """
    if value == float('inf'):
        return 'Inf'
    else:
        return str(value)

class TestRUtil(unittest.TestCase):

    def test_mk_call_str(self):
        observed = mk_call_str('wat', 'x', 'y', z='foo', w='bar')
        expected = 'wat(x, y, z=foo, w=bar)'
        self.assertEqual(observed, expected)

if __name__ == '__main__':
    unittest.main()
