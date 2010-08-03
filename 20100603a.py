"""Create a fungus .ind file for temperature.

Create an .ind temperature file from a .hud and a amdS_PCA_Info.csv file.
The .hud file provides the names of the OTUs.
The amdS_PCA_Info.csv file provides the 'case-control' status,
representing binarized location, temperature, or precipitation.
The output file is in eigenstrat format.
"""

from StringIO import StringIO
import csv

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import hud

g_tags = ['ztools']

g_default_hud_string = """
IC1 1 1 1 0
IC2 1 1 1 0
IC3 1 0 1 0
""".strip()

g_default_info_lines = [
        '"IC","Haplo","Location","Temp (C)","Precip (mm)","Species",'
            '"B1","B2","G1","G2","OMST"',
        '"1","H42","GA","15","600","Ap","+","+","+","+","-"',
        '"2","H42","GA","30","700","Ap","+","+","+","+","-"',
        '"3","*","GA","45","800","Ap","+","+","+","+","-"']

g_default_info_string = '\n'.join(g_default_info_lines)

class DataRowError(Exception):
    def __init__(self, row, e):
        lines = ['Error in data row ' + str(row), str(e)]
        msg = '\n'.join(lines)
        Exception.__init__(self, msg)

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'a list of OTUs names and binary character vectors',
                g_default_hud_string),
            Form.MultiLine('info',
                'amdS_PCA_Info.csv lines',
                g_default_info_string),
            Form.Float('threshold',
                    'temperature threshold (C)',
                    '22.0'),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.EigenstratInd('temperature')

def get_response_content(fs):
    out = StringIO()
    # extract name order from the .hud file
    names, hud_data = hud.decode(fs.hud.splitlines())
    # read the csv file
    rows = list(csv.reader(Util.get_stripped_lines(fs.info.splitlines())))
    header, data_rows = rows[0], rows[1:]
    cases, controls = get_temperature_info(data_rows, fs.threshold)
    # write the .ind file contents
    for name in names:
        gender = 'U'
        classification = 'Ignore'
        if name in cases:
            classification = 'Case'
        elif name in controls:
            classification = 'Control'
        row = [name, gender, classification]
        print >> out, '\t'.join(row)
    return out.getvalue()

def get_temperature_info(data_rows, threshold):
    """
    Asterisk is missing data.
    @param data_rows: rows of string elements
    @param threshold: temperature threshold in Celcius
    @return: temperature case and control OTU sets
    """
    cases = set()
    controls = set()
    for row in data_rows:
        try:
            otu = 'IC' + row[0]
            temperature = row[3]
            if temperature == '*':
                continue
            if float(temperature) < threshold:
                cases.add(otu)
            else:
                controls.add(otu)
        except Exception, e:
            raise DataRowError(row, e)
    return cases, controls
