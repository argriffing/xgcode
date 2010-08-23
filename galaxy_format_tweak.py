"""
Hack the output data type for Galaxy.

This is inspired by emboss_format_corrector.py which
hacks some emboss output formats.
Tell galaxy to make this substitution by referring to this script
from the xml file which specifies the app interface.
Add this only when there is only one output
and this output is an image whose format is specified
by the imageformat parameter.
"""

def exec_after_process(
        app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """
    @param inp_data: possibly an input data dict
    @param out_data: possibly an output data dict
    @param param_dict: parameter dictionary
    """
    items = out_data.items()
    if len(items) != 1:
        raise ValueError('found more than one output')
    name, data = items[0]
    fmt = param_dict['imageformat']
    data = app.datatypes_registry.change_datatype(data, fmt)
    app.model.context.add(data)
    app.model.context.flush()
