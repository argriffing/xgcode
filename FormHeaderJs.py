"""
Construct the javascript that goes into the html headers.

This includes four categories of javascript.
The first category is the reference to external javascript resources.
This includes the line or two to tell the client where to find mathjax.
The second category is javascript boilerplate which could
eventually go into its own external javascript file.
This includes utility functions for the javascript code behind the presets,
for example programmatically selecting radio buttons.
The third category is slightly less general javascript code
which relies on hardcoded identifiers of html elements,
but which is still general enough to not change between scripts.
The fourth category is script-specific preset actions.
These are functions that are called when the user clicks preset buttons,
and they are written programmatically by the server side python code
for the individual web scripts.
"""

from StringIO import StringIO

g_cat_2_default = """
<!-- these are javascript functions from the internet -->
<script type="text/javascript">
// this next function is from the internet
// set the radio button with the given value as being checked
// do nothing if there are no radio buttons
// if the given value does not exist, all the radio buttons
// are reset to unchecked
function setCheckedValue(radioObj, newValue) {
	if(!radioObj)
		return;
	var radioLength = radioObj.length;
	if(radioLength == undefined) {
		radioObj.checked = (radioObj.value == newValue.toString());
		return;
	}
	for(var i = 0; i < radioLength; i++) {
		radioObj[i].checked = false;
		if(radioObj[i].value == newValue.toString()) {
			radioObj[i].checked = true;
		}
	}
}
</script>
""".strip()

def get_header_script_text(form_objects=None, presets=None):
    """
    """
    out = StringIO()
    print >> out, _get_cat_1_text()
    if presets:
        print >> out, _get_cat_2_text()
        print >> out, _get_cat_3_text()
        print >> out, _get_cat_4_text(form_objects, presets)
    return out.getvalue().rstrip()

def _get_cat_1_text():
    """
    References to external javascript resources.
    @return: text
    """
    out = StringIO()
    mathjax_js = 'cdn.mathjax.org/mathjax/latest/MathJax.js'
    mathjax_params = 'config=TeX-AMS-MML_HTMLorMML'
    mathjax_url = 'http://' + mathjax_js + '?' + mathjax_params
    print >> out, '<!-- these are links to external javascript resources -->'
    print >> out, '<script type="text/javascript" src="%s">' % mathjax_url
    print >> out, '</script>'
    return out.getvalue().rstrip()

def _get_cat_2_text():
    return g_cat_2_default

def _get_cat_3_text():
    out = StringIO()
    print >> out, '<script type="text/javascript">'
    #
    print >> out, 'function wsfSetInner(id, value) {',
    print >> out, 'document.getElementById(id).innerHTML = value;',
    print >> out, '}'
    #
    print >> out, 'function wsfSetRadio(name, value) {',
    print >> out, 'setCheckedValue(',
    print >> out, 'document.forms["mainform"].elements[name], value);',
    print >> out, ')',
    print >> out, '}'
    #
    print >> out, '</script>'
    return out.getvalue().rstrip()

def _get_cat_4_text(form_objects, presets):
    return ''

