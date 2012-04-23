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

import Form

g_cat_2_default = """
<!-- these are javascript functions from the internet -->
<script type="text/javascript">
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

g_cat_3_default = """
<!-- these are custom javascript functions with hardcoded references -->
<script type="text/javascript">
function wsfSetInner(id, value) {
    document.getElementById(id).innerHTML = value;
}
function wsfSetValue(id, value) {
    document.getElementById(id).value = value;
}
function wsfSetRadio(name, value) {
    var myobj = document.forms["mainform"].elements[name];
    setCheckedValue(myobj, value);
}
function wsfSetChecks(name, value_to_checked) {
    var myobj = document.forms["mainform"].elements[name];
    for (var k in myobj) {
        myobj[k].checked = false;
    }
    for (var k in value_to_checked) {
        myobj[k].checked = value_to_checked[k];
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
    return g_cat_3_default

def _collection_to_javascript_literal(coll):
    out = StringIO()
    print >> out, '{',
    for elem in coll:
        print >> out, elem, ': true,',
    print >> out, '}'
    return out.getvalue().strip()

def _get_cat_4_text(form_objects, presets):
    label_to_form_object = dict((o.label, o) for o in form_objects)
    out = StringIO()
    print >> out, '<!-- these are hardcoded presets -->'
    print >> out, '<script type="text/javascript">'
    for i, preset in enumerate(presets):
        print >> out, 'function on_wsf_preset%d() {' % i
        for k, v in preset.d.items():
            # check the form object corresponding to this item
            # FIXME this is horrible typechecking that is not pythonic
            form_object = label_to_form_object[k]
            if isinstance(form_object, Form.RadioGroup):
                print >> out, 'wsfSetRadio("%s", "%s");' % (k, v)
            elif isinstance(form_object, Form.CheckGroup):
                js_literal = _collection_to_javascript_literal(v)
                print >> out, 'wsfSetChecks("%s", %s);' % (k, js_literal)
            elif isinstance(form_object, Form.Sequence):
                print >> out, 'wsfSetInner("%s", "%s");' % (k, '\\n'.join(v))
            elif any([
                isinstance(form_object, Form.Integer),
                isinstance(form_object, Form.Float),
                isinstance(form_object, Form.SingleLine)]):
                print >> out, 'wsfSetValue("%s", "%s");' % (k, v)
            else:
                # assume we want to set the inner html
                print >> out, 'wsfSetInner("%s", "%s");' % (k, v)
        print >> out, '}'
    print >> out, '</script>'
    return out.getvalue().rstrip()

