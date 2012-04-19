"""
Make an example interactive embedded svg using pydot.
"""

from StringIO import StringIO
import string

import pydot
import lxml
from lxml import etree

import Form
import FormOut
import Util
import graph

g_default_name_pairs = [
        ('beta1', 'gamma1'),
        ('beta1', 'gamma2'),
        ('beta2', 'gamma3'),
        ('alpha', 'beta1'),
        ('alpha', 'beta2'),
        ('alpha', 'beta3'),
        ('beta2', 'gamma4'),
        ('beta3', 'gamma5'),
        ('beta3', 'gamma6'),
        ('alpha', 'gamma6'),
        ]

g_script_init = """
// some constants
var NODE_PLAIN = 'node-plain';
var NODE_IMPLIED = 'node-implied';
var NODE_SELECTED = 'node-selected';
var EDGE_PLAIN = 'edge-plain';
var EDGE_IMPLIED = 'edge-implied';
// map the node state to the fill color
var node_state_to_fill = {};
node_state_to_fill[NODE_PLAIN] = 'white';
node_state_to_fill[NODE_IMPLIED] = 'white';
node_state_to_fill[NODE_SELECTED] = 'cyan';
// map the node state to the stroke color
var node_state_to_stroke = {};
node_state_to_stroke[NODE_PLAIN] = 'black';
node_state_to_stroke[NODE_IMPLIED] = 'blue';
node_state_to_stroke[NODE_SELECTED] = 'blue';
// map the node state to the stroke width
var node_state_to_stroke_width = {};
node_state_to_stroke_width[NODE_PLAIN] = '1px';
node_state_to_stroke_width[NODE_IMPLIED] = '2px';
node_state_to_stroke_width[NODE_SELECTED] = '2px';
// map the edge state to the stroke color
var edge_state_to_stroke = {};
edge_state_to_stroke[EDGE_PLAIN] = 'black';
edge_state_to_stroke[EDGE_IMPLIED] = 'blue';
// map the edge state to the stroke width
var edge_state_to_stroke_width = {};
edge_state_to_stroke_width[EDGE_PLAIN] = '1px';
edge_state_to_stroke_width[EDGE_IMPLIED] = '2px';
"""

g_script_guts = """
// initialize nodes to plain states
var node_to_state = {};
for (var i=0; i<topo_sorted_node_ids.length; i++) {
	node_to_state[topo_sorted_node_ids[i]] = NODE_PLAIN;
}
// initialize edges to plain states
var edge_to_state = {};
for (edge_id in edge_id_to_sink_id) {
	edge_to_state[edge_id] = EDGE_PLAIN;
}
function redraw() {
	var svg = document.getElementById('omgsvg');
	for (node_id in node_to_state) {
		var state = node_to_state[node_id];
		var elem = svg.getElementById(node_id);
		var children = elem.getElementsByTagName('polygon');
		for (var i=0; i<children.length; i++) {
			var child = children[i];
			//for (child in children) {
			child.setAttribute('fill',
					node_state_to_fill[state]);
			child.setAttribute('stroke',
					node_state_to_stroke[state]);
			child.setAttribute('stroke-width',
					node_state_to_stroke_width[state]);
		}
	}
	for (edge_id in edge_to_state) {
		var state = edge_to_state[edge_id];
		var elem = svg.getElementById(edge_id);
		var children = elem.getElementsByTagName('path');
		for (var i=0; i<children.length; i++) {
			var child = children[i];
			child.setAttribute('stroke',
					edge_state_to_stroke[state]);
			child.setAttribute('stroke-width',
					edge_state_to_stroke_width[state]);
		}
		var children = elem.getElementsByTagName('polygon');
		for (var i=0; i<children.length; i++) {
			var child = children[i];
			child.setAttribute('stroke',
					edge_state_to_stroke[state]);
			child.setAttribute('stroke-width',
					edge_state_to_stroke_width[state]);
		}
	}
}
function recompute_implications() {
	// the order is important here
	for (var i=0; i<topo_sorted_node_ids.length; i++) {
		// get the node id
		var node_id = topo_sorted_node_ids[i];
		// first go through the edges if any
		var found_edge = false;
		var edge_id_array = node_id_to_edge_ids[node_id];
		for (var j=0; j<edge_id_array.length; j++) {
			var edge_id = edge_id_array[j];
			var sink_id = edge_id_to_sink_id[edge_id];
			if (node_to_state[sink_id] == NODE_PLAIN) {
				edge_to_state[edge_id] = EDGE_PLAIN;
			} else {
				edge_to_state[edge_id] = EDGE_IMPLIED;
				found_edge = true;
			}
		}
		// next determine the node state
		if (node_to_state[node_id] == NODE_SELECTED) {
			;
		} else if (found_edge) {
			node_to_state[node_id] = NODE_IMPLIED;
		} else {
			node_to_state[node_id] = NODE_PLAIN;
		}
	}
}
function dosomethingcool(myevt) {
	var pid = myevt.target.parentNode.id;
	if (node_to_state[pid] == NODE_SELECTED) {
		node_to_state[pid] = NODE_PLAIN;
	} else {
		node_to_state[pid] = NODE_SELECTED;
	}
	recompute_implications();
	redraw();
}
function coolsubmission() {
    // accumulate the selected value ids into a single string
    var accum = '';
	for (node_id in node_to_state) {
		var state = node_to_state[node_id];
        if (state == NODE_SELECTED) {
            var name = node_id_to_name[node_id];
            if (accum == '') {
                accum = name;
            } else {
                accum = accum + ' ' + name;
            }
        }
    }
    // set the extraselections value to the selected values
    var field = document.getElementById('extraselections');
    field.setAttribute('value', accum);
}
"""

g_html_form = """
<form id="myform" action="http://does.not.exist" method="get"
onsubmit="coolsubmission();">
<input type="hidden" id="extraselections" name="extraselections" value="wat"/>
<input type="submit" name="submit" value="simulated form submission"/>
</form>
"""

def get_form():
    """
    @return: the body of a form
    """
    default_list = [' '.join(p) for p in g_default_name_pairs]
    form_objects = [
            Form.Float('width_inches', 'width inches',
                '3.0', low_inclusive=1, high_inclusive=10),
            Form.Sequence('strlist', 'directed edges', default_list)]
    return form_objects

def get_form_out():
    return FormOut.Html('interactive')

def modify_tree(tree, name_to_node, pair_to_edge):
    name_to_node_id = {}
    name_pair_to_edge_id = {}
    for elem in tree.iter():
        parent = elem.getparent()
        #print elem.tag, elem.text
        if elem.tag == '{http://www.w3.org/2000/svg}svg':
            elem.set('id', 'omgsvg')
        elif elem.tag == '{http://www.w3.org/2000/svg}g':
            if elem.get('class') == 'node':
                elem.set('cursor', 'pointer')
                elem.set('onclick', 'dosomethingcool(evt);')
        elif elem.tag == '{http://www.w3.org/2000/svg}text':
            elem.set('pointer-events', 'none')
        elif elem.tag == '{http://www.w3.org/2000/svg}polygon':
            if parent.get('class') == 'node':
                elem.set('fill', 'white')
            elif parent.get('class') == 'edge':
                elem.set('fill', 'black')
        elif elem.tag == '{http://www.w3.org/2000/svg}title':
            if parent.get('class') == 'node':
                name_to_node_id[elem.text] = parent.get('id')
            elif parent.get('class') == 'edge':
                source_name, sink_name = elem.text.split('->')
                name_pair = (source_name, sink_name)
                name_pair_to_edge_id[name_pair] = parent.get('id')
    return name_to_node_id, name_pair_to_edge_id

def _make_graph(pairs, width_inches, rankdir):
    """
    @param pairs: directed edges between node labels
    @param width_inches: width in inches
    @param rankdir: 'LR' or 'TD'
    @return: tree, name_to_node, pair_to_edge
    """
    # initialize the graph
    g = pydot.Dot(
            graph_type='digraph',
            size='%s, 8' % width_inches,
            rankdir=rankdir,
            ratio='compress')
    # create the nodes
    name_to_node = {}
    for pair in pairs:
        for name in pair:
            if name not in name_to_node:
                node = pydot.Node(name, shape='rect')
                g.add_node(node)
                name_to_node[name] = node
    # create the edges
    pair_to_edge = {}
    for pair in pairs:
        a, b = pair
        edge = pydot.Edge(name_to_node[a], name_to_node[b])
        g.add_edge(edge)
        pair_to_edge[pair] = edge
    # do the physical layout and create the svg string
    tmp_path = Util.create_tmp_file(data=None, prefix='tmp', suffix='.svg')
    g.write_svg(tmp_path)
    with open(tmp_path) as fin:
        svg_str = fin.read()
    # parse the svg as xml except not the first six lines
    svg_str = '\n'.join(svg_str.splitlines()[6:])
    tree = etree.parse(StringIO(svg_str))
    return tree, name_to_node, pair_to_edge

def make_graph(pairs, width_inches):
    """
    Choose the rankdir that gives more area.
    @param pairs: directed edges between node labels
    @param width_inches: width in inches
    @return: tree, name_to_node, pair_to_edge
    """
    quads = []
    for rankdir in ('TB', 'LR'):
        tree, name_to_node, pair_to_edge = _make_graph(
                pairs, width_inches, rankdir)
        area = None
        for elem in tree.iter():
            if elem.tag == '{http://www.w3.org/2000/svg}svg':
                widthpt = elem.get('width')
                heightpt = elem.get('height')
                width = float(widthpt[:-2])
                height = float(heightpt[:-2])
                area = width * height
                break
        quads.append((area, tree, name_to_node, pair_to_edge))
    area, tree, name_to_node, pair_to_edge = max(quads)
    return tree, name_to_node, pair_to_edge


def get_response_content(fs):
    # read the edges of the directed graph
    whitelist = set(string.uppercase + string.lowercase + string.digits)
    pairs = []
    for line in fs.strlist:
        a, b = line.split()
        if set(a) - whitelist:
            raise ValueError('invalid name: ' + a)
        if set(b) - whitelist:
            raise ValueError('invalid name: ' + b)
        pairs.append((a, b))
    # make the graph in convenient xml form
    tree, name_to_node, pair_to_edge = make_graph(pairs, fs.width_inches)
    # modify the node attributes of the svg and get the name id map
    name_to_node_id, name_pair_to_edge_id = modify_tree(
            tree, name_to_node, pair_to_edge)
    # get the topological sort of the node ids
    names = name_to_node.keys()
    topo_names = reversed(graph.topo_sort(names, pairs))
    topo_ids = [name_to_node_id[name] for name in topo_names]
    # get the map from node name to child edge ids
    name_to_edge_ids = dict((name, []) for name in names)
    for name_pair, edge_id in name_pair_to_edge_id.items():
        source_name, sink_name = name_pair
        name_to_edge_ids[source_name].append(edge_id)
    # get the map from edge id to sink id
    edge_id_to_sink_id = {}
    for name_pair, edge_id in name_pair_to_edge_id.items():
        source_name, sink_name = name_pair
        edge_id_to_sink_id[edge_id] = name_to_node_id[sink_name]
    # wrap the embedded svg into some html
    svg_str = etree.tostring(tree)
    out = StringIO()
    print >> out, '<html>'
    print >> out
    print >> out, '<head>'
    print >> out, "<script type='text/javascript'>"
    print >> out, g_script_init
    print >> out, 'var topo_sorted_node_ids = ['
    for node_id in topo_ids:
        print >> out, "'%s'," % node_id
    print >> out, '];'
    print >> out, 'var edge_id_to_sink_id = {'
    for edge_id, sink_id in edge_id_to_sink_id.items():
        print >> out, "%s: '%s'," % (edge_id, sink_id)
    print >> out, '};'
    print >> out, 'var node_id_to_edge_ids = {'
    for source_name, edge_ids in name_to_edge_ids.items():
        source_id = name_to_node_id[source_name]
        s = ', '.join("'%s'" % x for x in edge_ids)
        print >> out, '%s: [%s],' % (source_id, s)
    print >> out, '};'
    print >> out, 'var node_id_to_name = {'
    for name, node_id in name_to_node_id.items():
        print >> out, "%s: '%s'," % (node_id, name)
    print >> out, '};'
    print >> out, g_script_guts
    print >> out, '</script>'
    print >> out, '</head>'
    print >> out
    print >> out, '<body>'
    print >> out
    print >> out, '<fieldset>'
    print >> out, '<legend>selections</legend>'
    print >> out, svg_str
    print >> out, '</fieldset>'
    print >> out
    print >> out, g_html_form
    print >> out
    print >> out, '</body>'
    print >> out, '</html>'
    return out.getvalue()

