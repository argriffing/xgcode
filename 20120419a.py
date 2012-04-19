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
    for elem in tree.iter():
        parent = elem.getparent()
        print elem.tag, elem.text
        if elem.tag == '{http://www.w3.org/2000/svg}g':
            if elem.get('class') == 'node':
                elem.set('cursor', 'pointer')
                elem.set('onclick', 'dosomethingcool(evt);')
        elif elem.tag == '{http://www.w3.org/2000/svg}text':
            elem.set('pointer-events', 'none')
        elif elem.tag == '{http://www.w3.org/2000/svg}polygon':
            elem.set('fill', 'white')
        elif elem.tag == '{http://www.w3.org/2000/svg}title':
            if parent.get('class') == 'node':
                name_to_node_id[elem.text] = parent.get('id')
    return name_to_node_id

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
    # initialize the graph
    g = pydot.Dot(
            graph_type='graph',
            size='%s, 8' % fs.width_inches,
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
    # modify the node attributes of the svg and get the name id map
    name_to_node_id = modify_tree(tree, name_to_node, pair_to_edge)
    # get the topological sort of the node ids
    names = name_to_node.keys()
    topo_names = reversed(graph.topo_sort(names, pairs))
    topo_ids = [name_to_node_id[name] for name in topo_names]
    # wrap the embedded svg into some html
    svg_str = etree.tostring(tree)
    out = StringIO()
    print >> out, '<html>'
    print >> out
    print >> out, '<head>'
    print >> out, "<script type='text/javascript'>"
    print >> out, 'var topo_sorted_node_ids = ['
    for node_id in topo_ids:
        print >> out, node_id + ','
    print >> out, '];'
    print >> out, '</script>'
    print >> out, '</head>'
    print >> out
    print >> out, '<body>'
    print >> out
    print >> out, 'hello'
    print >> out, '<div>'
    print >> out, svg_str
    print >> out, '</div>'
    print >> out, 'goodbye'
    print >> out
    print >> out, '</body>'
    print >> out, '</html>'
    return out.getvalue()

