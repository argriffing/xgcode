"""
Make an example image with pydot.
"""

from StringIO import StringIO

import pydot

import Form
import FormOut
import Util

def get_form():
    """
    @return: the body of a form
    """
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Svg('dag')

def get_plain_graph():
    graph = pydot.Dot(graph_type='graph')
    for i in range(3):
        edge = pydot.Edge('king', 'lord%d' % i)
        graph.add_edge(edge)
        for j in range(2):
            vassal_num = i*2 + j
            edge = pydot.Edge('lord%d' % i, 'vassal%d' % vassal_num)
            graph.add_edge(edge)
    return graph

def get_fancy_graph():
    graph = pydot.Dot(
            graph_type='graph', size="4, 8",
            nodesep=str(0.25 + 0.234),
            ranksep=str(0.75 - 0.234))
    king = pydot.Node('omg', URL='"http://www.google.com"')
    graph.add_node(king)
    for i in range(3):
        edge = pydot.Edge(king, 'foo%d' % i)
        graph.add_edge(edge)
        for j in range(2):
            vassal_num = i*2 + j
            edge = pydot.Edge('foo%d' % i, 'bar%d' % vassal_num)
            graph.add_edge(edge)
    return graph

def get_response_content(fs):
    graph = get_fancy_graph()
    tmp_path = Util.create_tmp_file(data=None, prefix='tmp', suffix='.svg')
    #graph.write_svg(tmp_path, prog='neato')
    graph.write_svg(tmp_path)
    with open(tmp_path) as fin:
        svg_str = fin.read()
    return svg_str
