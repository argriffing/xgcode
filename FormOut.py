"""
This is how the output format of a snippet is specified.

Some metaprogramming is done for compatibility with Mobyle.
In Mobyle, all script output is to a file or to stdout.
"""

import Form

try:
    from lxml import etree
except ImportError as e:
    pass

#FIXME galaxy has trouble with varying output types
#so it thinks that all images are pngs

def _get_filename_metaprogram(fmt, interpolants):
    if interpolants:
        if len(interpolants) == 1:
            return '"%s" %% %s' % (fmt, interpolants[0])
        elif len(interpolants) > 1:
            return '"%s" %% (%s)' % (fmt, ', '.join(interpolants))
    else:
        return '"%s"' % fmt

class FormOut(object):

    def __init__(self,
            name_base='default', name_ext='out', name_interpolants=None):
        if name_ext:
            self.filename_format_string = name_base + '.' + name_ext
        else:
            self.filename_format_string = name_base
        self.filename_interpolants = name_interpolants

    def get_filename_metaprogram(self):
        """
        This is ugly but necessary metaprogramming for Mobyle compatibility.
        Mobyle wants python code embedded in XML.
        @return: a string that is a python expression
        """
        return _get_filename_metaprogram(
            self.filename_format_string, self.filename_interpolants)

    def get_filename(self, fs):
        """
        @param fs: fieldstorage
        """
        if self.filename_interpolants:
            rhs = tuple(getattr(fs, x) for x in self.filename_interpolants)
            return self.filename_format_string % rhs
        else:
            return self.filename_format_string

    def get_contenttype(self, fs):
        return 'text/plain'

    def get_response_headers(self, fs):
        response_headers = []
        response_headers.append(('Content-Type', self.get_contenttype(fs)))
        if hasattr(fs, 'contentdisposition'):
            filename = self.get_filename(fs)
            disposition = '%s; filename=%s' % (fs.contentdisposition, filename)
            response_headers.append(('Content-Disposition', disposition))
        return response_headers

    def get_galaxy_format(self):
        return 'txt'

    def get_mobyle_class(self):
        """
        Mobyle does not need to understand this name.
        Or maybe it does.
        For now, change everything to Text.
        """
        #FIXME why doesn't mobyle like custom names?
        #return self.__class__.__name__
        return 'Text'

    def get_mobyle_superclass(self):
        """
        This must be a name like 'Text' or 'Integer' that Mobyle understands.
        """
        return 'Text'

    def _add_mob_xml_datatypes(self, parent):
        etree.SubElement(parent, 'class').text = self.get_mobyle_class()
        superclass = self.get_mobyle_superclass()
        if superclass:
            etree.SubElement(parent, 'superclass').text = superclass

    def add_mob_xml(self, parent, short_name):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param short_name: the short name of the snippet module
        """
        desc = short_name + '_out'
        parameter = etree.SubElement(parent, 'parameter', isout='1')
        etree.SubElement(parameter, 'name').text = desc
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = desc
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        self._add_mob_xml_datatypes(datatype)
        filenames = etree.SubElement(parameter, 'filenames')
        mycode = etree.SubElement(filenames, 'code', proglang='python')
        mycode.text = self.get_filename_metaprogram()


class Image(FormOut):

    def __init__(self, base_format_string='img', base_interpolants=[]):
        self.filename_format_string = base_format_string + '.%s'
        self.filename_interpolants = base_interpolants + ['imageformat']

    def get_contenttype(self, fs):
        return Form.g_imageformat_to_contenttype[fs.imageformat]

    def get_mobyle_class(self):
        return 'Picture'

    def get_mobyle_superclass(self):
        return 'Binary'

    def get_galaxy_format(self):
        return 'png'


class Png(FormOut):

    def __init__(self, fmt='img', ext='png', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

    def get_contenttype(self, fs):
        return 'image/png'

    def get_mobyle_class(self):
        return 'Picture'

    def get_mobyle_superclass(self):
        return 'Binary'

    def get_galaxy_format(self):
        return 'png'


class Tikz(FormOut):

    def __init__(self, base_format_string='tkz', base_interpolants=[]):
        self.filename_format_string = base_format_string + '.%s'
        self.filename_interpolants = base_interpolants + ['tikzformat']

    def get_contenttype(self, fs):
        return Form.g_tikzformat_to_contenttype[fs.tikzformat]

    def get_mobyle_class(self):
        return 'Picture'

    def get_mobyle_superclass(self):
        return 'Binary'

    def get_galaxy_format(self):
        return 'png'


class Latex(FormOut):

    def __init__(self, base_format_string='latex', base_interpolants=[]):
        self.filename_format_string = base_format_string + '.%s'
        self.filename_interpolants = base_interpolants + ['latexformat']

    def get_contenttype(self, fs):
        return Form.g_latexformat_to_contenttype[fs.latexformat]

    def get_mobyle_class(self):
        return 'Picture'

    def get_mobyle_superclass(self):
        return 'Binary'

    def get_galaxy_format(self):
        return 'png'




class ContextDependent(FormOut):
    def __init__(self, fmt='out', ext='unknown', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class Report(FormOut):
    def __init__(self, fmt='out', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class Hud(FormOut):
    def __init__(self, fmt='out', ext='hud', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class RTable(FormOut):
    def __init__(self, fmt='out', ext='table', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class EigenstratInd(FormOut):
    def __init__(self, fmt='out', ext='ind', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class EigenstratPheno(FormOut):
    def __init__(self, fmt='out', ext='pheno', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class EigenstratGeno(FormOut):
    def __init__(self, fmt='out', ext='geno', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class EigenstratSnp(FormOut):
    def __init__(self, fmt='out', ext='snp', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class Newick(FormOut):
    def __init__(self, fmt='out', ext='newick', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Tree'

class Fasta(FormOut):
    def __init__(self, fmt='out', ext='fasta', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'

class Phylip(FormOut):
    def __init__(self, fmt='out', ext='phy', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'

class TransitionMatrix(FormOut):
    def __init__(self, fmt='transition_matrix', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class NucleotideRateMatrix(FormOut):
    def __init__(self, fmt='nt_rate_matrix', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class CodonRateMatrix(FormOut):
    def __init__(self, fmt='codon_rate_matrix', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class Html(FormOut):
    def __init__(self, fmt='out', ext='html', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_contenttype(self, fs):
        return 'text/html'

class BeastXml(FormOut):
    def __init__(self, fmt='out', ext='xml', interpolants=None):
        FormOut.__init__(self, fmt, interpolants)
    def get_contenttype(self, fs):
        return 'text/xml'

class RateMatrix(FormOut):
    def __init__(self, fmt='rate_matrix', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class NucleotideFasta(FormOut):
    def __init__(self, fmt='nt_alignment', ext='fasta', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'

class CodonFasta(FormOut):
    def __init__(self, fmt='codon_alignment', ext='fasta', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'

class Nexus(FormOut):
    def __init__(self, fmt='out', ext='nex', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'

class Stockholm(FormOut):
    def __init__(self, fmt='out', ext='stockholm', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'

class Matrix(FormOut):
    def __init__(self, fmt='matrix', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class CorrelationMatrix(FormOut):
    def __init__(self, fmt='correlation_matrix', ext='txt', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class Csv(FormOut):
    def __init__(self, fmt='out', ext='csv', interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)

class Gff(FormOut):
    def __init__(self, fmt='out', ext='gff', interpolants=None):
        FormOut.__init__(self, fmt, interpolants)

class Alignment(FormOut):
    def __init__(self, fmt='alignment', ext=None, interpolants=None):
        FormOut.__init__(self, fmt, ext, interpolants)
    def get_mobyle_class(self):
        return 'Alignment'
