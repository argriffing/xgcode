"""
Run BEAST on a sub-alignment.

The full alignment and the BEAST analysis parameters are hardcoded.
It uses the twelve-taxon primate analysis from the BEAST tutorial.
Running this should take a few seconds
before showing moments of posterior distributions
of a few rate statistics.
"""

from StringIO import StringIO
import math
import random
import os
import subprocess

import numpy as np

import Form
import FormOut
import const
import Util
import Fasta

g_fasta_string = const.read('20120405a').strip()
g_ntax = 12
g_nchar = 898
g_beast_root = os.path.expanduser('~/svn-repos/beast-mcmc-read-only')

#TODO use lxml


g_xml_pre_alignment = """
<?xml version="1.0" standalone="yes"?>
<beast>
	<!-- The list of taxa analyse (can also include dates/ages). -->
	<taxa id="taxa">
		<taxon id="Tarsius_syrichta"/>
		<taxon id="Lemur_catta"/>
		<taxon id="Homo_sapiens"/>
		<taxon id="Pan"/>
		<taxon id="Gorilla"/>
		<taxon id="Pongo"/>
		<taxon id="Hylobates"/>
		<taxon id="Macaca_fuscata"/>
		<taxon id="M_mulatta"/>
		<taxon id="M_fascicularis"/>
		<taxon id="M_sylvanus"/>
		<taxon id="Saimiri_sciureus"/>
	</taxa>
	<taxa id="Human-Chimp">
		<taxon idref="Homo_sapiens"/>
		<taxon idref="Pan"/>
	</taxa>
	<taxa id="ingroup">
		<taxon idref="Gorilla"/>
		<taxon idref="Homo_sapiens"/>
		<taxon idref="Hylobates"/>
		<taxon idref="M_fascicularis"/>
		<taxon idref="M_mulatta"/>
		<taxon idref="M_sylvanus"/>
		<taxon idref="Macaca_fuscata"/>
		<taxon idref="Pan"/>
		<taxon idref="Pongo"/>
		<taxon idref="Saimiri_sciureus"/>
		<taxon idref="Tarsius_syrichta"/>
	</taxa>
	<taxa id="HomiCerco">
		<taxon idref="Gorilla"/>
		<taxon idref="Homo_sapiens"/>
		<taxon idref="Hylobates"/>
		<taxon idref="M_fascicularis"/>
		<taxon idref="M_mulatta"/>
		<taxon idref="M_sylvanus"/>
		<taxon idref="Macaca_fuscata"/>
		<taxon idref="Pan"/>
		<taxon idref="Pongo"/>
	</taxa>
""".strip()

g_xml_post_alignment = """
	<yuleModel id="yule" units="substitutions">
		<birthRate>
			<parameter id="yule.birthRate" value="1.0"
              lower="0.0" upper="Infinity"/>
		</birthRate>
	</yuleModel>
	<constantSize id="initialDemo" units="substitutions">
		<populationSize>
			<parameter id="initialDemo.popSize" value="100.0"/>
		</populationSize>
	</constantSize>
	<!-- Generate a random starting tree under the coalescent process -->
	<coalescentTree id="startingTree">
		<constrainedTaxa>
			<taxa idref="taxa"/>
			<tmrca monophyletic="false">
				<taxa idref="Human-Chimp"/>
			</tmrca>
			<tmrca monophyletic="true">
				<taxa idref="ingroup"/>
			</tmrca>
			<tmrca monophyletic="false">
				<taxa idref="HomiCerco"/>
			</tmrca>
		</constrainedTaxa>
		<constantSize idref="initialDemo"/>
	</coalescentTree>
	<!-- Generate a tree model -->
	<treeModel id="treeModel">
		<coalescentTree idref="startingTree"/>
		<rootHeight>
			<parameter id="treeModel.rootHeight"/>
		</rootHeight>
		<nodeHeights internalNodes="true">
			<parameter id="treeModel.internalNodeHeights"/>
		</nodeHeights>
		<nodeHeights internalNodes="true" rootNode="true">
			<parameter id="treeModel.allInternalNodeHeights"/>
		</nodeHeights>
	</treeModel>
	<!-- Generate a speciation likelihood for Yule or Birth Death -->
	<speciationLikelihood id="speciation">
		<model>
			<yuleModel idref="yule"/>
		</model>
		<speciesTree>
			<treeModel idref="treeModel"/>
		</speciesTree>
	</speciationLikelihood>
	<!--
      The uncorrelated relaxed clock
      (Drummond, Ho, Phillips & Rambaut, 2006)
    -->
	<discretizedBranchRates id="branchRates">
		<treeModel idref="treeModel"/>
		<distribution>
			<logNormalDistributionModel meanInRealSpace="true">
				<mean>
					<parameter id="ucld.mean" value="0.033"
                      lower="0.0" upper="1.0"/>
				</mean>
				<stdev>
					<parameter id="ucld.stdev" value="0.3333333333333333"
                      lower="0.0" upper="Infinity"/>
				</stdev>
			</logNormalDistributionModel>
		</distribution>
		<rateCategories>
			<parameter id="branchRates.categories" dimension="22"/>
		</rateCategories>
	</discretizedBranchRates>
	<rateStatistic id="meanRate" name="meanRate"
      mode="mean" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
	<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation"
      mode="coefficientOfVariation" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
	<rateCovarianceStatistic id="covariance" name="covariance">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateCovarianceStatistic>
	<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985) -->
	<HKYModel id="firsthalf.hky">
		<frequencies>
			<frequencyModel dataType="nucleotide">
				<frequencies>
					<parameter id="firsthalf.frequencies"
                      value="0.25 0.25 0.25 0.25"/>
				</frequencies>
			</frequencyModel>
		</frequencies>
		<kappa>
			<parameter id="firsthalf.kappa"
              value="2.0" lower="0.0" upper="Infinity"/>
		</kappa>
	</HKYModel>
	<!-- site model -->
	<siteModel id="firsthalf.siteModel">
		<substitutionModel>
			<HKYModel idref="firsthalf.hky"/>
		</substitutionModel>
		<gammaShape gammaCategories="4">
			<parameter id="firsthalf.alpha"
              value="0.5" lower="0.0" upper="1000.0"/>
		</gammaShape>
	</siteModel>
	<treeLikelihood id="firsthalf.treeLikelihood" useAmbiguities="false">
		<patterns idref="firsthalf.patterns"/>
		<treeModel idref="treeModel"/>
		<siteModel idref="firsthalf.siteModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</treeLikelihood>
	<tmrcaStatistic id="tmrca(Human-Chimp)" includeStem="false">
		<mrca>
			<taxa idref="Human-Chimp"/>
		</mrca>
		<treeModel idref="treeModel"/>
	</tmrcaStatistic>
	<tmrcaStatistic id="tmrca(ingroup)" includeStem="false">
		<mrca>
			<taxa idref="ingroup"/>
		</mrca>
		<treeModel idref="treeModel"/>
	</tmrcaStatistic>
	<monophylyStatistic id="monophyly(ingroup)">
		<mrca>
			<taxa idref="ingroup"/>
		</mrca>
		<treeModel idref="treeModel"/>
	</monophylyStatistic>
	<tmrcaStatistic id="tmrca(HomiCerco)" includeStem="false">
		<mrca>
			<taxa idref="HomiCerco"/>
		</mrca>
		<treeModel idref="treeModel"/>
	</tmrcaStatistic>
	<!-- Define operators -->
	<operators id="operators">
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="firsthalf.kappa"/>
		</scaleOperator>
		<deltaExchange delta="0.01" weight="0.1">
			<parameter idref="firsthalf.frequencies"/>
		</deltaExchange>
		<scaleOperator scaleFactor="0.75" weight="0.1">
			<parameter idref="firsthalf.alpha"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="ucld.mean"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="ucld.stdev"/>
		</scaleOperator>
		<subtreeSlide size="0.9" gaussian="true" weight="15">
			<treeModel idref="treeModel"/>
		</subtreeSlide>
		<narrowExchange weight="15">
			<treeModel idref="treeModel"/>
		</narrowExchange>
		<wideExchange weight="3">
			<treeModel idref="treeModel"/>
		</wideExchange>
		<wilsonBalding weight="3">
			<treeModel idref="treeModel"/>
		</wilsonBalding>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="treeModel.rootHeight"/>
		</scaleOperator>
		<uniformOperator weight="30">
			<parameter idref="treeModel.internalNodeHeights"/>
		</uniformOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="yule.birthRate"/>
		</scaleOperator>
		<upDownOperator scaleFactor="0.75" weight="3">
			<up>
				<parameter idref="ucld.mean"/>
			</up>
			<down>
				<parameter idref="treeModel.allInternalNodeHeights"/>
			</down>
		</upDownOperator>
		<swapOperator size="1" weight="10" autoOptimize="false">
			<parameter idref="branchRates.categories"/>
		</swapOperator>
		<randomWalkIntegerOperator windowSize="1" weight="10">
			<parameter idref="branchRates.categories"/>
		</randomWalkIntegerOperator>
		<uniformIntegerOperator weight="10">
			<parameter idref="branchRates.categories"/>
		</uniformIntegerOperator>
	</operators>
	<!-- Define MCMC -->
	<mcmc id="mcmc" chainLength="8000" autoOptimize="true">
		<posterior id="posterior">
			<prior id="prior">
				<booleanLikelihood>
					<monophylyStatistic idref="monophyly(ingroup)"/>
				</booleanLikelihood>
				<normalPrior mean="6.0" stdev="0.5">
					<statistic idref="tmrca(Human-Chimp)"/>
				</normalPrior>
                <normalPrior mean="24.0" stdev="0.5">
					<statistic idref="tmrca(HomiCerco)"/>
				</normalPrior>
				<logNormalPrior mean="1.0" stdev="1.25"
                  offset="0.0" meanInRealSpace="false">
					<parameter idref="firsthalf.kappa"/>
				</logNormalPrior>
				<exponentialPrior mean="0.3333333333333333" offset="0.0">
					<parameter idref="ucld.stdev"/>
				</exponentialPrior>
				<speciationLikelihood idref="speciation"/>
			</prior>
			<likelihood id="likelihood">
				<treeLikelihood idref="firsthalf.treeLikelihood"/>
			</likelihood>
		</posterior>
		<operators idref="operators"/>
		<!-- write log to screen -->
        <!--
		<log id="screenLog" logEvery="10000">
			<column label="Posterior" dp="4" width="12">
				<posterior idref="posterior"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			<column label="rootHeight" sf="6" width="12">
				<parameter idref="treeModel.rootHeight"/>
			</column>
			<column label="ucld.mean" sf="6" width="12">
				<parameter idref="ucld.mean"/>
			</column>
		</log>
        -->
""".rstrip()

class BeastLogFileError(Exception): pass

def get_log_xml(log_loc):
    s = """
		<!-- write log to file -->
		<log id="fileLog" logEvery="1" fileName="%s" overwrite="false">
            <!--
			<posterior idref="posterior"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="treeModel.rootHeight"/>
			<tmrcaStatistic idref="tmrca(Human-Chimp)"/>
			<tmrcaStatistic idref="tmrca(ingroup)"/>
			<tmrcaStatistic idref="tmrca(HomiCerco)"/>
			<parameter idref="yule.birthRate"/>
			<parameter idref="firsthalf.kappa"/>
			<parameter idref="firsthalf.frequencies"/>
			<parameter idref="firsthalf.alpha"/>
			<parameter idref="ucld.mean"/>
			<parameter idref="ucld.stdev"/>
			<treeLikelihood idref="firsthalf.treeLikelihood"/>
			<speciationLikelihood idref="speciation"/>
            -->
			<rateStatistic idref="meanRate"/>
			<rateStatistic idref="coefficientOfVariation"/>
			<rateCovarianceStatistic idref="covariance"/>
		</log>
		<!-- write tree log to file -->
        <!--
		<logTree id="treeFileLog" logEvery="200" nexusFormat="true"
          fileName="primates.trees" sortTranslationTable="true">
			<treeModel idref="treeModel"/>
			<discretizedBranchRates idref="branchRates"/>
			<posterior idref="posterior"/>
		</logTree>
        -->
            </mcmc>
            <!--
            <report>
                <property name="timer">
                    <mcmc idref="mcmc"/>
                </property>
            </report>
            -->
        </beast>
        """ % log_loc
    return s

def make_xml(start_pos, stop_pos):
    """
    @return: location of xml file, location of log file
    """
    out = StringIO()
    print >> out, g_xml_pre_alignment
    print >> out, """
        <!-- The sequence alignment (each sequence refers to a taxon above). -->
        <alignment id="alignment" dataType="nucleotide">
    """
    lines = g_fasta_string.splitlines()
    for header, seq in Fasta.gen_header_sequence_pairs(lines):
        print >> out, '<sequence>'
        print >> out, '<taxon idref="%s"/>' % header
        print >> out, seq
        print >> out, '</sequence>'
    print >> out, '</alignment>'
    print >> out, """
        <patterns id="firsthalf.patterns" from="%d" to="%d">
            <alignment idref="alignment"/>
        </patterns>
    """ % (start_pos, stop_pos)
    print >> out, g_xml_post_alignment
    log_loc = Util.get_tmp_filename(prefix='beast', suffix='.log')
    print >> out, get_log_xml(log_loc)
    xml_loc = Util.create_tmp_file(
            out.getvalue(), prefix='beast', suffix='.xml')
    return xml_loc, log_loc

def run_beast(xml_loc):
    args = (
            'java',
            '-jar',
            os.path.join(g_beast_root, 'build', 'dist', 'beast.jar'),
            xml_loc,
            )
    subprocess.call(args)
    

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('start',
                'sub-sequence start position (1-%d)' % g_nchar,
                1, low=1, high=g_nchar),
            Form.Integer('stop',
                'sub-sequence stop position (1-%d)' % g_nchar,
                g_nchar, low=1, high=g_nchar)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # init the response and get the user variables
    start_pos = fs.start
    stop_pos = fs.stop
    if stop_pos < start_pos:
        raise ValueError('the stop pos must be after the start pos')
    out = StringIO()
    # create the xml describing the analysis
    xml_loc, log_loc = make_xml(start_pos, stop_pos)
    # run beast
    run_beast(xml_loc)
    # read the log file
    with open(log_loc) as fin:
        lines = [line.strip() for line in fin.readlines()]
        iines = [line for line in line if line]
    # check the number of non-whitespace lines
    expected = 8004
    observed = len(lines)
    if expected != observed:
        msg= 'expected %d lines but observed %d' % (expected, observed)
        raise BeastLogFileError(msg)
    # check the first line
    expected = '# BEAST'
    if not lines[0].startswith(expected):
        msg = 'expected the first line to start with ' + expected
        raise BeastLogFileError(msg)
    # check the second line
    expected = '# Generated'
    if not lines[1].startswith(expected):
        msg = 'expected the second line to start with ' + expected
        raise BeastLogFileError(msg)
    # check the third line
    values = lines[2].split()
    if len(values) != 4:
        msg = 'expected the third line to have four column labels'
        raise BeastLogFileError(msg)
    if values != ['state', 'meanRate', 'coefficientOfVariation', 'covariance']:
        msg = 'unexpected column labels on the third line'
        raise BeastLogFileError(msg)
    # read the rest of the lines
    means = []
    variations = []
    covariances = []
    for line in lines[3:]:
        s1, s2, s3, s4 = line.split()
        state = int(s1)
        means.append(float(s2))
        variations.append(float(s3))
        covariances.append(float(s4))
    print >> out, 'summaries of posterior distributions of rate statistics'
    print >> out
    print >> out, 'statistic: mean rate among branches'
    print >> out, 'posterior mean:', np.mean(means)
    print >> out, 'posterior stdev:', np.std(means)
    print >> out
    print >> out, 'statistic: coefficient of variation of rates among branches'
    print >> out, 'posterior mean:', np.mean(variations)
    print >> out, 'posterior stdev:', np.std(variations)
    print >> out
    print >> out, 'statistic: correlation of parent and child branch rates'
    print >> out, 'posterior mean:', np.mean(covariances)
    print >> out, 'posterior stdev:', np.std(covariances)
    print >> out
    # return the response
    return out.getvalue()

