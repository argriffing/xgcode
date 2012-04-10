"""
This module is related to the BEAST primate tutorial.
"""

from StringIO import StringIO

import const
import mcmc
import Fasta
import RUtil

g_fasta_string = const.read('20120405a').strip()
g_ntax = 12
g_nchar = 898

g_headers = (
        'sequence.length',
        'midpoint',
        'mean.low', 'mean.mean', 'mean.high',
        'var.low', 'var.mean', 'var.high',
        'cov.low' ,'cov.mean', 'cov.high')

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

class BeastLogFileError(Exception): pass

def get_xml_post_alignment(nsamples):
    return """
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
	<mcmc id="mcmc" chainLength="%d" autoOptimize="true">
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
    """ % nsamples

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

def get_xml_string(start_pos, stop_pos, nsamples, log_path):
    """
    @param start_pos: start position within the hardcoded alignment
    @param stop_pos: stop position within the hardcoded alignment
    @param nsamples: run the mcmc chain for this many samples
    @param log_path: tell beast to put its posterior sample log here
    @return: multiline xml string
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
    print >> out, get_xml_post_alignment(nsamples)
    print >> out, get_log_xml(log_path)
    return out.getvalue().rstrip()

def get_html(values_name_pairs):
    """
    Web based only.
    """
    out = StringIO()
    #
    #print >> out, 'statistic:', name
    #print >> out, 'mean:', corr_info.mean
    #print >> out, 'standard error of mean:', corr_info.stdErrorOfMean
    #print >> out, 'auto correlation time (ACT):', corr_info.ACT
    #print >> out, 'standard deviation of ACT:', corr_info.stdErrOfACT
    #print >> out, 'effective sample size (ESS):', corr_info.ESS
    #print >> out, 'posterior density interval (0.95): [%f, %f]' % hpd
    #print >> out
    print >> out, '<html>'
    # write the html head
    print >> out, '<head>'
    print >> out, '<script'
    print >> out, '  type="text/javascript"'
    print >> out, '  src="https://www.google.com/jsapi">'
    print >> out, '</script>'
    print >> out, '<script type="text/javascript">'
    print >> out, "  google.load('visualization', '1', {packages:['table']});"
    print >> out, "  google.setOnLoadCallback(drawTable);"
    print >> out, "  function drawTable() {"
    print >> out, "    var data = new google.visualization.DataTable();"
    # add columns
    print >> out, "    data.addColumn('string', 'description');"
    print >> out, "    data.addColumn('number', '95% HPD low');"
    print >> out, "    data.addColumn('number', 'mean');"
    print >> out, "    data.addColumn('number', '95% HPD high');"
    print >> out, "    data.addColumn('number', 'ACT');"
    print >> out, "    data.addColumn('number', 'ESS');"
    # add rows
    print >> out, "    data.addRows(3);"
    # add entries
    for i, (values, name) in enumerate(values_name_pairs):
        corr_info = mcmc.Correlation()
        corr_info.analyze(values)
        hpd_low, hpd_high = mcmc.get_hpd_interval(0.95, values)
        print >> out, "    data.setCell(%d, 0, '%s');" % (i, name)
        print >> out, "    data.setCell(%d, 1, %f);" % (i, hpd_low)
        print >> out, "    data.setCell(%d, 2, %f);" % (i, corr_info.mean)
        print >> out, "    data.setCell(%d, 3, %f);" % (i, hpd_high)
        print >> out, "    data.setCell(%d, 4, %f);" % (i, corr_info.ACT)
        print >> out, "    data.setCell(%d, 5, %f);" % (i, corr_info.ESS)
    print >> out, "    var table = new google.visualization.Table("
    print >> out, "      document.getElementById('table_div'));"
    print >> out, "    table.draw(data, {showRowNumber: false});"
    print >> out, "  }"
    print >> out, "</script>"
    print >> out, '</head>'
    # write the html body
    print >> out, '<body><div id="table_div"></div></body>'
    # end the html
    print >> out, '</html>'
    # return the html string
    return out.getvalue().rstrip()

def read_log(log_loc, nsamples_expected):
    """
    @param log_loc: path to the log file
    @param nsamples_expected: expected number of mcmc posterior samples
    @return: means, variations, covariances
    """
    with open(log_loc) as fin:
        lines = [line.strip() for line in fin.readlines()]
        iines = [line for line in line if line]
    # check the number of non-whitespace lines
    expected = nsamples_expected + 3 + 1
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
    # skip the first three lines
    # skip the initial state
    # skip ten percent of the remaining states
    nburnin = nsamples_expected / 10
    for line in lines[3 + 1 + nburnin:]:
        s1, s2, s3, s4 = line.split()
        state = int(s1)
        means.append(float(s2))
        variations.append(float(s3))
        covariances.append(float(s4))
    return means, variations, covariances

def get_R_tick_cmd(axis, positions):
    """
    @param axis: 1 for x, 2 for y
    @param positions: a sequence of positions
    @return: a single line R command to draw the ticks
    """
    s = 'c(' + ', '.join(str(x) for x in positions) + ')'
    return RUtil.mk_call_str('axis', axis, at=s)

def get_ggplot2_x_tick_cmd(positions):
    s = 'c(' + ', '.join(str(x) for x in positions) + ')'
    return RUtil.mk_call_str('scale_x_discrete', breaks=s)

def get_ggplot2_legend_cmd():
    s_labels = "c('57', '114', '228', '456')"
    return RUtil.mk_call_str('scale_colour_discrete',
            labels=s_labels)

def get_ggplot2_scripts(nsamples, sequence_lengths, midpoints):
    scripts = []
    # get the plot for the mean
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'ggplot', 'my.table',
            RUtil.mk_call_str(
                'aes',
                x='midpoint',
                y='mean.mean')), '+'
    print >> out, RUtil.mk_call_str(
            'geom_errorbar',
            RUtil.mk_call_str(
                'aes',
                ymin='mean.low',
                ymax='mean.high',
                colour='factor(sequence.length)'),
            width='20'), '+'
    print >> out, "opts(title='mcmc chain length %d') +" % nsamples
    print >> out, "geom_point() + xlab('midpoint') + ylab('mean of rates') +"
    print >> out, "scale_color_discrete('length') +"
    print >> out, get_ggplot2_x_tick_cmd(midpoints)
    scripts.append(out.getvalue().rstrip())
    # get the plot for the coefficient of variation
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'ggplot', 'my.table',
            RUtil.mk_call_str(
                'aes',
                x='midpoint',
                y='var.mean')), '+'
    print >> out, RUtil.mk_call_str(
            'geom_errorbar',
            RUtil.mk_call_str(
                'aes',
                ymin='var.low',
                ymax='var.high',
                colour='factor(sequence.length)'),
            width='20'), '+'
    print >> out, "geom_point() + xlab('midpoint') +"
    print >> out, "ylab('coefficient of variation of rates') +"
    print >> out, "scale_color_discrete('length') +"
    print >> out, get_ggplot2_x_tick_cmd(midpoints)
    scripts.append(out.getvalue().rstrip())
    # get the plot for the correlation
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'ggplot', 'my.table',
            RUtil.mk_call_str(
                'aes',
                x='midpoint',
                y='cov.mean')), '+'
    print >> out, RUtil.mk_call_str(
            'geom_errorbar',
            RUtil.mk_call_str(
                'aes',
                ymin='cov.low',
                ymax='cov.high',
                colour='factor(sequence.length)'),
            width='20'), '+'
    print >> out, "geom_point() + xlab('midpoint') +"
    print >> out, "ylab('parent child correlation of rates') +"
    print >> out, "scale_color_discrete('length') +"
    print >> out, get_ggplot2_x_tick_cmd(midpoints)
    scripts.append(out.getvalue().rstrip())
    return scripts

def get_table_string_and_scripts(nsamples):
    """
    Command-line only.
    """
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_pos, stop_pos in g_start_stop_pairs:
        sequence_length = stop_pos - start_pos + 1
        means, variations, covs = get_value_lists(
                start_pos, stop_pos, nsamples)
        midpoint = (start_pos + stop_pos) / 2.0
        row = [sequence_length, midpoint]
        for values in means, variations, covs:
            corr_info = mcmc.Correlation()
            corr_info.analyze(values)
            hpd_low, hpd_high = mcmc.get_hpd_interval(0.95, values)
            row.extend([hpd_low, corr_info.mean, hpd_high])
        data_arr.append(row)
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table string
    table_string = RUtil.get_table_string(data_arr, g_headers)
    # get the scripts
    scripts = get_ggplot2_scripts(nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, scripts

def get_table_string_and_scripts_from_logs(
        start_stop_pairs, log_paths, nsamples):
    """
    This is for analysis of remote execution.
    """
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_stop_pair, log_path in zip(
            start_stop_pairs, log_paths):
        start_pos, stop_pos = start_stop_pair
        sequence_length = stop_pos - start_pos + 1
        means, variations, covs = read_log(log_path, nsamples)
        midpoint = (start_pos + stop_pos) / 2.0
        row = [sequence_length, midpoint]
        for values in means, variations, covs:
            corr_info = mcmc.Correlation()
            corr_info.analyze(values)
            hpd_low, hpd_high = mcmc.get_hpd_interval(0.95, values)
            row.extend([hpd_low, corr_info.mean, hpd_high])
        data_arr.append(row)
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table string
    table_string = RUtil.get_table_string(data_arr, g_headers)
    # get the scripts
    scripts = get_ggplot2_scripts(nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, scripts

