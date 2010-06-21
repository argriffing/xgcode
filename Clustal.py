"""Minimally deal with clustalw formatted multiple sequence alignment.

This module assumes that the clustalw file fits easily in memory.
"""

import unittest
from StringIO import StringIO

import Util


class ClustalError(Exception):
    pass


def get_headers_and_sequences(raw_lines):
    """
    @param raw_lines: raw lines of a clustalw multiple sequence alignment
    @return: a (headers, sequences) pair
    """
    taxa = None
    sequences_of_chunks = None
    paragraphs = Util.get_paragraphs(raw_lines)
    for paragraph in paragraphs:
        # identify and skip the header paragraph
        first_line = paragraph[0]
        if len(paragraph) == 1 and first_line.startswith('CLUSTAL'):
            continue
        # break each line into a taxon and a sequence chunk
        taxon_chunk_pairs = [line.split() for line in paragraph]
        for taxon_chunk_pair in taxon_chunk_pairs:
            if len(taxon_chunk_pair) != 2:
                msg_a = 'each sequence line '
                msg_b = 'should be two whitespace separated strings'
                raise ClustalError(msg_a + msg_b)
        # validate the paragraph (taxon, chunk) pairs
        paragraph_taxa, chunks = zip(*taxon_chunk_pairs)
        if len(set(len(chunk) for chunk in chunks)) > 1:
            msg_a = 'within each paragraph each sequence chunk '
            msg_b = 'should be the same length'
            raise ClustalError(msg_a + msg_b)
        if taxa:
            if len(taxa) != len(paragraph_taxa):
                msg = 'each paragraph should have the same number of taxa'
                raise ClustalError(msg)
            if taxa != paragraph_taxa:
                msg_a = 'each paragraph should have the same taxa '
                msg_b = 'in the same order'
                raise ClustalError(msg_a + msg_b)
        else:
            taxa = paragraph_taxa
        # append the chunks to the sequences of chunks
        if not sequences_of_chunks:
            sequences_of_chunks = [[] for taxon in taxa]
        for sequence_of_chunks, chunk in zip(sequences_of_chunks, chunks):
            sequence_of_chunks.append(chunk)
    # for each taxon merge the sequence of chunks into a single string
    sequences = [''.join(x) for x in sequences_of_chunks]
    return taxa, sequences

# This is how someone decided that the sample taxa are supposed to be grouped
g_sample_cluster_labels = (
        1,1,1,1,2,3,3,3,3,3,3,3,4,1,1,1,1,1,1,1,1,3,3,3,3,5,5,5,3)

# Sample data from http://pbil.univ-lyon1.fr/mva/pco.php
g_sample_data = """
CLUSTAL W (1.7) multiple sequence alignment


S.enterica          -------------------------------------------------------MTDNR
S.typhimurium       -------------------------------------------------------MTDNR
E.coli              -------------------------------------------------------MTNYR
Y.pestis            -------------------------------------------------------MLNGR
D.radiodurans       -------------------------------------MDSRSPSLPDDRPDPPEQHLDAR
B.subtilis          ------------------------------------------------MMNKDITAQSPR
C.acetobutylicum_1  --------------------------------------------MDKLSIREDNYNIKHK
M.bovis             -------------------------------------------------LSHLCS-RGSK
M.tuberculosis      ------------------------------------------MAGEFRLLSHLCS-RGSK
M.leprae            ---------------------------------------------MFVNLTCLSEGNGSR
L.monocytogenes     ------------------------MKKDKTERTKQSW---------RKAQ----NAPSLS
L.innocua           ------------------------MKKEKTERTKQSW---------RKAQ----NAPSLS
Anabaena            ------------------------------------M---------TPPE----NRPSLP
A.tumefasciens      ---------MFFMRMICNCINSISMKAQKMDKPVFGW---------RRNG----DDLSLS
P.aeruginosa_1      ------------------------MSAKDTPAPQA-------------------AEP--G
X.fastidiosa        MSCGIANAIRWDVSLLIHRSMSSSMSSVQPTSPVS-------------------DSPSLG
R.loti              ------------------------MSDAEATAPRSSW---------RFAGRDEDDQPSLR
P.aeruginosa_2      ---------------------------MKYSLPTTAT-----------------APFCPS
P.aeruginosa_3      ------------------------------------------------------------
B.solanacearum      -------------------------MFRLPALPRIPT-----------------APFCPS
B.cepacia           ------------------------------------------------------------
S.aureus_1          ----------------------------------MNN---------KRHSTNE--QLSLD
S.aureus_2          ----------------------------------MNN---------KRHSTNE--QLSLD
L.brevis            ------------------MKEGIDMRESVAEEPKHKF---------FQSTEGE--NKSLD
L.lactis            ------------------------MND--KEHQRHRL---------IYHTNG----KSLE
S.cerevisiae_1      ------------------------------------------------------------
S.cerevisiae_2      --------------MTSQEYEPIQWSDESQTNNDSVN-------DAYADVNTTHESRRRT
S.cerevisiae_3      --MVNVGPSHAAVAVDASEARKRNISEEVFELRDKKDSTVVIEGEAPVRTFTSSSSNHER
C.acetobutylicum_2  ----------------------------------------------------------MR


S.enterica          VENSSG-----RAARKLRLALMGPAFIAAIGYIDPGNFATNIQAGASFGYQLLWVVVWAN
S.typhimurium       VENSSG-----RAARKLRLALMGPAFIAAIGYIDPGNFATNIQAGASFGYQLLWVVVWAN
E.coli              VESSSG-----RAARKMRLALMGPAFIAAIGYIDPGNFATNIQAGASFGYQLLWVVVWAN
Y.pestis            AVDTSR-----RPLRKIKLSLMGPAFIAAIAYIDPGNFATNIQAGATFGYTLLWVVVWAN
D.radiodurans       AGATLRGTAGPRGVRRI-LPFLGPAVIASIAYMDPGNFATNIEGGARYGYSLLWVILAAN
B.subtilis          SKAVQDALDGKIRGFRGLLPFLGPAFIAAIAYIDPGNFATNISAGSKYGYMLLWVILFSN
C.acetobutylicum_1  TISLNRASS-VTRNLKGLLKFLGPAFVVSVAYIDPGNFATNISGGSSFNYNLIWVILWSN
M.bovis             VGELAQDT---RTSLKTSWYLLGPAFVAAIAYVDPGNVAANVSSGAQFGYLLLWVIVAAN
M.tuberculosis      VGELAQDT---RTSLKTSWYLLGPAFVAAIAYVDPGNVAANVSSGAQFGYLLLWVIVAAN
M.leprae            GKALAEST---QAELKTSWYLLGPAFVAAIAYVDPGNVAANVSSGAQFGYLQLWVVVVAN
L.monocytogenes     EVNNSVAIPKNAKFFRKLFAFMGPGALIAVGYVDPGNWATSIAGGSEFGYTLLSVILISN
L.innocua           EVNNSVAIPKNAKFFRKLFAFMGPGALIAVGYVDPGNWATSIAGGSEFGYTLLSVILISN
Anabaena            EVHRSIRVPNSNSFWRKMLAYAGPGYLVSVGYIDPGNWATDIAGGSKFGYTLLTVILLSN
A.tumefasciens      DVHSSIRIKPDASTFRRAMAFFGPGYLVAVGYMDPGNWATSLAGGSKFGYTLLAVALVSN
P.aeruginosa_1      RAQP-VEVPGGGHWWRRLLAFAGPGYLVAVGYMDPGNWATDVAGGAQLGYLLLSVILLSS
X.fastidiosa        EMHASVAVSRRGHWGFRLLAFLGPGYMVSVGYMDPGNWATGLAGGSRFGYMLLSVILLSN
R.loti              EVNSTIAVPSSGVWFRRLFAFMGPGYMVSVGYMDPGNWATDLAGGAQFGYTLLFVIMLSN
P.aeruginosa_2      AVSHSVAVPADASPLRKLALFVGPGLLVSVGYMDPGNWATAIEAGSRFGYALLFVVVLAS
P.aeruginosa_3      -----MAVPADASPLRKLALFVGPGLLVSVGYMDPGNWATAIEAGSRFGYALLFVVVLAS
B.solanacearum      EVRGCVAIPPDLPLWKKLLRFAGPGLLVSVGYMDPGNWATDIEAGSRYGYSLLFVVLLSS
B.cepacia           -------------------------------------WATDIQAGSQFGYSLLWVVAFSS
S.aureus_1          EINNTIKFDHRSSNKQKFLSFLGPGLLVAVGYMDPGNWITSMQGGAQYGYTLLFVILISS
S.aureus_2          EINNTIKFDHRSSNKQKFLSFLGPGLLVAVGYMDPGNWITSMQGGAQYGYTLLFVILISS
L.brevis            EVNGSVKVPKNAGFWRTLFAYTGPGILIAVGYMDPGNWITSIAGGAQFKYSLLSVILISS
L.lactis            EINGTVEVPKNVGFFKMLLTYSGPGALVAVGYMDPGNWSTSITGGQNFQYLLISVILMSS
S.cerevisiae_1      -------MRSYMQILQKFAKFIGPGILVSVAYMDPGNYATSVSGGAQYKYTLLFSIFISN
S.cerevisiae_2      TLQPNSTSQSMIGTLRKYARFIGPGLMVSVSYMDPGNYSTAVAAGSAHRYKLLFSVLVSN
S.cerevisiae_3      EDTYVSKRQVMRDIFAKYLKFIGPGLMVSVAYIDPGNYSTAVDAGASNQFSLLCIILLSN
C.acetobutylicum_2  KINIFKGKHNPQRTAIDFFKYVGPGLLVTVGFIDPGNWASNVSAGSDYGYSLLWIVTLST


S.enterica          LMAMLIQILSAKLGIATGKNLAEQIRDHYPRPVVWFYWVQAEIIAMATDLAEFIGAAIGF
S.typhimurium       LMAMLIQILSAKLGIATGKNLAEQIRDHYPRPVVWFYWVQAEIIAMATDLAEFIGAAIGF
E.coli              LMAMLIQILSAKLGIATGKNLAEQIRDHYPRPVVWFYWVQAEIIAMATDLAEFIGAAIGF
Y.pestis            VMAMLVQLLSAKLGIATGKNLAEHIRDRFPRPVVWAYWVQAEIIVMATDLAEFIGAAIGF
D.radiodurans       LMAMVIQNLSANLGIASGRNLPELIRERWPRPLVWFYWIQAELVAMATDLAEFLGAALAI
B.subtilis          IMALLIQSLSAKLGIATGKNLPEVAREEFPKPVSIGLWIQGELVIIATDLAEFIGAALGL
C.acetobutylicum_1  LMAIFLQTMSAKLGIATGCSLPEMCAKVFSKRANWIFWIVGELGAMATDLAEFIGGTLGL
M.bovis             VMAALVQYLSAKLGLVTGRSLPEAIGKRMGRPARLAYWAQAEIVAMATDVAEVIGGAIAL
M.tuberculosis      VMAALVQYLSAKLGLVTGRSLPEAIGKRMGRPARLAYWAQAEIVAMATDVAEVIGGAIAL
M.leprae            VLAGLVQYLSAKLGLVTGQSLPQAISKQMSHPFRLGFWLQAELVAMATDVAEIVGGAIAF
L.monocytogenes     ILAVLLQSLASKLGIVTGRDLAQASSDHFSKPFGFVLWILAELAIIATDIAEVIGSAIAL
L.innocua           ILAVLLQSLASKLGIVTGRDLAQASSDHFSKPFGFVLWILAELAIIATDIAEVIGSAIAL
Anabaena            LMAILLQSLCVRLGVATGRDLAQACRDYFSPKVSFCLWVLCEIAIAACDLAELLGSAIAL
A.tumefasciens      IMAIVLQSLCARLAIASGRDLAQACRDAYPKPVAMVLWLLAEIAIIATDIAEVIGTAIGL
P.aeruginosa_1      LMAMLLQALSARLGIASGLDLAQACRERYSPSTCRLLWLACETAIIACDLAEVIGTAIAL
X.fastidiosa        VMAIVLQALAARLGIASDMDLAQACRARYSRGTTLALWVVCELAIIACDLAEVIGTAIAL
R.loti              LMAILLQALAARLGIATGRDLAQACRAYYPRPVNFVLWIACELAIIACDLAEVIGTAIAL
P.aeruginosa_2      FSGMLLQSLCSRLGIATGRDLAQLSRERYRPGVARGQWLLAELSIVATDLAEVLGAALAF
P.aeruginosa_3      FSGMLLQSLCSRLGIATGRDLAQLSRERYRPGVARGQWLLAELSIVATDLAEVLGAALAF
B.solanacearum      LAAMVLQCLSARLGIVTGKDLARLSRERYRPGAVRVQWLLAELSIVACDLAEVLGCALAF
B.cepacia           LAAIFLQMLAARLGLVAGKDLAQASYERYGRFGRVVQWITAEVSIIACDIAEVLGCALAF
S.aureus_1          LSAMLLQSMTVRLGIATGMDLAQMTRHYLSRPIAIIFWIIAELAIIATDIAEVIGSAIAL
S.aureus_2          LSAMLLQSMTVRLGIATGMDLAQMTRHYLSRPIAIIFWIIAELAIIATDIAEVIGSAIAL
L.brevis            LIAMLLQSMAARLGIVTGKDLAQLTRERTSKTMGIILWLITESAIMATDVAEIIGSGIAI
L.lactis            LIAMLLQYMSAKLGIVSQMDLAQAIRARTSKSLGIVLWILTELAIMATDIAEVIGAAIAL
S.cerevisiae_1      IFAVLLQCLCVKLGTITGYDLAENCRHNLPKKLNYTLYLFAEVAIIATDLAEVVGTAIAL
S.cerevisiae_2      FMAAFWQYLCARLGAVTGLDLAQNCKKHLPFGLNITLYILAEMAIIATDLAEVVGTAISL
S.cerevisiae_3      FIAIFLQCLCIKLGSVTGLDLSRACREYLPRWLNWTLYFFAECAVIATDIAEVIGTAIAL
C.acetobutylicum_2  IMLIILQHNAAHLGIVTGDCLSEAVSKKFKPWLTNIILYSAIGAAISTAMAELLGGAIAL


S.enterica          KLILGVSLLQGAVLTGIATFLILMLQ-RRGQ------KPLEKVIGGLLLFVAAAYIVELF
S.typhimurium       KLILGVSLLQGAVLTGIATFLILMLQ-RRGQ------KPLEKVIGGLLLFVAAAYIVELF
E.coli              KLILGVSLLQGAVLTGIATFLILMLQ-RRGQ------KPLEKVIGGLLLFVAAAYIVELI
Y.pestis            KLLFGVTLLQGAVLTGIATFLILMLQ-NRGQ------KPLELVIGGLLLFVAAAYIVELI
D.radiodurans       QLLTGLPMFWGAVVTGVVTFWLLNLQ-KRGT------RPLELAVGAFVLMIGVAYLVQVV
B.subtilis          YLLFGIPMLEASIIAAIGSFAILELQ-RRGY------RSLEAGIAGMLFVVVIAFALQTF
C.acetobutylicum_1  YLLFRIPMIYAGLLTGVLTFIIVYME-KYGQ------KMVETIIAALIAVICVAYTIELF
M.bovis             RIMFNLPLPIGGIITGVVSLLLLTIQDRRGQ------RLFERVITALLLVIAIGFTASFF
M.tuberculosis      RIMFNLPLPIGGIITGVVSLLLLTIQDRRGQ------RLFERVITALLLVIAIGFTASFF
M.leprae            HILFRVSLLLGGVITGTVSLLLLMVKDRRGQ------LLFERVITGLLFVIVVGFTSSFF
L.monocytogenes     NLLFGIPLIWGVCITALDIFLVLFLQ-HKGF------RYIEVIVITLMVTILVCFGAEMV
L.innocua           NLLFGIPLIWGVCITALDIFLVLFLQ-HKGF------RYIEVIVITLMVTILVCFGAEMV
Anabaena            QLLFVIPLIWGVCITALDVLVLLFLQ-HKGF------RYTEALVIMLVATVGICFTAEIL
A.tumefasciens      NLIFGIPLELGVLITALDVFLILYLQ-KLGF------RWVEALVITLLGVIAVCFAIQLA
P.aeruginosa_1      KLLFGLPLAWGALLCVGDALLVLVLI-GRGQ------RPLEAFVVALLTLIFACFAVQLL
X.fastidiosa        NLLLGVPIIWGVVITAVDVVLVLLLM-HRGF------RALEAFVIALLLVIFGCFVVQIV
R.loti              KLLFGIPLIGGAILTALDAFLVLLLM-NKGF------RYLEAFVIALLIIIFSCFAIQIF
P.aeruginosa_2      HLLLGVSITTGVVLTAFDTLIVLALQ-GANF------RRLEAIVLGLIATIGACFFVELV
P.aeruginosa_3      HLLLGVSITTGVVLTAFDTLIVLALQ-GANF------RRLEAIVLGLIATIGACFFVELV
B.solanacearum      HLLLGVPILGGVALTALDTLIVLGLK-GKNF------RQLEAIVLGLILTIGLCYFVELA
B.cepacia           KLLLGVPLAWGVVLTALDTVIVLGLQ-GKGF------RQIEAIVLGLIATMAFCFVAQVA
S.aureus_1          NLLFNIPLIVGALITVLDVFLLLFIM-KYGF------RKIEAIVGTFIFTVLFIFIFEVY
S.aureus_2          NLLFNIPLIVGALITVLDVFLLLFIM-KYGF------RKIEAIVGTLIFTVLFIFIFEVY
L.brevis            KLLFNIPLVVGILITTADVLILLLLM-KLGF------RKIEAIVATLVAVILLVFTYEVF
L.lactis            YLLFNIPLILAVFITVLDVFLLLLLT-RVGF------RKIEALVICLILVILVIFAYQVA
S.cerevisiae_1      QILFKIPLTWGVLLTVLDVLVILMFY-TPNGQSLKKVRVFEFGVGILVIGTCICFVLELF
S.cerevisiae_2      NILFHIPLALGVILTVVDVLIVLLAY-KPNG-SMKGIRIFEAFVSLLVVLTVVCFTVELF
S.cerevisiae_3      NILIKVPLPAGVAITVVDVFLIMFTY-KPGASSIRFIRIFECFVAVLVVGVCICFAIELA
C.acetobutylicum_2  NMLFRLPIKIGTLIMLSLVIWMLFSS---SY------KKLERWIIGFVSIIGLSFIFELL


S.enterica          FSQPD--MAQLGKGMVIP------ALPNP---EAVFLAAGVLGATIMPHVIYLHSSLTQ-
S.typhimurium       FSQPD--MAQLGKGMVIP------ALPNP---EAVFLAAGVLGATIMPHVIYLHSSLTQ-
E.coli              FSQPN--LAQLGKGMVIP------SLPTS---EAVFLAAGVLGATIMPHVIYLHSSLTQ-
Y.pestis            FSQPD--IAALGRGMLIP------NLPDG---NAVFLAAGVLGATIMPHVIYLHSALTQ-
D.radiodurans       LARPD--LAAVGAG-FVP------RLQGP---GSAYLAVGIIGATVMPHVIYLHSALTQG
B.subtilis          FAKPD--AVSVMKGLFVP------AFHGT---DSVLLAAGILGATVMPHAIYLHSALTQR
C.acetobutylicum_1  LARPA--WTQVGMHTLIP------SLPNG---EAVLIAVGMLGATVMPHVIYLHSELVQH
M.bovis             VVTPP--PNAVLGG-LAP------RFQGT---ESVLLAAAIMGATVMPHAVYLHSGLARD
M.tuberculosis      VVTPP--PNAVLGG-LAP------RFQGT---ESVLLAAAIMGATVMPHAVYLHSGLARD
M.leprae            VATPS--PEDMVNG-LLP------RFQGT---ESVLLAAAIIGATVMPHAVYLHSGLALD
L.monocytogenes     MSHPD--MQAIAKG-FIP------QSEIVTNPAMLYIALGILGATVMPHNLYLHSSIVQT
L.innocua           MSHPD--MQAIAKG-FIP------QSEIVTNPAMLYIALGILGATVMPHNLYLHSSIVQT
Anabaena            FSRPD--MGGILLG-YLP------KKEILQNPEMLYIAIGILGATVMPHNLYLHSSIVQT
A.tumefasciens      LADPD--WGQVILG-FAP------TTEIVTNPDMLYLALGILGATVMPHNLYLHSGIVQT
P.aeruginosa_1      LSRPE--LGEVLQG-FLP------SPRVLSDPAALYLAIGIVGATVMPHNLYLHSSLVQS
X.fastidiosa        LAAPP--LQEVLGG-FVP------RWQVVADPQALYLAIGIVGATVMPHNLYLHSSIVQT
R.loti              VAAPP--AGTILHSMFVP------SSEIVTNPAMLYIAIGIIGATVMPHNLYLHSSIVQT
P.aeruginosa_2      LIGPY--WPDVAAG-LRP------SWDTLSSQEPLYLAIGILGATVMPHNLYLHSSVVQT
P.aeruginosa_3      LIGPY--WPDVAAG-LRP------SWDTLSSQEPLYLAIGILGATVMPHNLYLHSSVVQT
B.solanacearum      LIRPH--WPSVAGA-LVP------SWQALSAREPLYLAIGILGATVMPHNLYLHSSIVQT
B.cepacia           ITPPD--WHAVVGG-LVP------GDPGHDRKDAIVLALGIVGATI--------------
S.aureus_1          ISSPQ--LNAVLNG-FIP------HSEIITNNGILYIALGIIGATIMPHNLYLHSSIVQS
S.aureus_2          ISSPQ--LNAVLNG-FIP------HSEIITNNGILYIALGIIGATIMPHNLYLHSSIVQS
L.brevis            LAGPQ--LDQMFAG-YMP------TKDIVTNKSMLYLALGIVGATVMPHDLYLGSSISQT
L.lactis            LSNPD--WKGIIEG-FIPNSRTFASSPTVAGMSPLTGALGIIGATVMPHNLYLHSSISQS
S.cerevisiae_1      KVSIP-DKAELFKG-FLP------SNIIFKEQQALYISLGILGATVMPHSLYLGSSIVKP
S.cerevisiae_2      YAKLG-PAKEIFSG-FLP------SKAVF-EGDGLYLSLAILGATVMPHSLYLGSGVVQP
S.cerevisiae_3      YIPKSTSVKQVFRG-FVP------SAQMF-DHNGIYTAISILGATVMPHSLFLGSALVQP
C.acetobutylicum_2  LVKVN--WQVAAVSFVSP------SFPKN----SMPIIMSVLGAVVMPHNLFLHSEIIQS


S.enterica          -----------------------------------HLHGGTRQQR---YSATKWDVAIAM
S.typhimurium       -----------------------------------HLHGGTRQQR---YSATKWDVAIAM
E.coli              -----------------------------------HLHGGSRQQR---YSATKWDVAIAM
Y.pestis            -----------------------------------TGGEESKTER---YASTKFDVAIAM
D.radiodurans       -----------------------------------RIQTDTTEEKRRLVRLNRVDVIAAM
B.subtilis          -----------------------------------RVVGKTDAERKKIFRFEFIDILIAM
C.acetobutylicum_1  -----------------------------------RNTNSSDKEKLHHLKMEKIDILIAM
M.bovis             R----------------------------------HGHPDPGPQRRRLLRVTRWDVGLAM
M.tuberculosis      R----------------------------------HGHPDPGPQRRRLLRVTRWDVGLAM
M.leprae            R----------------------------------HGHPHAGRSRRRLLRVTRLDVILAM
L.monocytogenes     R-----------------------------Q-----YAR-TKEGKKEAIRFSFIDSTFSL
L.innocua           R-----------------------------Q-----YAR-TKEGKREAIRFSFIDSTFSL
Anabaena            R-----------------------------D-----WQP-TTEKRWEAIKFGTIDSTFAL
A.tumefasciens      R-----------------------------E-----IGP-TIAEKREALKFATLDSTIAL
P.aeruginosa_1      R-----------------------------A-----YPR-SLAGKRKALRWAVADSSLAL
X.fastidiosa        R-----------------------------A-----YPR-TPVGRRSALRWAVADSTLAL
R.loti              R-----------------------------A-----YER-TEKGKRDAIKWATTDSTIAL
P.aeruginosa_2      R-----------------------------V-----SGD-DAASKRSAIRFSRLDTIGSL
P.aeruginosa_3      R-----------------------------V-----SGD-DAASKRSAIRFSRLDTIGSL
B.solanacearum      R-----------------------------V-----VSE-TEPARREAVGLSRLDTIVSL
B.cepacia           ------------------------------------------------------------
S.aureus_1          R-----------------------------T-----YSRHNNEEKAQAIKFATIDSNIQL
S.aureus_2          R-----------------------------T-----YSRHNNEEKAQAIKFATIDSNIQL
L.brevis            R-----------------------------A-----VDRHDRQDVAKAIKFTTIDSNLQL
L.lactis            R-----------------------------K-----INHQDKSDVARAVRFSTWDSNIQL
S.cerevisiae_1      RLHDYDLKK-----------------------YGKVNARPSLSAIKYTLNYAYAELIISL
S.cerevisiae_2      RLREYDIKNGHYLPDAND------------MDNNHDNYRPSYEAISETLHFTITELLISL
S.cerevisiae_3      RLLDYDVKHGNYTVSEEQDKVKKSKSTEEIMEEKYFNYRPTNAAIKYCMKYSMVELSITL
C.acetobutylicum_2  R----------------------------------QWNVEKEDVIKRQLRYEFLDTIISM


S.enterica          T-IAGFVNLAMMATAAAAFHFSGHTG-IADLDQAYLTLEPLLSHAAAT------VFGLSL
S.typhimurium       T-IAGFVNLAMMATAAAAFHFSGHTG-IADLDQAYLTLEPLLSHAAAT------VFGLSL
E.coli              T-IAGFVNLAMMATAAAAFHFSGHTG-VADLDEAYLTLQPLLSHAAAT------VFGLSL
Y.pestis            T-IAGFVNLAMMATAAAAFHFNGYEN-IAEIEEAYITLQPLLGNAAAT------VFGLSL
D.radiodurans       G-LAGLINMSMLAVAAATFHGKNVEN-AGDLTTAYQTLTPLLGPAASV------LFAVAL
B.subtilis          L-IAGAINASMLIVAAALFFKNGLF--VEDLDVAFQQFGHLVSPMSAA------LFGIGL
C.acetobutylicum_1  N-IAFVVNAAMVIVSAAVFFKHGIK--VSTIEEAHRSLQPLLGNLSSG------AFGIAL
M.bovis             L-IAGGVNAAMLLVAALNMRGRGDT---ASIEGAYHAVHDTLGATIAV------LFAVGL
M.tuberculosis      L-IAGGVNAAMLLVAALNMRGRGDT---ASIEGAYHAVHDTLGATIAV------LFAVGL
M.leprae            T-IAGIVNTAMLLVAAINLQHHQVT---AYIEGTYTAIQDTLGATIAM------LFAIGL
L.monocytogenes     T-IALLINASILILAAAAFYTTGQ-HNVAGIEDAYKLLNPTLGSSIAS-----TVFAVAL
L.innocua           T-IALLINASILILAAAAFYTTGQ-HNVAGIEDAYKLLNPTLGSSIAS-----TVFAVAL
Anabaena            S-LALFINSAILIVSAATFHFSGN-QNVAEIQDAYKLLSPLLGVSAAS-----AIFGIAL
A.tumefasciens      M-FALLINASILILAAATFNKTGQ-TNVAELGEAHSLLAPLLGLAIAP-----TLFGVAL
P.aeruginosa_1      T-LALLVNAAILIVAASVFHRNGH-TEVVDIEQAHALLSPLLGLELAS-----LLFAVAL
X.fastidiosa        M-LALFINASILILAAAVFHAQHH-FDVEEIEQAYQLLAPVLGVGVAA-----TLFATAL
R.loti              M-LALFVNAAILIVSAVAFHNTGH-QDVAEIDQAFELLSPLLGLGIAS-----ILFAVAL
P.aeruginosa_2      S-LALLVNAAILILAAAAFHGSGH-TEVVEIQDAYHLLDPLVGGALAS-----FLFGFAL
P.aeruginosa_3      S-LALLVNAAILILAAAAFHGSGH-TEVVEIQDAYHLLDPLVGGALAS-----FLFGFAL
B.solanacearum      S-LALLVNGAILVLAAAAFHANGH-QDVADIQDAHRLLEPIVGTAVAG-----VLFGIAL
B.cepacia           ------------------------------------------------------------
S.aureus_1          S-IAFVVNCLLLVLGASLFFNSNA-DDLGGFYDLYHALKTEPVLGATMGAIMSTLFAVAL
S.aureus_2          S-IAFVVNCLLLVLGASLFFNSNA-DDLGGFYDLYHALKTEPVLGATMGAIMSTLFAVAL
L.brevis            T-IAFIVNSLLLILGAALFFGTNS-T-VGRFVDLFNSLNNSHIVGAIASPMLSMLFAVAL
L.lactis            T-VAFVVNSLLLIMGVAVFK-TG----V---IDDPSFFNHSFVDVD--------------
S.cerevisiae_1      FLIATFVNSAILIVAGATLSGQPE-AEDADLLSIYKLLVHYISPAAGL------IFALAM
S.cerevisiae_2      FTVALFVNCAILIVSGATLYGSTQNAEEADLFSIYNLLCSTLSKGAGT------VFVLAL
S.cerevisiae_3      FTLALFVNCAILVVAGSTLYNSPE-ADGADLFTIHELLSRNLAPAAGT------IFMLAL
C.acetobutylicum_2  I-VGFAINSAMVLLAASTFFKMHVS--VSELSQAQALLKPLLGGSASV------VFALAL


S.enterica          VAAGLSSTVVGTLAGQVVMQGFVRFHIPLWVRRTIT----MLP-SFIVILMG---LDPTR
S.typhimurium       VAAGLSSTVVGTLAGQVVMQGFVRFHIPLWVRRTIT----MLP-SFIVILMG---LDPTR
E.coli              VAAGLSSTVVGTLAGQVVMQGFIRFHIPLWVRRTVT----MLP-SFIVILMG---LDPTR
Y.pestis            IAAGLSSTVVGTLAGQVVMQGFVRFYIPMWVRRIVT----MLP-SFIVILAG---MDATQ
D.radiodurans       LASGLSSSAVGTMAGDVIMQGFMGFHIPLWLRRLIT----MLP-AFIVILLG---MDPSS
B.subtilis          LVAGLSSSSVGTLSGDVIMQGFINYRIPLYVRRFIT----IIP-PILIIASG---VNPTT
C.acetobutylicum_1  LASGLSSSAVGTMAGQTIMKGFVNLSIPINLRRIIT----MLP-ALIIIALG---INPMR
M.bovis             LASGLASSSVGAYAGAMIMQGLLHWSVPMLVRRLIT----LGP-ALAILTLG---FDPTR
M.tuberculosis      LASGLASSSVGAYAGAMIMQGLLHWSVPMLVRRLIT----LGP-ALAILTLG---FDPTR
M.leprae            LASSLASASVGAYAGALIMQGLLQRSIPMLIRRLIT----LCP-AIAILALG---FDPTR
L.monocytogenes     LASGQNSTLTGTLAGQIVMEGFLNIRLKPVVRRLLTRVLAIVPAVIITALYG--ANGINE
L.innocua           LASGQNSTLTGTLAGQIVMEGFLNIRLKPVVRRLLTRVLAIVPAVIITALYG--ANGINE
Anabaena            LASGQSSTLTATLAGQIVMEGFLQFRLPSWLRRLITRLLAIIPALITIILFG--ENSTSS
A.tumefasciens      LCCGINSTVTATLAGQIVMEGFLKMRLAPWLRRLITRAIAIVPAAGVTIFYG--DSGTGQ
P.aeruginosa_1      LASGLNSTVTTTLAGQIVMEGFLRLRLAPWARRLLTRGVAVLPVLLVTLLYG--EDGTAR
X.fastidiosa        LASGINSTVTATLAGQIVMEGFLRLRLRPWLRRVLTRGLAIVPVIVVVALYG--EQGTGR
R.loti              LASGLNSTVTATLAGQIIMEGFLRLRIPNWARRLLTRGLAIVPVVVVTALYG--EKGTGQ
P.aeruginosa_2      LAAGQSSTFTGTIAGQVVMEGFLRAKIPCWQRRLITRGLALVPALIGVLWLG--EAAVGK
P.aeruginosa_3      LAAGQSSTFTGTIAGQVVMEGFLRAKIPCWQRRLITRGLALVPALIGVLWLG--EAAVGK
B.solanacearum      LAAGQSSTFTGTIAGQILMEGFLELRIPCWQRRLATRALALIPAFVGVAMLG--DHAIGR
B.cepacia           ------------------------------------------------------------
S.aureus_1          LASGQNSTITGTLAGQIVMEGFLRLHIPNWLRRLITRSLAVIPVIVCLIIFKGNAAKIEQ
S.aureus_2          LASGQNSTITGTLAGQIVMEGFLRLHIPNWLRRLITRSLAVIPVIVCLIIFKGNAAKIEQ
L.brevis            LSSGQSSTITGTLSGQIIMEGFIHLRMPLWAQRLLTRLLSVTPVLIFAIYYHGNEAKIEN
L.lactis            ------------------------------------------------------------
S.cerevisiae_1      LCSGQSAGIICTLAGQIVSEGFLQWSLPPWATRLCTRLIAIVPCLFVTLTMG--EKGISD
S.cerevisiae_2      LFSGQSAGIVCTLSGQMVSEGFLNWTVSPALRRSATRAVAITPCLILVLVAG--RSGLSG
S.cerevisiae_3      LLSGQSAGVVCTMSGQIVSEGHINWKLQPWQRRLATRCISIIPCLVISICIG--REALSK
C.acetobutylicum_2  LFSGIASTVTAGMAGGSIFAGIYKEPYDIEDSHTQIG-VIITMVLAAVIIFF--IKDPFK


S.enterica          ILVMSQVLLSFGIALALVPLLIFTSNATLMGELVN-------------------------
S.typhimurium       ILVMSQVLLSFGIALALVPLLIFTSNATLMGELVN-------------------------
E.coli              ILVMSQVLLSFGIALALVPLLIFTSDSKLMGDLVN-------------------------
Y.pestis            ILVMSQVLLSFGIALALVPLLVFTGNKELMGELVD-------------------------
D.radiodurans       VLILSQVILCFGVPFALVPLLLFTARRDVMGALVT-------------------------
B.subtilis          ALVLSQVVLSFGIAFALIPLIMFTSNKRIMGSLIN-------------------------
C.acetobutylicum_1  VLVLSQVALSFILPFPIIQMLLIAGRKDLMGILVN-------------------------
M.bovis             TLVLSQVVLSFGIPFAVLPLVKLTGSPAVMGGDTN-------------------------
M.tuberculosis      TLVLSQVVLSFGIPFAVLPLVKLTGSPAVMGGDTN-------------------------
M.leprae            ALVLSQIVLSFGIPFAVLPLVKLTNNRGLMGNDTN-------------------------
L.monocytogenes     LLIFSQVILSMQLSFAVIPLVMFTSDKQKMGEFV--------------------------
L.innocua           LLIFSQVILSMQLSFAVIPLVMFTSDKQKMGEFV--------------------------
Anabaena            LIVLSQVILSLQLPFAVIPLVMFTSNRRLMGEFV--------------------------
A.tumefasciens      LLILTQVVLSLQLSFAVFPLVMFTSDKAKMGELR--------------------------
P.aeruginosa_1      LLIFSQVILSMQLPLAVIPLLQFVSDRRLMGPLA--------------------------
X.fastidiosa        LLLLSQVILSMQLPFAVIPLLRCVADRKVMGALV--------------------------
R.loti              LLVFSQVILSMQLPFAVVPLVQFVSDKKKMGNLA--------------------------
P.aeruginosa_2      LLVLSQVVLSLQLPFALWPLIRFSSDRGLMGEFV--------------------------
P.aeruginosa_3      LLVLSQVVLSLQLPFALWPLIRFSSDRGLMGEFV--------------------------
B.solanacearum      LLVISQVVLGFQLPFAMFPLIRMTGDRALMGTFA--------------------------
B.cepacia           ------------------------------------------------------------
S.aureus_1          LLVFSQVFLSIALPFCLIPLQLATSNKDLMGPFY--------------------------
S.aureus_2          LLVFSQVFLSIALPFCLIPLQLATSNKDLMGPFY--------------------------
L.brevis            LLTLSQVFLSVALPFAIVPLVKFTSSKELMGEFV--------------------------
L.lactis            ------------------------------------------------------------
S.cerevisiae_1      ILNFSQVVLSLILPIVSAPLIYFTANRKLMVVHD-ENGVVRAPADVNAIADETTP-----
S.cerevisiae_2      ALNASQVVLSLLLPFVSAPLLYFTSSKKIMRVQLNRTKELSRTTDKKPVADRTEDDETIE
S.cerevisiae_3      ALNASQVVLSIVLPFLVAPLIFFTCKKSIMKTEITVDHTEEDSHNHQNNNDRSAG-----
C.acetobutylicum_2  GLIYSQMILSIQLPITIFTQIYLTSSKKVMGKFSN-------------------------


S.enterica          ---------------------------TRRVKQIGWIIVVLVVALNIWLLVGTVMGLS--
S.typhimurium       ---------------------------TRRVKQVGWIIVVLVVALNIWLLVGTVMGLS--
E.coli              ---------------------------SKRVKQTGWVIVVLVVALNIWLLVGTALGL---
Y.pestis            ---------------------------TKTTQILGKLVVLIVVGLNAYLLISLL------
D.radiodurans       ---------------------------RRSFTVIGWVIAVIIIALNGYLLWELLGG----
B.subtilis          ---------------------------AKWITVVSWLIAVLIVALNVFLIVDTFR-----
C.acetobutylicum_1  ---------------------------KKFTKIVGFIIATMIILLNIILLYLTFTGQT--
M.bovis             ---------------------------HRATTWVGWVVAVMVSLLNVMLI----------
M.tuberculosis      ---------------------------HRATTWVGWVVAVMVSLLNVMLIYLTVTG----
M.leprae            ---------------------------HPATTVLGWAVAILVSLLNVVLIYLTVTS----
L.monocytogenes     ------------N--------------PTWLKIISWAVAIFIAVLNIYLLFYTLTSL---
L.innocua           ------------N--------------SPWLKIVSWSVAIFIAFLNIYLLFYTLTSL---
Anabaena            ------------N--------------PLWLKSLAWLVAIVIVGLNAWLLLQSLWGWLLQ
A.tumefasciens      ------------S--------------PLWLSAIAWLIAVVIAALNVKLLMDFMG-----
P.aeruginosa_1      ------------I--------------GAGTRWLAWAVALAIVGLNLQLLADFAFG----
X.fastidiosa        ------------A--------------PRWLMVVAWLIAGVIVVLNVKLLGDYAVHLMVG
R.loti              ------------I--------------PRGVAALAWVVAAIILVLNFKLLYDTLFG--VG
P.aeruginosa_2      ------------N--------------PRWVSALAWSLFGLISAANLTLLYFWFG-----
P.aeruginosa_3      ------------N--------------PRWVSALAWSLFGLISAANLTLLYFWFG-----
B.solanacearum      ------------N--------------GRLTSAVAWCLFAVISVANLWLVWQVLAG----
B.cepacia           ------------------------------------------------------------
S.aureus_1          ------------N--------------KTWVNIISWTLIIILSILNVYLIVQTFQELQS-
S.aureus_2          ------------N--------------KTWVNIISWTLIIILSILNVYLIVQTFQELQS-
L.brevis            ------------N--------------KAWVKYSAWVATVVLVSLNIYLILQTVGVIG--
L.lactis            ------------------------------------------------------------
S.cerevisiae_1      --------------LNSKHSKIVDFTNSRLLTYTSVFVWALIGSLNCYLVISYLLG-ADI
S.cerevisiae_2      LEEMGIGSSSQERSLVSPAPEYKDMSNGMIVTVLAIIVWLIISGLNFYMLLGFTTG-KEV
S.cerevisiae_3      -SVIEQDGSSGMEIENGKDVKIVYMANNWIITVIAIIVWLFLSLLNVYAIVQLGMSHGDI
C.acetobutylicum_2  ---------------------------SLLDRVLLGIIAVIVTALNIALLVSYF------


S.enterica          ---
S.typhimurium       ---
E.coli              ---
Y.pestis            ---
D.radiodurans       ---
B.subtilis          ---
C.acetobutylicum_1  ---
M.bovis             ---
M.tuberculosis      ---
M.leprae            ---
L.monocytogenes     ---
L.innocua           ---
Anabaena            VPS
A.tumefasciens      ---
P.aeruginosa_1      ---
X.fastidiosa        VSD
R.loti              ---
P.aeruginosa_2      ---
P.aeruginosa_3      ---
B.solanacearum      ---
B.cepacia           ---
S.aureus_1          ---
S.aureus_2          ---
L.brevis            ---
L.lactis            ---
S.cerevisiae_1      HF-
S.cerevisiae_2      HL-
S.cerevisiae_3      S--
C.acetobutylicum_2  ---

"""


class TestClustal(unittest.TestCase):

    def test_clustal(self):
        """
        Parse the sample multiple alignment.
        """
        # parse the raw lines of text
        raw_data_lines = StringIO(g_sample_data).readlines()
        headers, sequences = get_headers_and_sequences(raw_data_lines)
        # check for basic inconsistencies
        self.assertEqual(len(headers), len(sequences))
        self.assertEqual(len(headers), len(g_sample_cluster_labels))
        # Compare the header and sequence of the third item
        # to the expected values.
        expected_taxon = 'E.coli'
        expected_sequence = (
                '-------------------------------------------------------MTNYR'
                'VESSSG-----RAARKMRLALMGPAFIAAIGYIDPGNFATNIQAGASFGYQLLWVVVWAN'
                'LMAMLIQILSAKLGIATGKNLAEQIRDHYPRPVVWFYWVQAEIIAMATDLAEFIGAAIGF'
                'KLILGVSLLQGAVLTGIATFLILMLQ-RRGQ------KPLEKVIGGLLLFVAAAYIVELI'
                'FSQPN--LAQLGKGMVIP------SLPTS---EAVFLAAGVLGATIMPHVIYLHSSLTQ-'
                '-----------------------------------HLHGGSRQQR---YSATKWDVAIAM'
                'T-IAGFVNLAMMATAAAAFHFSGHTG-VADLDEAYLTLQPLLSHAAAT------VFGLSL'
                'VAAGLSSTVVGTLAGQVVMQGFIRFHIPLWVRRTVT----MLP-SFIVILMG---LDPTR'
                'ILVMSQVLLSFGIALALVPLLIFTSDSKLMGDLVN-------------------------'
                '---------------------------SKRVKQTGWVIVVLVVALNIWLLVGTALGL---'
                '---')
        self.assertEqual(expected_taxon, headers[2])
        self.assertEqual(expected_sequence, sequences[2])

if __name__ == '__main__':
    unittest.main()
