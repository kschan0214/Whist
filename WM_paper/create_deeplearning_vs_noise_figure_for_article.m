clear
close all


loss(1, :) = [0.17175695270796618, 0.13106820748373865, 0.11835515262807408, 0.11156319030125936, 0.10714281022921204, 0.1036791924238205, 0.10096599850803614, 0.09859991363311807, 0.09663867319747806, 0.09494291415934761, 0.09351901936779419, 0.09207456689948837, 0.09082911949480574, 0.08974624820798635, 0.08873015727351109, 0.08774292489513755, 0.08683838049322366, 0.08603243186449011, 0.08522582525884112, 0.08445392028987407, 0.08379305142660935, 0.08310892755289873, 0.08250990175952513, 0.0819320538242658, 0.0813456176618735, 0.08074446869144837, 0.08023616481820742, 0.07980699828267097, 0.0792837530101339, 0.07878733686109383, 0.078394211375465, 0.07796370585635304, 0.07752223065371315, 0.0771684679997464, 0.07673569291457534, 0.07636935281132658, 0.0760602826885879, 0.07568628142153223, 0.07530884752919277, 0.07499201961482564];
loss(2, :) = [0.20677172106561206, 0.1693283399320784, 0.15494955731460025, 0.1470076146750223, 0.1417489738657361, 0.1380937023809978, 0.13526646114417484, 0.13299154129255386, 0.13115769295692445, 0.12950186515876225, 0.12801974383308773, 0.12675026755560012, 0.12560755904402052, 0.12454861305668241, 0.12356827075765246, 0.12263541248298827, 0.1218379644655046, 0.12107872943707874, 0.12028790130615234, 0.11964690584682283, 0.119062885368438, 0.11842462663026083, 0.11781319834618341, 0.11728153826055072, 0.11677095705781664, 0.11628133637110392, 0.11583004493259248, 0.11534777188471386, 0.11493943161850884, 0.11452880527235212, 0.114103874145803, 0.11372490855171567, 0.11334171476477668, 0.11306217373723075, 0.11270675939548583, 0.11233563977025804, 0.11201668843201229, 0.11170862214281445, 0.11142289702154341, 0.11114932047185444];
loss(3, :) = [0.21641174866699037, 0.1809471282175609, 0.16963995158785866, 0.1626189982527778, 0.15805955619925544, 0.15491174162683033, 0.15245558829534622, 0.15044705957912263, 0.14865989773727598, 0.14716844973450616, 0.14586337960788182, 0.14467911795207433, 0.1436771626165935, 0.14266655084178562, 0.14180775527613504, 0.14100592116741906, 0.14021153537432351, 0.139532122363363, 0.13882484826701028, 0.13829718294484275, 0.13764841210955664, 0.1371054504212879, 0.13653402777285803, 0.13604278022902352, 0.13561149304367248, 0.13512049031030565, 0.13467753178846267, 0.13428786183311825, 0.13386702964646477, 0.13351642541431247, 0.13313199023973374, 0.13280165745190212, 0.1324541105236326, 0.1320773253009433, 0.13184552876268113, 0.1314723581961223, 0.13118329141253518, 0.13088410946528117, 0.13065526800042107, 0.13036363118376051];
loss(4, :) = [0.23889724267721177, 0.20877331976890565, 0.19658005444606144, 0.18976597041289012, 0.18531927086114883, 0.18239966626962026, 0.1800189746260643, 0.17760722883145014, 0.17571395383675892, 0.17411846752961477, 0.17285464284817378, 0.17181340150435764, 0.17080238041877746, 0.17005873963832854, 0.16924471834500632, 0.16844895120064418, 0.16777130628029505, 0.16714315638542177, 0.16645128346681595, 0.16579737254778545, 0.16510833816130957, 0.16458219501574833, 0.16397970622380575, 0.16343334392706554, 0.16291655864715576, 0.16248258643547694, 0.1619871325890223, 0.16148929001887638, 0.16107567734320957, 0.16069199203252793, 0.16028237388928732, 0.15990190333127977, 0.15956807385285696, 0.15919600078662235, 0.15880092107057572, 0.1584332579255104, 0.15824185002644856, 0.1579098454038302, 0.1576491945664088, 0.15725203816890718];
loss(5, :) = [0.24638170011838278, 0.22053744817972185, 0.21396643701791762, 0.2087904941479365, 0.20502170034646988, 0.20232561569213867, 0.20022150921026866, 0.19866895178953806, 0.19739258242050806, 0.19632373164892197, 0.19533518299659094, 0.19439438966910044, 0.1933724791566531, 0.19256208546161652, 0.19180849371353786, 0.19116609596014023, 0.19067779725790024, 0.19015551414489745, 0.1896956755042076, 0.1893543222506841, 0.18888286412159602, 0.1886109324336052, 0.18825611166159312, 0.18791308792034786, 0.18755063457886378, 0.18725201752583187, 0.18697435064315795, 0.18669430214564006, 0.1863514504154523, 0.18605072348912557, 0.18581593073209127, 0.185450832593441, 0.1851552099307378, 0.18487656616369882, 0.18460370386044184, 0.18433255481322605, 0.18401499707301458, 0.1838493396202723, 0.183520700999101, 0.18325556213061014];








val_lost(1,:) = [0.12598090418676536, 0.1256287168065707, 0.11016786222159862, 0.09354773058493931, 0.0902646634404858, 0.08685967759291331, 0.08032015028595925, 0.08450620961934328, 0.08366523606826862, 0.07993031193315983, 0.07592930718511343, 0.08057555122176806, 0.08152315337210894, 0.07387962449590366, 0.0755090057651202, 0.0731505712022384, 0.0665117098564903, 0.07304557579010724, 0.06521140954146783, 0.06689010488986968, 0.07608922380705674, 0.06711332640051841, 0.06957449675599733, 0.06786161730438471, 0.07018186881393194, 0.06421899056682984, 0.06575387086222569, 0.06553255305935939, 0.07332394805550575, 0.0642288728033503, 0.06335351091374954, 0.06958961219340563, 0.06708995365103086, 0.061740563576420146, 0.06167365163813035, 0.05991661816835404, 0.05954123201221228, 0.0601393416399757, 0.06122053165237109, 0.06143884856502215];
val_lost(2,:) = [0.17148941822052002, 0.15316043315728506, 0.1432251616080602, 0.14203970566590626, 0.12971011520226797, 0.127405424284935, 0.12599133655230205, 0.12012808817227681, 0.12208513791163762, 0.11489318896929424, 0.11257161951462427, 0.1169859027703603, 0.11329868632157644, 0.11039381082057953, 0.10739976547161738, 0.10705758362611136, 0.1089965986331304, 0.1115350493947665, 0.11307317452828089, 0.11329978526830674, 0.11228980070352554, 0.10522172153393428, 0.10498455913464229, 0.1020896887143453, 0.1073804470618566, 0.10305965943733851, 0.10243346723715464, 0.10368256841897965, 0.1033185606678327, 0.10209848102728526, 0.10243638725280761, 0.1018171763976415, 0.10213015331427257, 0.10036350762844086, 0.09801770521799723, 0.10212694668769837, 0.10753585063616435, 0.09955813217957815, 0.09624635602633158, 0.09854659387270609];
val_lost(3,:) = [0.18723516092300416, 0.1796586712439855, 0.15919430799484252, 0.15645080188910165, 0.15456234142780303, 0.14703033929665885, 0.1440857818365097, 0.13954355906645458, 0.13999342908064524, 0.13503959543704985, 0.1322941511074702, 0.13099623448053996, 0.1334685054620107, 0.13279561659495034, 0.13233314122358958, 0.13008998934427898, 0.12718715903759004, 0.12816682583093644, 0.1303101567029953, 0.12678768185774486, 0.12592016541957854, 0.129045689813296, 0.12435629683335622, 0.1264129345536232, 0.12825670335292816, 0.123113991355896, 0.12230006280342738, 0.1240440494219462, 0.12387540321350098, 0.1259557019074758, 0.12539126069148382, 0.12095569180250168, 0.12302548898061116, 0.12449175111452739, 0.1243597288052241, 0.12533584338823955, 0.11873556129535039, 0.1224897422115008, 0.12040665776729584, 0.11973101702928543];
val_lost(4,:) = [0.22109912548065186, 0.20314967203140258, 0.19297968684037525, 0.1877348781267802, 0.18571410868962607, 0.18363214441140494, 0.17782071950435638, 0.17840314993063608, 0.17551592203776042, 0.17382751574516297, 0.17071620790163675, 0.17156510332425434, 0.16495499114195505, 0.16863751296202342, 0.16460202105045318, 0.16735789198875428, 0.16648158167203267, 0.16365830500125886, 0.16105073006947834, 0.16129257018566132, 0.15941689813931784, 0.16266558566888173, 0.1563265276193619, 0.1591515388647715, 0.15842586030960082, 0.16043498734633127, 0.1550199997742971, 0.15847909479141237, 0.15718410081068676, 0.1569009710232417, 0.15550568070411683, 0.15426190898418427, 0.15510865452289582, 0.15503633207480114, 0.15122501514752706, 0.15377741734981537, 0.1504656809091568, 0.15275701944828032, 0.15032602953910829, 0.15171411587397257];
val_lost(5,:) = [0.23001228213310243, 0.22071895221074422, 0.21320784695943196, 0.2065045249859492, 0.20752520790894827, 0.20501928894519805, 0.200215243212382, 0.19834393701553343, 0.2018178415854772, 0.19777225233713785, 0.19630413693586984, 0.19339856181144716, 0.1931287004550298, 0.19278885895411174, 0.19657925546964009, 0.18879714336395265, 0.18919882816473643, 0.19292804333368938, 0.190015189743042, 0.18823314673105876, 0.1847286056836446, 0.1874588452974955, 0.18676933155059813, 0.188026287261645, 0.18391194232304892, 0.1866545171578725, 0.18275740514596303, 0.18208019999663036, 0.1822321426788966, 0.18462763907114665, 0.1825910884062449, 0.181021288951238, 0.1807718030691147, 0.18132456168333688, 0.18277459886868794, 0.17862919745445252, 0.180620472184817, 0.1793598038037618, 0.17736973434289297, 0.1794336346467336];

colors = linspecer(5)

figure
hold on
for k = 1:5
%     plot(lost(k, 1:10), '-', 'Color', colors(k, :), 'LineWidth', 2);
    h(k) = plot(val_lost(k, :), '-', 'Color', colors(k, :), 'LineWidth', 2);
    plot(loss(k, :), '--', 'Color', colors(k, :), 'LineWidth', 2);

end
leg = legend(h, '0', '0.5', '1', '2', '4');
title(leg, 'noise %')
xlabel('epochs')
ylabel('val loss')
ylim([0 0.25])
% set(gca, 'FontSize', 24, 'FontWeight','bold')
 
    
    