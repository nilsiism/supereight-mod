#ifdef __CUDACC__
#define CONST __constant__
#else
#define CONST static const 
#endif
CONST float bspline_num_samples = 1000;
CONST float bspline_lookup[1000] = {0.0,4.51352704507e-09,3.61082163605e-08,1.21865230217e-07,2.88865730884e-07,5.64190880633e-07,9.74921841735e-07,1.54813977646e-06,2.31092584707e-06,3.29036121585e-06,4.51352704507e-06,6.00750449698e-06,7.79937473388e-06,9.91621891801e-06,1.23851182117e-05,1.52331537771e-05,1.84874067766e-05,2.21749583724e-05,2.63228897268e-05,3.09582820021e-05,3.61082163605e-05,4.17997739644e-05,4.80600359759e-05,5.49160835573e-05,6.2394997871e-05,7.05238600792e-05,7.93297513441e-05,8.88397528281e-05,9.90809456933e-05,0.000110080411102,0.000121865230217,0.0001344624842,0.000147899254213,0.000162202621419,0.000177399666979,0.000193517472057,0.000210583117815,0.000228623685414,0.000247666256017,0.000267737910786,0.000288865730884,0.000311076797473,0.000334398191715,0.000358856994772,0.000384480287807,0.000411295151982,0.000439328668459,0.0004686079184,0.000499159982968,0.000531011943325,0.000564190880633,0.000598723876055,0.000634638010753,0.000671960365889,0.000710718022625,0.000750938062123,0.000792647565547,0.000835873614057,0.000880643288817,0.000926983670989,0.000974921841735,0.00102448488222,0.0010756998736,0.00112859389704,0.0011831940337,0.00123952736475,0.00129762097135,0.00135750193466,0.00141919733583,0.00148273425605,0.00154813977646,0.00161544097823,0.00168466494252,0.00175583875049,0.00182898948331,0.00190414422214,0.00198133004814,0.00206057404247,0.00214190328629,0.00222534486077,0.00231092584707,0.00239867332636,0.00248861437978,0.00258077608852,0.00267518553372,0.00277186979655,0.00287085595818,0.00297217109976,0.00307584230246,0.00318189664743,0.00329036121585,0.00340126308888,0.00351462934767,0.00363048707339,0.0037488633472,0.00386978525026,0.00399327986374,0.0041193742688,0.0042480955466,0.0043794707783,0.00451352704507,0.00465029142806,0.00478979100844,0.00493205286738,0.00507710408602,0.00522497174555,0.00537568292711,0.00552926471187,0.005685744181,0.00584514841565,0.00600750449698,0.00617283950617,0.00634118052437,0.00651255463275,0.00668698891246,0.00686451044467,0.00704514631054,0.00722892359123,0.00741586936791,0.00760601072174,0.00779937473388,0.00799598848549,0.00819587905773,0.00839907353177,0.00860559898877,0.0088154825099,0.0090287511763,0.00924543206916,0.00946555226962,0.00968913885885,0.00991621891801,0.0101468195283,0.0103809677708,0.0106186907267,0.0108600154772,0.0111049691035,0.0113535786867,0.0116058713079,0.0118618740484,0.0121216139892,0.0123851182117,0.0126524137968,0.0129235278258,0.0131984873799,0.0134773195401,0.0137600513878,0.0140467100039,0.0143373224698,0.0146319158665,0.0149305172752,0.0152331537771,0.0155398524533,0.0158506403851,0.0161655446535,0.0164845923397,0.016807810525,0.0171352262903,0.017466866717,0.0178027588862,0.018142929879,0.0184874067766,0.0188362166602,0.0191893866109,0.0195469437098,0.0199089150383,0.0202753276773,0.0206462087081,0.0210215852119,0.0214014842698,0.0217859329629,0.0221749583724,0.0225685875795,0.0229668476654,0.0233697657112,0.0237773687981,0.0241896840072,0.0246067384197,0.0250285591167,0.0254551731795,0.0258866076891,0.0263228897268,0.0267640463737,0.027210104711,0.0276610918198,0.0281170347814,0.0285779606767,0.0290438965871,0.0295148695937,0.0299909067776,0.03047203522,0.0309582820021,0.031449674205,0.03194623891,0.032448003198,0.0329549941504,0.0334672388483,0.0339847643728,0.0345075978051,0.0350357662264,0.0355692967178,0.0361082163605,0.0366525522357,0.0372023314245,0.037757581008,0.0383183280675,0.0388845996841,0.039456422939,0.0400338249133,0.0406168326882,0.0412054733448,0.0417997739644,0.042399761628,0.0430054634169,0.0436169064121,0.044234117695,0.0448571243465,0.045485953448,0.0461206320805,0.0467611873252,0.0474076462633,0.0480600359759,0.0487183835442,0.0493827160494,0.0500530605726,0.050729444195,0.0514118939977,0.052100437062,0.0527951004689,0.0534959112997,0.0542028966354,0.0549160835573,0.0556354991466,0.0563611704843,0.0570931246517,0.0578313887299,0.0585759898,0.0593269549433,0.0600843112409,0.0608480857739,0.0616183056236,0.062394997871,0.0631781895974,0.0639679078839,0.0647641798117,0.0655670324619,0.0663764929156,0.0671925882542,0.0680153455587,0.0688447919102,0.06968095439,0.0705238600792,0.0713735360589,0.0722300094104,0.0730933072148,0.0739634565533,0.0748404845069,0.0757244181569,0.0766152845845,0.0775131108708,0.0784179240969,0.0793297513441,0.0802486196935,0.0811745562262,0.0821075880234,0.0830477421663,0.0839950457361,0.0849495258138,0.0859112094807,0.086880123818,0.0878562959067,0.0888397528281,0.0898305216633,0.0908286294934,0.0918341033997,0.0928469704633,0.0938672577654,0.0948949923871,0.0959302014096,0.096972911914,0.0980231509815,0.0990809456933,0.100146323131,0.101219310374,0.102299934506,0.103388222607,0.104484201757,0.105587899039,0.106699341533,0.107818556321,0.108945570484,0.110080411102,0.111223105258,0.112373680031,0.113532162505,0.114698579758,0.115872958874,0.117055326932,0.118245711014,0.119444138202,0.120650635575,0.121865230217,0.123087949207,0.124318819627,0.125557868558,0.126805123081,0.128060610277,0.129324357228,0.130596391014,0.131876738718,0.133165427419,0.1344624842,0.13576793614,0.137081810323,0.138404133827,0.139734933736,0.14107423713,0.142422071089,0.143778462696,0.145143439032,0.146517027177,0.147899254213,0.149290147221,0.150689733281,0.152098039476,0.153515092887,0.154940920594,0.156375549679,0.157819007222,0.159271320306,0.160732516011,0.162202621419,0.16368166361,0.165169669665,0.166666666667,0.168172668155,0.169687633507,0.171211508562,0.172744239157,0.17428577113,0.175836050319,0.17739502256,0.178962633692,0.180538829553,0.182123555979,0.18371675881,0.185318383882,0.186928377033,0.188546684101,0.190173250923,0.191808023338,0.193450947182,0.195101968294,0.196761032511,0.198428085671,0.200103073612,0.20178594217,0.203476637185,0.205175104493,0.206881289933,0.208595139341,0.210316598556,0.212045613415,0.213782129757,0.215526093417,0.217277450236,0.219036146049,0.220802126695,0.222575338011,0.224355725835,0.226143236004,0.227937814358,0.229739406732,0.231547958965,0.233363416894,0.235185726357,0.237014833192,0.238850683237,0.240693222328,0.242542396304,0.244398151003,0.246260432262,0.248129185918,0.25000435781,0.251885893776,0.253773739652,0.255667841276,0.257568144487,0.259474595121,0.261387139017,0.263305722012,0.265230289944,0.267160788651,0.26909716397,0.271039361738,0.272987327794,0.274941007976,0.27690034812,0.278865294065,0.280835791648,0.282811786708,0.28479322508,0.286780052604,0.288772215117,0.290769658457,0.292772328461,0.294780170967,0.296793131813,0.298811156836,0.300834191874,0.302862182765,0.304895075346,0.306932815455,0.30897534893,0.311022621608,0.313074579327,0.315131167925,0.31719233324,0.319258021109,0.321328177369,0.323402747859,0.325481678416,0.327564914877,0.329652403082,0.331744088866,0.333839918068,0.335939836526,0.338043790077,0.340151724559,0.342263585809,0.344379319665,0.346498871966,0.348622188548,0.350749215249,0.352879897907,0.35501418236,0.357152014444,0.359293339999,0.361438104862,0.363586254869,0.36573773586,0.367892493671,0.370050474141,0.372211623107,0.374375886406,0.376543209877,0.378713539356,0.380886820682,0.383062999693,0.385242022226,0.387423834118,0.389608381208,0.391795609333,0.393985464331,0.396177892039,0.398372838295,0.400570248937,0.402770069802,0.404972246728,0.407176725554,0.409383452115,0.411592372251,0.413803431799,0.416016576596,0.41823175248,0.42044890529,0.422667980861,0.424888925033,0.427111683643,0.429336202528,0.431562427527,0.433790304476,0.436019779214,0.438250797579,0.440483305407,0.442717248537,0.444952572806,0.447189224052,0.449427148112,0.451666290825,0.453906598028,0.456148015559,0.458390489255,0.460633964954,0.462878388493,0.465123705711,0.467369862445,0.469616804533,0.471864477812,0.47411282812,0.476361801295,0.478611343174,0.480861399595,0.483111916397,0.485362839415,0.487614114489,0.489865687455,0.492117504152,0.494369510417,0.496621652088,0.498873875002,0.501126124998,0.503378347912,0.505630489583,0.507882495848,0.510134312545,0.512385885511,0.514637160585,0.516888083603,0.519138600405,0.521388656826,0.523638198705,0.52588717188,0.528135522188,0.530383195467,0.532630137555,0.534876294289,0.537121611507,0.539366035046,0.541609510745,0.543851984441,0.546093401972,0.548333709175,0.550572851888,0.552810775948,0.555047427194,0.557282751463,0.559516694593,0.561749202421,0.563980220786,0.566209695524,0.568437572473,0.570663797472,0.572888316357,0.575111074967,0.577332019139,0.57955109471,0.58176824752,0.583983423404,0.586196568201,0.588407627749,0.590616547885,0.592823274446,0.595027753272,0.597229930198,0.599429751063,0.601627161705,0.603822107961,0.606014535669,0.608204390667,0.610391618792,0.612576165882,0.614757977774,0.616937000307,0.619113179318,0.621286460644,0.623456790123,0.625624113594,0.627788376893,0.629949525859,0.632107506329,0.63426226414,0.636413745131,0.638561895138,0.640706660001,0.642847985556,0.64498581764,0.647120102093,0.649250784751,0.651377811452,0.653501128034,0.655620680335,0.657736414191,0.659848275441,0.661956209923,0.664060163474,0.666160081932,0.668255911134,0.670347596918,0.672435085123,0.674518321584,0.676597252141,0.678671822631,0.680741978891,0.68280766676,0.684868832075,0.686925420673,0.688977378392,0.69102465107,0.693067184545,0.695104924654,0.697137817235,0.699165808126,0.701188843164,0.703206868187,0.705219829033,0.707227671539,0.709230341543,0.711227784883,0.713219947396,0.71520677492,0.717188213292,0.719164208352,0.721134705935,0.72309965188,0.725058992024,0.727012672206,0.728960638262,0.73090283603,0.732839211349,0.734769710056,0.736694277988,0.738612860983,0.740525404879,0.742431855513,0.744332158724,0.746226260348,0.748114106224,0.74999564219,0.751870814082,0.753739567738,0.755601848997,0.757457603696,0.759306777672,0.761149316763,0.762985166808,0.764814273643,0.766636583106,0.768452041035,0.770260593268,0.772062185642,0.773856763996,0.775644274165,0.777424661989,0.779197873305,0.780963853951,0.782722549764,0.784473906583,0.786217870243,0.787954386585,0.789683401444,0.791404860659,0.793118710067,0.794824895507,0.796523362815,0.79821405783,0.799896926388,0.801571914329,0.803238967489,0.804898031706,0.806549052818,0.808191976662,0.809826749077,0.811453315899,0.813071622967,0.814681616118,0.81628324119,0.817876444021,0.819461170447,0.821037366308,0.82260497744,0.824163949681,0.82571422887,0.827255760843,0.828788491438,0.830312366493,0.831827331845,0.833333333333,0.834830330335,0.83631833639,0.837797378581,0.839267483989,0.840728679694,0.842180992778,0.843624450321,0.845059079406,0.846484907113,0.847901960524,0.849310266719,0.850709852779,0.852100745787,0.853482972823,0.854856560968,0.856221537304,0.857577928911,0.85892576287,0.860265066264,0.861595866173,0.862918189677,0.86423206386,0.8655375158,0.866834572581,0.868123261282,0.869403608986,0.870675642772,0.871939389723,0.873194876919,0.874442131442,0.875681180373,0.876912050793,0.878134769783,0.879349364425,0.880555861798,0.881754288986,0.882944673068,0.884127041126,0.885301420242,0.886467837495,0.887626319969,0.888776894742,0.889919588898,0.891054429516,0.892181443679,0.893300658467,0.894412100961,0.895515798243,0.896611777393,0.897700065494,0.898780689626,0.899853676869,0.900919054307,0.901976849018,0.903027088086,0.90406979859,0.905105007613,0.906132742235,0.907153029537,0.9081658966,0.909171370507,0.910169478337,0.911160247172,0.912143704093,0.913119876182,0.914088790519,0.915050474186,0.916004954264,0.916952257834,0.917892411977,0.918825443774,0.919751380307,0.920670248656,0.921582075903,0.922486889129,0.923384715415,0.924275581843,0.925159515493,0.926036543447,0.926906692785,0.92776999059,0.928626463941,0.929476139921,0.93031904561,0.93115520809,0.931984654441,0.932807411746,0.933623507084,0.934432967538,0.935235820188,0.936032092116,0.936821810403,0.937605002129,0.938381694376,0.939151914226,0.939915688759,0.940673045057,0.9414240102,0.94216861127,0.942906875348,0.943638829516,0.944364500853,0.945083916443,0.945797103365,0.9465040887,0.947204899531,0.947899562938,0.948588106002,0.949270555805,0.949946939427,0.950617283951,0.951281616456,0.951939964024,0.952592353737,0.953238812675,0.95387936792,0.954514046552,0.955142875653,0.955765882305,0.956383093588,0.956994536583,0.957600238372,0.958200226036,0.958794526655,0.959383167312,0.959966175087,0.960543577061,0.961115400316,0.961681671932,0.962242418992,0.962797668576,0.963347447764,0.963891783639,0.964430703282,0.964964233774,0.965492402195,0.966015235627,0.966532761152,0.96704500585,0.967551996802,0.96805376109,0.968550325795,0.969041717998,0.96952796478,0.970009093222,0.970485130406,0.970956103413,0.971422039323,0.971882965219,0.97233890818,0.972789895289,0.973235953626,0.973677110273,0.974113392311,0.974544826821,0.974971440883,0.97539326158,0.975810315993,0.976222631202,0.976630234289,0.977033152335,0.97743141242,0.977825041628,0.978214067037,0.97859851573,0.978978414788,0.979353791292,0.979724672323,0.980091084962,0.98045305629,0.980810613389,0.98116378334,0.981512593223,0.981857070121,0.982197241114,0.982533133283,0.98286477371,0.983192189475,0.98351540766,0.983834455347,0.984149359615,0.984460147547,0.984766846223,0.985069482725,0.985368084134,0.98566267753,0.985953289996,0.986239948612,0.98652268046,0.98680151262,0.987076472174,0.987347586203,0.987614881788,0.987878386011,0.988138125952,0.988394128692,0.988646421313,0.988895030896,0.989139984523,0.989381309273,0.989619032229,0.989853180472,0.990083781082,0.990310861141,0.99053444773,0.990754567931,0.990971248824,0.99118451749,0.991394401011,0.991600926468,0.991804120942,0.992004011515,0.992200625266,0.992393989278,0.992584130632,0.992771076409,0.992954853689,0.993135489555,0.993313011088,0.993487445367,0.993658819476,0.993827160494,0.993992495503,0.994154851584,0.994314255819,0.994470735288,0.994624317073,0.994775028254,0.994922895914,0.995067947133,0.995210208992,0.995349708572,0.995486472955,0.995620529222,0.995751904453,0.995880625731,0.996006720136,0.99613021475,0.996251136653,0.996369512927,0.996485370652,0.996598736911,0.996709638784,0.996818103353,0.996924157698,0.9970278289,0.997129144042,0.997228130203,0.997324814466,0.997419223911,0.99751138562,0.997601326674,0.997689074153,0.997774655139,0.997858096714,0.997939425958,0.998018669952,0.998095855778,0.998171010517,0.99824416125,0.998315335057,0.998384559022,0.998451860224,0.998517265744,0.998580802664,0.998642498065,0.998702379029,0.998760472635,0.998816805966,0.998871406103,0.998924300126,0.998975515118,0.999025078158,0.999073016329,0.999119356711,0.999164126386,0.999207352434,0.999249061938,0.999289281977,0.999328039634,0.999365361989,0.999401276124,0.999435809119,0.999468988057,0.999500840017,0.999531392082,0.999560671332,0.999588704848,0.999615519712,0.999641143005,0.999665601808,0.999688923203,0.999711134269,0.999732262089,0.999752333744,0.999771376315,0.999789416882,0.999806482528,0.999822600333,0.999837797379,0.999852100746,0.999865537516,0.99987813477,0.999889919589,0.999900919054,0.999911160247,0.999920670249,0.99992947614,0.999937605002,0.999945083916,0.999951939964,0.999958200226,0.999963891784,0.999969041718,0.99997367711,0.999977825042,0.999981512593,0.999984766846,0.999987614882,0.999990083781,0.999992200625,0.999993992496,0.999995486473,0.999996709639,0.999997689074,0.99999845186,0.999999025078,0.999999435809,0.999999711134,0.999999878135,0.999999963892,0.999999995486,1.0};
