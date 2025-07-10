use pyo3::prelude::*;
use pyo3::types::PyInt;
use pyo3::prelude::PyAnyMethods;
use num_bigint::BigInt;
use num_bigint::RandBigInt;
use std::str::FromStr;
use num_traits::ToPrimitive;
use crossbeam::channel;
use std::collections::HashMap;
use hex;
use num_integer::Integer;

type Vector = Vec<BigInt>;

// DDH Struct
#[pyclass]
#[derive(Clone)]
struct DdhParams {
    l:      BigInt, // Length of input vectors
    bound:  BigInt, // Bound for the values
    g:      BigInt, // Generator of the cyclic group
    p:      BigInt, // Modulus
    q:      BigInt, //Order of Generator (always P-1)
}

#[pyclass]
pub struct Ddh {
    params: DdhParams,
}

#[pyclass]
#[derive(Clone)]
pub struct UniformRange {
	min: BigInt,
	max: BigInt,
}

#[pyclass]
#[derive(Clone)]
pub struct CalcZp {
    p: BigInt,
    bound: BigInt,
    m: BigInt,
    neg: bool,
}

#[pyfunction]
fn check_bound(v: Vector, bound: Bound<'_, PyInt>) -> PyResult<bool> {
    let bound: BigInt = bound.extract()?;
    for value in v.iter() {
        if value > &bound {
            return Ok(false);
        }
    }
    Ok(true)
}

#[pyfunction]
fn dot(v1: Vector, v2: Vector, py: Python<'_>) -> Bound<'_, PyInt> {
    if v1.len() != v2.len() {
        panic!("Vectors must be of the same length");
    }

    let mut result = BigInt::from(0);
    for (a, b) in v1.iter().zip(v2.iter()) {
        result += a * b;
    }

    result.into_pyobject(py).unwrap()
}

#[pymethods]
impl UniformRange {
    #[new]
    fn new(min: Bound<'_, PyInt>, max: Bound<'_, PyInt>) -> PyResult<UniformRange> {
            let min:BigInt = min.extract().unwrap();
            let max:BigInt = max.extract().unwrap();
            Ok(UniformRange {max: max, min:min})
    }

    fn sample<'py>(&self, py: Python<'py>) -> Bound<'py, PyInt> {
        let max:BigInt = self.max.clone();
        let min:BigInt = self.min.clone();
        let range:BigInt = &max - &min;
        let mut rng = rand::thread_rng();
        let random_value = rng.gen_bigint_range(&BigInt::from(0u32), &range) + &min;
        random_value.into_pyobject(py).unwrap()
    }
}

#[pyfunction]
fn in_zp(p: BigInt, order: BigInt) -> CalcZp {
    let one = BigInt::from(1);
    let bound = order.clone();
    let m = bound.sqrt() + one;
    CalcZp {
        p: p,
        bound: bound,
        m: m,
        neg: true,
    }
}

pub fn baby_step_giant_step(calc: CalcZp, h: BigInt, g: BigInt, ) -> Result<BigInt, String> {    
    let zero = BigInt::from(0);
    let one = BigInt::from(1);
    let two = BigInt::from(2);
    let (tx, rx) = channel::unbounded::<Result<BigInt, String>>();

    // First attempt with g
    let mut found = false;  // Add this to track if we found a solution
    let mut t: HashMap<String, BigInt> = HashMap::new();
    let mut x = one.clone();
    let mut y = h.clone();
    let mut z = g.modinv(&calc.p).expect("Inverse does not exist");
    z = z.modpow(&two, &calc.p);
    let bitlen: u32 = calc.bound.bits().to_u32().expect("bitlen too large for u32");  // Use calc.bound instead of x

    let key = hex::encode(x.clone().to_bytes_be().1); 
    t.insert(key.clone(), zero.clone());
    x = (x * g.clone()).mod_floor(&calc.p);
    let mut j = zero.clone();
    let mut giant_step = BigInt::from(0);
    let mut bound = BigInt::from(0);

    // First attempt loop
    for i in 0..bitlen {
        if found { break; }  // Exit loop if we found a solution

        giant_step = two.pow(i+1);
        if giant_step > calc.m {
            giant_step = calc.m.clone();
            z = g.modinv(&calc.p).expect("Inverse does not exist");
            z = z.clone().modpow(&calc.m, &calc.p);
        }

        let val = two.pow(i);
        let mut k = val.clone();
        while k < giant_step {
            let key = hex::encode(x.clone().to_bytes_be().1); 
            t.insert(key, k.clone());
            x = (x*g.clone()).mod_floor(&calc.p);
            k = k + one.clone();
        }

        bound = two.pow(2*(i+1));
        while j < bound {
            let current_key = hex::encode(y.clone().to_bytes_be().1);
            if let Some(e) = t.get(&current_key) {
                let result = j.clone() + e.clone();
                tx.send(Ok(result)).unwrap();
                found = true;  // Mark as found
                break;
            }
            y = (y.clone() * z.clone()).mod_floor(&calc.p);
            j = j + giant_step.clone();
        }
        
        if !found {
            z = z.modpow(&two, &calc.p);
        }
    }
    
    // Only send error if we didn't find anything
    if !found {
        tx.send(Err("failed to find the discrete logarithm within bound".to_string())).unwrap();
    }

    // Second attempt with g^-1
    found = false;  // Reset flag
    let mut ginv = g.clone().modinv(&calc.p).expect("Inverse does not exist");
    
    let mut t: HashMap<String, BigInt> = HashMap::new();
    x = one.clone();
    y = h.clone();
    z = ginv.modinv(&calc.p).expect("Inverse does not exist");
    z = z.clone().modpow(&two, &calc.p);
    let bitlen:u32 = x.bits().to_u32().expect("bitlen too large for u32");

    let key = hex::encode(x.to_bytes_be().1); 
    t.insert(key.clone(), zero.clone());
    x = (x * ginv.clone()).mod_floor(&calc.p);
    let mut j = zero.clone();
    let mut giant_step = BigInt::from(0);
    let mut bound = BigInt::from(0);

    for i in 0..bitlen {
        if found { break; }  // Exit loop if we found a solution

        giant_step = two.pow(i+1);
        if giant_step>calc.m{
            giant_step = calc.m.clone();
            z = ginv.modinv(&calc.p).expect("Inverse does not exist");
            z = z.modpow(&calc.m, &calc.p);
        }

        let i_big = BigInt::from(i);
        let val = two.pow(i);
        let mut k = val.clone();
        while k < giant_step {
            let key = hex::encode(x.to_bytes_be().1); 
            t.insert(key, k.clone());
            x = (x.clone()*ginv.clone()).mod_floor(&calc.p.clone());
            k = k.clone() + one.clone();
        }

        bound = two.pow(2*(i+1));
        while j < bound {
            let current_key = hex::encode(y.clone().to_bytes_be().1);
            if let Some(e) = t.get(&current_key) {
                let result = j.clone() + e;
                tx.send(Ok(result)).unwrap();
                found = true;  // Mark as found
                break;
            }
            y = (y.clone() * z.clone()).mod_floor(&calc.p);
            j = j + giant_step.clone();
        }
        
        if !found {
            z = z.modpow(&two, &calc.p);
        }
    }
    
    // Only send error if we didn't find anything
    if !found {
        tx.send(Err("failed to find the discrete logarithm within bound".to_string())).unwrap();
    }

    // Process results
    let res1 = rx.recv().unwrap();
    match res1 {
        Ok(val) => {
            // Test if g^val = h mod p
            let g_pow = g.modpow(&val, &calc.p);
            let diff = if g_pow == h { 
                val 
            } else { 
                -val 
            };
            return Ok(diff);
        },
        Err(_) => {
            // If first failed, try second
            let res2 = rx.recv().unwrap();
            match res2 {
                Ok(val) => {
                    let g_pow = g.modpow(&val, &calc.p);
                    return Ok(if g_pow == h { val } else { -val });
                },
                Err(e2) => {
                    return Err(format!("Both attempts failed: {}", e2));
                }
            }
        }
    }
}

// Create a new ddh struct.
#[pymethods]
impl Ddh {
    #[new]
    fn new(l: Bound<'_, PyInt>, bound: Bound<'_, PyInt>, g: Bound<'_, PyInt>, p: Bound<'_, PyInt>, q: Bound<'_, PyInt>) -> PyResult<Self> {
        let l: BigInt = l.extract()?;
        let bound: BigInt = bound.extract()?;
        let g: BigInt = g.extract()?;
        let p: BigInt = p.extract()?;
        let q: BigInt = q.extract()?;

        Ok(Ddh {
            params: DdhParams {l, bound, g, p, q,},
        })
    }

    #[staticmethod]
    fn new_ddh_precomp(py: Python<'_>, l: Bound<'_, PyInt>, modulus_length:usize, bound: Bound<'_, PyInt>) -> PyResult<Self> {
        let one = BigInt::from(1);
        let two = BigInt::from(2);
        let (g,p) = match modulus_length {
            1024 => (BigInt::from_str("34902160241479276675539633849372382885917193816560610471607073855548755350834003692485735908635894317735639518678334280193650806072183057417077181724192674928134805218882803812978345229222559213790765817899845072682155064387311523738581388872686127675360979304234957611566801734164757915959042140104663977828").expect("invalid big int literal"), 
            BigInt::from_str("166211269243229118758738154756726384542659478479960313411107431885216572625212662756677338184675400324411541201832214281445670912135683272416408753424543622705770319923251281963485084208425069817917631106045349238686234860629044433560424091289406000897029571960128048529362925472176997104870527051276406995203").expect("invalid big int literal")),
            1536 => (BigInt::from_str("676416913692519694440150163403654362412279108516867264953779609011365998625435399420336578530015558254310139891236630566729665914687641028600402606957815727025192669238117788115237116562468680376464346714542467465836552396661693422160454402926392749202926871877212792118140354124110927269910674002861908621272286950597240072605316784317536178700101838123530590145680002962405974024190384775185108002307650499125333676880320808656556635493186351335151559453463208").expect("invalid big int literal"),
            BigInt::from_str("1851297899986638926486011430658634631676522135433726749065856232802142091866650774719879427474637700607873256035038534449089405369134066444876856913629831069906506096279113968447116822488133963417347136141052507685108634240736100862550194947326287783557220764070479431781692630708747550712729778398000353165406458520850089303530985563143326919073190605085889925484113854496074216626577246143598303709289292397203458923541841135799203967503522114881404128535647507").expect("invalid big int literal")),
            2048 => (BigInt::from_str("4006960929413042209594165215465319088439374252008797022450541422457034721098828647078778605657155669917104962611933792130890703423519992986737966991597160684973795472419962788730248050852176194215699504914899438223683843401963466624139534923052671383315398134823370041633710463630745156269175253639670460050105594663691338308037509280576148624454011047879615100156717631945194107791315234171086603775159708325087759679758438868772220133433497821899045165244202228696902434100209752952701657306825368599999359102329396520012735146260911352901326915877502873633420811221206110021993351144711002138373506576799781061829").expect("invalid big int literal"),
            BigInt::from_str("28884237504713658990682089080899862128005980675308910325841161962760155725037929764087367167449843609136681034352509183117742758446654629096509285354423361556493020266963222508540306384896802796001914743293196010488452478370041404523014215612960481024232879327123268440037633547483165934132901270561772860319969916305482525766132307669097012989986613879246932730824899649301621408341438037745468033187743673001187803377254713546325789438300798311106106322698517805307792059495696632070953526611920926003483451787562399452650878943515646786958216714025307572678422373120397225912926110031401983688860264234966561627699").expect("invalid big int literal")),
            2560 => (BigInt::from_str("283408881721750179985507845260248881237898607313021593637913985490973593382472548378053368228310040484438920416918021737085067966444840449068073874437662089659563355479608930182263107110969154912883370207053052795964658868443319273167870311282045348396320159912742394374241864712808383029954025256232806465551969466207671603658677963161975454703127476120201164519187150268352527923664649275471494757270139533433456630363925187498055211365480086561354743681517539297815712218419607006668655891574362066382949706266666189227897710299445185100212256741698216505337617571970963008519334554537811591236478130526432239803909461119767954934793813410765013072006162612226471775059215628326278458577643374735250370115470812597459244082296191871275203831471332697557979904062571849").expect("invalid big int literal"),
            BigInt::from_str("403126381462544353337185732949672984701995200694926494258456907009253008321275627278199160008680481542251933923210212978304159426174307419863781623411302777276318286800690060166638633211627521530173324015527027650548199185539205697958056639406116068885574865579676651743820636201007864067569576455725489531113260031526605827601510665037511961715114944815619491261828558745083975042855063688267346905844510423020844412350570902289599734320004108557140241966071165594059732527795488131297017205383953304055105007982366596746708951250486384299368612656872813778220074826250625689603663742175288397398948456522281031888042417278385238985218731264092285591879578299600853004336936458454638992426900228708418575870946630137618851131144232868141478901063119847104013555395370887").expect("invalid big int literal")),
            3072 => (BigInt::from_str("3696261717511684685041982088526264061294500298114921057816954458306445697150348237912729036967670872345042594238223317055749478029025374644864924550052402546275985983344583674703146236623453822520422465163020824494790581472736649085281450730460445260696087738043872307629635997875332076478424042345012004769107421873566499123042621978973433575500345010912635742477932006291250637245855027695943163956584173316781442078828050076620331751405548730676363847526959436516279320074682721438642683731766682502490935962962293815202487144775533102010333956118641968798500514719248831145108532912211817219793191951880318961073149276914867129023978524587935704313755469570162971499124682746476415187933097132047611840762510892175328320025164466873845777990557296853549970943298347924080102740724512079409979152285019931666423541870247789529268168448010024121369388707140296446100906359619586133848407970098685310317291828335700424602208").expect("invalid big int literal"),
            BigInt::from_str("4387756306134544957818467663802660683665166110605728231080818705443663402154316615145921798856363268744945754470238000282108344905251127487705736550297997444150840902348669718478564904142834154197029830975532074167513046443903186309497214496864577129616824062991068960005865144004932069025136224356325248036029606434443391988386519658751798077031844645051726026696307027395796695909035405241040411794836124123435225690961994089776517262574417789067836840997650095451062948856617211542724543995145259735683916440579956961657374517806591607068842498749297993409884001044324428640569001916341503645559748760311343179943896427393009949062735145363544745972252566600994034655540841225414736222780096833045470605544717177880459300618917961703559234544541206877026518430276932498602360341258899345739335298856394124351357206871568254540730107127298623178526868418799471896060015463201459762913197633841160710893895836663035998106119").expect("invalid big int literal")),
            4096 => (BigInt::from_str("51665588681810560577916524923861643358980285220048008212528567741884121491554604183472728540139463099618903178110360757930742372390027135064809646425064896539133721148335557788263239281487173350543811713890328584918216783142094297306639941000480756707312457878765754357205186485080839623690156744636468433787780205323460166423447602447200754978133176713947189000663528355089645281397174452923418212485422962705227706103188302892660448134233848971142570881089940852441776074246332915421265800026335300100610273942459340241610730244726628211914068945587128124478812632725838440727321816905181830592204023095726270782834020990986443265625389712733369116937470448592846480352222814297792606318850361699893703272484112273500581408730519942517586496563772194165844831300501908379990979449691597045730512107756238377635183257797115883839801779086058652272455400286891699445584526719648220045380141260347316315487340493029966105973850214850475440630205768783542021741101804842248602349004364816943429122368563644935802417389995380389429997320053299323220481603252879925927515844929958940305561718295197935926645561977544440676439150126025681320050786964708227836328341875446457912905977470123640014345655062829575775837287500880054558386787").expect("invalid big int literal"),
            BigInt::from_str("1022249395832567838406986294560330159176972202126664245047364146720891252715766488477689126342364655087193411078517616569887825896401401223927363505007778278205623713273194552498760148834874746839752870298152746450585455651115247220867383465863156721401567161663838310658875672995951663020449772454232797368263754624173026584111779206080723120076751471597509403139249260220696195263597156452889920392585797464801375940661326779247976331028637271512085826066667631423502199894046717721786935806581428328491087482664043743281068318459302242239861275878019857365021173868449409246193470959347916848019032536247915451026158871684654213802886886213841729258073333569276986893577214659899227179735448593265633219968622571880602115519942763955551007919826002851866939641065270816032435114864853636918330698605282572789904941484540512478406984407320963402583009124880812235841866246441862987563989772424040933513333746472128494254253767426962063553015635240386636751473945937412527996558505231385625318878887383161350102080329744822052478052004574860361461762694379860797225344866320388590336321515376486033237159694567932935601775209663052272120524337888258857351777348841323194553467226791591208931619058871750498804369190487499494069660723").expect("invalid big int literal")),
            _ => panic!("Modulus length should be one of values 1024, 1536, 2048, 2560, 3072, or 4096")
        };

        let q = &p - one;
        let q = q / two;
        let bigl:BigInt = l.extract()?;
        let bigbound:BigInt = bound.extract()?;
        let calc:BigInt = bigbound.pow(2u32);

        if (2 * bigl) * calc > q {
            panic!("2 * l * bound^2 should be smaller than group order")
        }

        let g = g.into_pyobject(py).unwrap();
        let p = p.into_pyobject(py).unwrap();
        let q = q.into_pyobject(py).unwrap();

        Self::new(l, bound, g, p, q)
    }

    fn generate_master_keys(&self, py: Python<'_>) -> (Vector, Vector) {
        let mut master_sec_key: Vector = Vec::new();
        let mut master_pub_key: Vector = Vec::new();
        let q = (self.params.q.clone()).into_pyobject(py).unwrap();
        let sampler = UniformRange::new(PyInt::new(py, 2), q).unwrap();

        // Convert BigInt to usize (if possible)
        let l_usize = self.params.l.to_usize().expect("l does not fit in usize");
        for i in 0..l_usize {
            let x: BigInt = UniformRange::sample(&sampler.clone(), py).extract::<BigInt>().unwrap();
            master_sec_key.push(x.clone());
            let y: BigInt = (&self.params.g).modpow(&x, &self.params.p);
            master_pub_key.push(y);
        }
        (master_sec_key, master_pub_key)
    }

    fn derive_key<'py>(&self, master_sec_key: Vector, y: Vector, py: Python<'py>) -> Bound<'py, PyInt> {
        let bound = (self.params.bound.clone()).into_pyobject(py).unwrap();
        if check_bound(y.clone(), bound).unwrap() == false {
            panic!("y values are not within the specified bound");
        }

        let derived_key = dot(master_sec_key.clone(), y.clone(), py).into_pyobject(py).unwrap();
        derived_key
    }

    fn encrypt<'py>(&self, x: Vector, master_pub_key: Vector, py: Python<'py>) -> Vector {
        let bound = (self.params.bound.clone()).into_pyobject(py).unwrap();
        if check_bound(x.clone(), bound).unwrap() == false {
            panic!("y values are not within the specified bound");
        }

        let zero = BigInt::from(0);
        let q = (self.params.q.clone()).into_pyobject(py).unwrap();
        let sampler = UniformRange::new(PyInt::new(py, 2), q).unwrap();

        let r: BigInt = UniformRange::sample(&sampler, py).extract::<BigInt>().unwrap();

        let mut ciphertext: Vector = Vec::with_capacity(x.len() + 1);

        ciphertext.push((&self.params.g).modpow(&r, &self.params.p));
        for i in 0..x.len() {
            let t1 = (master_pub_key[i].clone()).modpow(&r, &self.params.p);
            let t2:BigInt = if x[i] < zero {
                let x  = -(x[i].clone());
                (&self.params.g).modpow(&x, &self.params.p).modinv(&self.params.p).expect("Inverse does not exist")
            } else {
                (&self.params.g).modpow(&x[i], &self.params.p)
            };
            let t = t1 * t2;
            ciphertext.push(t.mod_floor(&self.params.p));
        }

        ciphertext
    }

    fn decrypt<'py>(&self, ciphertext: Vector, key: Bound<'py, PyInt>, y: Vector, py: Python<'py>) -> Bound<'py, PyInt>{
        let bound = (self.params.bound.clone()).into_pyobject(py).unwrap();
        let key: BigInt = key.extract().unwrap();
        let zero = BigInt::from(0);

        // Fix for vector length mismatch - skip the first element in ciphertext
        let ciphertext_without_first = ciphertext.iter().skip(1).cloned().collect::<Vec<_>>();
        println!("Ciphertext length: {}, y length: {}", ciphertext.len(), y.len());
        println!("Expected inner product: {}", dot(y.clone(), ciphertext_without_first, py));
        
        let mut num = BigInt::from(1);
        for (i, ct) in ciphertext.iter().skip(1).enumerate() {
            let t1 = if y[i] < zero {
                let y = -(y[i].clone());
                ct.modpow(&y, &self.params.p).modinv(&self.params.p).expect("Inverse does not exist")
            } else {
                ct.modpow(&y[i], &self.params.p)
            }; 
            num = (&num * &t1).mod_floor(&self.params.p.clone());
        }

        let denom:BigInt = if key < zero {
            let x  = -(key.clone());
            (ciphertext[0]).modpow(&x, &self.params.p).modinv(&self.params.p).expect("Inverse does not exist")
        } else {
            (ciphertext[0]).modpow(&key, &self.params.p)
        };
        let denom_inv = denom.modinv(&self.params.p).expect("Inverse does not exist");
        let r = (num * denom_inv).mod_floor(&self.params.p);
        
        let c = in_zp(self.params.p.clone(), self.params.q.clone());
        let decrypted = baby_step_giant_step(c, r, self.params.g.clone()).expect("Failed to find discrete logarithm");
        decrypted.into_pyobject(py).unwrap()
    }
}

#[pymodule]
fn ddhlib(_py: Python<'_>, m: Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Ddh>()?;
    m.add_class::<UniformRange>()?;
    Ok(())
}