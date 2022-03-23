use std::error::Error;
use std::io::{BufRead, BufReader};
use std::{self,fs};

const pi            : f64 = 3.141592653589793238;
const m_N_c2        : f64 = 939.1 * 1.602e-13;
const m_N           : f64 = 1.6715e027;
const k_B           : f64 = 1.380649e-23;
const hbar          : f64 = 1.054571e-34;
const hbar_c        : f64 = hbar * 2.998e8;
const mu            : f64 = 1.660539e-27;
const NA            : f64 = 6.02214076e26;



struct ElementProfile {
    Element: String,
    Z: f64,
    A: f64,
    N: f64,
    Mass: f64,
    SolarMF: f64,
    MF: f64,
    IMF: f64,
}

impl ElementProfile {
    fn new() -> Self {
        Self {Element: "".to_string(), Z: 0., A: 0., N: 0., Mass: 0., SolarMF: 0., MF: 0., IMF: 0.,}
    }
}

struct Condition {
    Ye: f64,
    Temperature: f64,
    Density: f64,
}    

impl Condition {
    fn Ye_str (&self) -> String {
        let mut s = self.Ye.to_string();
        s.remove(1);
        s
    }    
}    

enum Status {
    NSE,
    Freezeout,
    Last,
    Information,
}

struct ElementData {
    kind: Status,
    profile: Vec<ElementProfile>,
}

impl ElementData {
    fn calc_entropy (&self, condition: &Condition) -> Entropy {
        let T = condition.Temperature;
        let rho = condition.Density;
        let Ye = condition.Ye;
        
        let S_rad = 11. * pi * pi * T * T * T / (45. * rho);

        let mu_e = (3. * pi * pi * Ye * rho).powf(1./3.)
        let S_deg = mu_e.powf(2.) * T / (3 * rho)

        let S_ideal = 5./(2. * m_N) + (2 * (2 * pi * m_N * k_B * T).powf(3./2.) / )
        
        let data = Elementdata.profile;
        Entropy{sum: 0., deg: 0., ideal: 0., rad: 0.}
    }

    fn is_sumMF_less1 (&self) -> bool{
        false
    }

    fn is_cnt_less100 (&self) -> bool{
        true
    }
}

struct Entropy {
    sum: f64,
    deg: f64,
    ideal: f64,
    rad: f64,
}

impl Entropy {
    fn new() -> Self {
        Self {sum: 0., deg: 0., ideal: 0., rad: 0.}
    }

    fn calc(condition: &Condition, element_data: &Vec<ElementProfile>) -> Self {
        let A: f64 = k_B.powf(3.) * 7.0 * pi * pi * m_N / (hbar_c.powf(3.) * 45.);
        let B: f64 = (3. * pi.powf(2.)).powf(2./3.) * k_B * m_N.powf(1./3.) / (3. * hbar_c);
        let Ye = condition.Ye;
        let T = condition.Temperature;
        let rho = condition.Density;

        let S_rad: f64 = A * T.powf(3.) / rho;
        let S_dec: f64 = B * Ye.powf(2./3.) * T / rho.powf(1./3.);

        let mut S_ideal: f64 = 0.;

        for data in element_data{
            // println!("{:?}", data.Element);
            let Xi = data.MF;
            let mi = data.Mass * mu;
            let mi_c2 = mi * 9e16;
            let Yi = Xi / (mi * NA);
            let ni = Xi * rho / mi;
            let LOG_NAKAMI = (1./ni) * ((mi_c2 * k_B * T) / (2. * pi * hbar_c.powf(2.))).powf(1.5);

            let Sci = if Xi == 0. || Xi < 1e-300 {0.} else {Yi * (5./2. + LOG_NAKAMI.log10())};
            S_ideal += Sci;
        }
        Entropy{rad: S_rad, deg: S_dec, ideal: S_ideal, sum: S_rad + S_dec + S_ideal}
    }
}

fn main() {
    let status = ["NSE", "freezeout", "last", "Information"];
    
    let Ye0 = [0.01, 0.03, 0.09, 0.10, 0.13, 0.16, 0.19, 0.20, 0.23, 0.26, 0.29, 0.30, 0.36, 0.39, 0.40, 0.43, 0.46, 0.49, 0.50];
    let T0 = [4e9, 7e9, 1e10, 4e10, 7e10, 1e11, 4e11, 7e11, 1e12];
    let rho0 = [1e10, 4e10, 7e10, 1e11, 4e11, 7e11, 1e12, 4e12, 7e12, 1e13, 4e13];
    
    for Ye0 in Ye0.iter() {
        for T0 in T0.iter() {
            for rho0 in rho0.iter() {
                for status in status.iter() {
                    // 各状況でconditionという名前をつける。
                    let condition = Condition{Ye: *Ye0, Temperature: *T0, Density: *rho0};

                    // このconditionでのデータを読み込み。
                    let data = load(&condition, status.to_string());
                    
                    // 総MFを計算。

                    if !data.is_sumMF_less1() {
                        panic!("sumMF is larger than 1");
                    }

                    let entropy = data.calc_entropy(&condition);
                    println!(
                        "S_rad = {}, S_deg = {}, S_ideal = {}, S_sum = {}", 
                        entropy.rad, entropy.deg, entropy.ideal, entropy.sum
                    )
                    // calc_MF_sum(&data);

                    // InformationDataを読み込む。
                    // load_Infomartion()

                    // エントロピーを計算。
                    // calc_entropy(&data);

                    // 図１４を再現するのに、重み付きの和を取る。
                    return
                }
            }
        }
    }
}

fn load (condition: &Condition, status: String) -> ElementData{
    // 読み込むファイルの名前とパスを指定。
    let input_data = format!("../../../10228to0307/original/Ye{1}/data/{0}_Ye_{1}_T0_{2:e}_rho0_{3:e}.txt",
    status, condition.Ye_str(), condition.Temperature, condition.Density
    );

    // load_data関数で得られるVec<Element>型データを取り出して、返り値にする。
    load_data(&input_data, status).unwrap()
}


fn load_data (filename: &str, status: String) -> Result<ElementData, Box<std::error::Error>> {
    // 返り値となるVec<Element>の型を定義。
    let mut ret: Vec<ElementProfile> = Vec::new();
    
    // BufReaderを用いてファイルを開き、各行をStringとして、Vec<String>に落とし込む。
    let lines: Vec<String> = BufReader::new(fs::File::open(filename)?)
                            .lines()
                            .map(|line| line.expect("Failed to read from BufReader"))
                            .collect();
                            println!("load_data");
    // println!("{:?}", lines);

    // 行ごとに処理していく。
    for (i, line) in lines.iter().enumerate() {
        // 空白で区切り、Vecとして保持。
        let line = line.trim().split_whitespace().collect::<Vec<_>>();
        // println!("a: {:?}", line);

        // Vec<String>型として取り出されたある行のデータをElementData型のタプルとする。
        let now = ElementProfile {
            Element     : line[0].to_string(),
            Z           : line[1].parse().unwrap(),
            A           : line[2].parse().unwrap(),
            N           : line[3].parse().unwrap(),
            Mass        : line[4].parse().unwrap(),
            SolarMF     : line[5].parse().unwrap(),
            MF          : line[6].parse().unwrap(),
            IMF         : line[7].parse().unwrap(),
        };

        // 返り値のVecにpushする。
        ret.push(now);
    }
    let kind = match &*status {
        "NSE" => Status::NSE,
        "freezeout" => Status::Freezeout,
        "last" => Status::Last,
        "Information" => Status::Information,
        _ => unreachable!(),
    };
    //読み終えたら、Result型でreturnする。

    Ok(ElementData{kind: kind, profile: ret})
}