//! Collection of da/dn equations for calculating crack growth.
//!
//! So far these subroutines are history independent and only depend
//! on $\Delta K$ and $R$.

use crate::table;
use log::{error, info, warn};
use std::collections::BTreeMap;
use std::f64::consts::FRAC_PI_2;
use std::{f64, fmt, process};

/// Shared settings
#[derive(Debug, Clone)]
pub struct Options {
    pub rmax: f64,
    pub rmin: f64,
    pub deltak_th: f64,
    pub kneg: bool,
    pub kmax_th: f64,
    pub kmax_th_follows_deltak_th: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            rmax: 0.99,
            rmin: -1.0,
            deltak_th: 0.0,
            kneg: false,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        }
    }
}

/// Defines the state of the crack for any dadn equation that requires
/// some sort of memory.
pub struct CrackState {
    // length of the crack
    pub a: f64,
}

/// A standard definition of variable names used in dadn equations.
///
/// The rust style is being intentionally disregarded here, as
/// case may matter in some equations.
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ParameterLabel {
    a,
    b,
    c,
    d,
    e,
    f,
    m,
    n,
    p,
    q,
    a_intr,
    alpha,
    c1,
    c2,
    cth,
    cth_minus,
    deltak_th,
    deltak0,
    k_crit,
    k_ut,
    kic,
    kf,
    smax_on_sigma0,
}

impl ParameterLabel {
    /// Return a printable value for the label.
    // When adding new entries, it is important to note that this
    // text value will be presented to the user as required input
    // into the --parameters cli argument.
    // Ensure that it is an easily typed string, to maintain a
    // good user experience.
    fn text(&self) -> &'static str {
        match *self {
            ParameterLabel::a => "a",
            ParameterLabel::b => "b",
            ParameterLabel::c => "c",
            ParameterLabel::d => "d",
            ParameterLabel::e => "e",
            ParameterLabel::f => "f",
            ParameterLabel::m => "m",
            ParameterLabel::n => "n",
            ParameterLabel::p => "p",
            ParameterLabel::q => "q",
            ParameterLabel::a_intr => "a_intr",
            ParameterLabel::alpha => "alpha",
            ParameterLabel::c1 => "c1",
            ParameterLabel::c2 => "c2",
            ParameterLabel::cth => "cth",
            ParameterLabel::cth_minus => "cth_minus",
            ParameterLabel::deltak_th => "deltak_th",
            ParameterLabel::deltak0 => "deltak0",
            ParameterLabel::k_crit => "k_crit",
            ParameterLabel::k_ut => "k_ut",
            ParameterLabel::kic => "kic",
            ParameterLabel::kf => "kf",
            ParameterLabel::smax_on_sigma0 => "smax_on_sigma0",
        }
    }

    pub fn from_text(input: &str) -> Option<ParameterLabel> {
        match input {
            "a" => Some(ParameterLabel::a),
            "b" => Some(ParameterLabel::b),
            "c" => Some(ParameterLabel::c),
            "d" => Some(ParameterLabel::d),
            "e" => Some(ParameterLabel::e),
            "f" => Some(ParameterLabel::f),
            "m" => Some(ParameterLabel::m),
            "n" => Some(ParameterLabel::n),
            "p" => Some(ParameterLabel::p),
            "q" => Some(ParameterLabel::q),
            "a_intr" => Some(ParameterLabel::a_intr),
            "alpha" => Some(ParameterLabel::alpha),
            "c1" => Some(ParameterLabel::c1),
            "c2" => Some(ParameterLabel::c2),
            "cth" => Some(ParameterLabel::cth),
            "cth_minus" => Some(ParameterLabel::cth_minus),
            "deltak_th" => Some(ParameterLabel::deltak_th),
            "deltak0" => Some(ParameterLabel::deltak0),
            "k_crit" => Some(ParameterLabel::k_crit),
            "k_ut" => Some(ParameterLabel::k_ut),
            "kic" => Some(ParameterLabel::kic),
            "kf" => Some(ParameterLabel::kf),
            "smax_on_sigma0" => Some(ParameterLabel::smax_on_sigma0),
            _ => None,
        }
    }
}

// pub enum DadnEqn {
//     Nasgro,
//     Forman,
//     Paris,
//     Walker,
//     Burchill,
//     Hartman,
//     White,
//     Kujawski,
// }

/// data for White equation
/// $$ dadn =
#[derive(Debug, Clone)]
pub struct White {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    kic: f64,
    options: Options,
}

/// data for Kujawski equation
/// $$ dadn =
#[derive(Debug, Clone)]
pub struct Kujawski {
    c: f64,
    m: f64,
    alpha: f64,
    options: Options,
}

// Burchill equation
// dadn = C Kmax^m - D Kmin^n
#[derive(Debug, Clone)]
pub struct Burchill {
    c: f64,
    m: f64,
    d: f64,
    n: f64,
    options: Options,
}

/// data for Paris equation
/// $$ dadn = C \Delta K^m $$
#[derive(Debug, Clone)]
pub struct Paris {
    c: f64,
    m: f64,
    options: Options,
}

/// Forman equation
/// $$ da/dn = C \Delta K^n / ((1-R) * `K_f`  - \Delta K) $$
#[derive(Debug, Clone)]
pub struct Forman {
    c: f64,
    n: f64,
    kf: f64,
    options: Options,
}

/// Forman equation - Modified to behave the same as in AFGROW
/// $$ da/dn = C \Delta K^n / ((1-R) * `K_f`  - \Delta K) $$
#[derive(Debug, Clone)]
pub struct FormanAFGROW {
    c: f64,
    n: f64,
    kf: f64,
    options: Options,
}

/// Walker Equation
/// $$ dadn = C \Delta K^m  $$
#[derive(Debug, Clone)]
pub struct Walker {
    c: f64,
    m: f64,
    n: f64,
    options: Options,
}

/// Hartman-Schijve Variant (Jones-Molent)
/// The influence of cyclic stress intensity threshold on fatigue life scatter
/// L. Molent and R. Jones
/// International Journal of Fatigue
///Volume 82, Part 3, January 2016, Pages 748–756
#[derive(Debug, Clone)]
pub struct Hartman {
    d: f64,
    a: f64,
    alpha: f64,
    options: Options,
}

#[derive(Debug, Clone)]
pub struct Nasgro {
    // ratio of maximum far field stress to flow stress sigma0
    smax_on_sigma0: f64,
    // constraint factor
    alpha: f64,
    // fracture toughness
    k_crit: f64,
    // delta k threshold at R=0
    deltak0: f64,
    // curve control coefficient for different values of R, equals 0 for negative R, equals 1 for R>=0.
    cth: f64,
    p: f64,
    // emperical crack coefficient
    q: f64,
    // emperical crack coefficient
    c: f64,
    // emperical crack coefficient
    n: f64,
    // intrinsic crack size (m)
    a_intr: f64,
    options: Options,
}

#[derive(Clone)]
pub struct Tabular {
    table: table::Table,
    options: Options,
}

/// Linear Lower Threshold 5 equation
/// b (\Delta K - \Delta K_th)^a (1 - R)^c
#[derive(Debug, Clone)]
pub struct LinLowerTh5 {
    a: f64,
    b: f64,
    c: f64,
    options: Options,
}

/// Linear Lower Threshold 5A equation
/// b (\Delta K - \Delta K_th)^a (1 - R)^c (1 - (\Delta K / (K_ut (1 - R))))^d
#[derive(Debug, Clone)]
pub struct LinLowerTh5A {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    k_ut: f64,
    options: Options,
}

/// Linear Lower Threshold 5B equation
/// 
/// For R >= 0:
/// b (\Delta K - \Delta K_th)^a (1 - R)^c1 (1 - (Kmax / K_ut))^d
/// 
/// For R < 0:
/// b (\Delta K - \Delta K_th)^a (1 - R)^c2 (1 - (Kmax / K_ut))^d
#[derive(Debug, Clone)]
pub struct LinLowerTh5B {
    a: f64,
    b: f64,
    c1: f64,
    c2: f64,
    d: f64,
    k_ut: f64,
    options: Options,
}

/// Linear Lower Threshold 6 equation
/// b (\Delta K - \Delta K_th)^a (1 - R)^(c + d * \Delta K)
#[derive(Debug, Clone)]
pub struct LinLowerTh6 {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    options: Options,
}

/// Linear Lower Threshold 6A equation
/// b (\Delta K - \Delta K_th)^a (1 - R)^(c + d * \Delta K) (1 - (\Delta K / (K_ut (1 - R))))^e
#[derive(Debug, Clone)]
pub struct LinLowerTh6A {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    k_ut: f64,
    options: Options,
}

/// No R-effects Lower Threshold 1 equation
/// 
/// For use with the picc growth method
/// 
/// b (\Delta K - \Delta K_th)^a
#[derive(Debug, Clone)]
pub struct NoRLowerTh1 {
    a: f64,
    b: f64,
    options: Options,
}

/// No R-effects Lower Threshold 2 equation
/// 
/// For use with the picc growth method
/// 
/// b (\Delta K)^a (1 - (\Delta K_th / \Delta K)^c)
#[derive(Debug, Clone)]
pub struct NoRLowerTh2 {
    a: f64,
    b: f64,
    c: f64,
    options: Options,
}

/// No R-effects Lower and Upper Threshold 1 equation
/// 
/// For use with the picc growth method
/// 
/// b (\Delta K - \Delta K_th)^a (1 - (Kmax / K_ut))^c
#[derive(Debug, Clone)]
pub struct NoRLowerUpperTh1 {
    a: f64,
    b: f64,
    c: f64,
    k_ut: f64,
    options: Options,
}

/// No R-effects Lower and Upper Threshold 2 equation
/// 
/// For use with the picc growth method
/// 
/// b (\Delta K)^a (1 - (\Delta K_th / \Delta K)^c) (1 / (1 - (Kmax / K_ut)^d))
#[derive(Debug, Clone)]
pub struct NoRLowerUpperTh2 {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    k_ut: f64,
    options: Options,
}



/// Trait for `DaDn` function that can be used for crack growth
pub trait DaDn {
    /// Display the crack growth equation
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result;

    /// dadn equation gives the instaneous rate of crack growth
    ///
    /// Currently we also pass it the crack length which is currently
    /// not used by all equations but it could be generalise to be a
    /// general memory state parameter.
    fn dadn(&self, kmin: f64, kmax: f64, state: CrackState) -> f64;

    /// Update the equation parameters
    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>);

    fn get_name(&self) -> &str;

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync>;
}

impl Clone for Box<dyn DaDn + Send + Sync> {
    fn clone(&self) -> Self {
        self.inner_clone()
    }
}

impl fmt::Display for dyn DaDn + Send + Sync {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt(f)
    }
}

impl fmt::Display for ParameterLabel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.text())
    }
}

/// Apply limits and calculate r and delta_k values.
/// Returns a tuple of the form: (r, delta_k)
fn get_r_deltak(mut kmin: f64, kmax: f64, rmin: f64, rmax: f64, kneg: bool) -> (f64, f64) {
    let r = if kmin / kmax < rmin {
        kmin = kmax * rmin;
        rmin
    } else {
        (kmin / kmax).min(rmax)
    };

    let delta_k = if kmin >= 0.0 || kneg {
        kmax - kmin
    } else {
        kmax
    };

    (r, delta_k)
}

fn update_kmax_th(options: &mut Options, deltak_th: Option<f64>) {
    if options.kmax_th_follows_deltak_th && deltak_th.is_some() {
        options.kmax_th = deltak_th.unwrap();
    }
}

/// White equation.
/// An extension of the Forman equation to create a dadn curve made from
/// a combination of a cubic with an asymptotic curve for kic effect.
/// Used in the specification for MSMP3.

impl White {
    pub const NAME: &'static str = "white";
    const CITE: &'static str = "unknown";
    pub const UNITS: &'static str = "m";
    // A list of the parameters/variables/coefficients used in the equation
    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::d,
        ParameterLabel::e,
        ParameterLabel::f,
        ParameterLabel::kic,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> White {
        // Check first that we have all the required parameters

        White {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            d: params[&ParameterLabel::d],
            e: params[&ParameterLabel::e],
            f: params[&ParameterLabel::f],
            kic: params[&ParameterLabel::kic],
            options,
        }
    }
}

impl DaDn for White {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:#?}", self);
        write!(
            f,
            "#  da/dN ({units}) = exp[({a:e} * ΔKeff^3 - {b:e} * ΔKeff^2 + {c:e} * ΔKeff - {d:e}) 
#              + (dkic - ΔK)^{f}] [White14:{cite}]
#  where dkic = {kic:e} (1 - R) and ΔKeff = ΔK / (1 - R)^{e}",
            units = White::UNITS,
            a = self.a,
            b = self.b,
            c = self.c,
            d = self.d,
            e = self.e,
            kic = self.kic,
            f = -self.f,
            cite = White::CITE
        )
    }

    // note the signs of the b,d and f variables has been change to make all the coefficients positive
    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        // this function does not work well if the rmax variable is
        // included in the optimisation so we set it directly
        // let rmax = coeffs[7];
        let rmax = 0.8; // hardwired

        let delta_k = kmax - kmin.max(0.0);
        let r = (kmin.max(0.0) / kmax).min(rmax);

        let dkic = self.kic * (1.0 - r);
        let d_ke = (delta_k / (1.0 - r).powf(self.e)).ln();

        // simple cubic curve
        let main = self.a * d_ke.powi(3) - self.b * d_ke.powi(2) + self.c * d_ke - self.d;

        // extra bit due to approaching fracture toughness
        let bit = (dkic - delta_k).max(0.01).powf(-self.f);

        if delta_k > 0.0 {
            (main + bit).exp()
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.d = params[&ParameterLabel::d];
        self.e = params[&ParameterLabel::e];
        self.f = params[&ParameterLabel::f];
        self.kic = params[&ParameterLabel::kic];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        White::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Kujawski equation
/// parameters = [c, n, k]
impl Kujawski {
    pub const NAME: &'static str = "kujawski";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] =
        &[ParameterLabel::c, ParameterLabel::m, ParameterLabel::alpha];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Kujawski {
        Kujawski {
            c: params[&ParameterLabel::c],
            m: params[&ParameterLabel::m],
            alpha: params[&ParameterLabel::alpha],
            options,
        }
    }
}

impl DaDn for Kujawski {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * [Kmax ^ {alpha} * ΔK^ {{1 - {alpha}}}]^{m}  [Kujawski01]",
            units = Kujawski::UNITS,
            c = self.c,
            m = self.m,
            alpha = self.alpha,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = kmax.max(0.0) - kmin.max(0.0);
        self.c * (kmax.max(0.0).powf(self.alpha) * delta_k.powf(1.0 - self.alpha)).powf(self.m)
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.c = params[&ParameterLabel::c];
        self.m = params[&ParameterLabel::m];
        self.alpha = params[&ParameterLabel::alpha];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Kujawski::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Burchill equation
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
impl Burchill {
    pub const NAME: &'static str = "burchill";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::c,
        ParameterLabel::m,
        ParameterLabel::d,
        ParameterLabel::n,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Burchill {
        Burchill {
            c: params[&ParameterLabel::c],
            m: params[&ParameterLabel::m],
            d: params[&ParameterLabel::d],
            n: params[&ParameterLabel::n],
            options,
        }
    }
}

impl DaDn for Burchill {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * Kmax ^ {m} - {d:e} Kmin ^ {n}]  [Burchill17]",
            units = Burchill::UNITS,
            c = self.c,
            m = self.m,
            d = self.d,
            n = self.n
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        self.c * kmax.max(0.0).powf(self.m) - self.d * kmin.max(0.0).powf(self.n)
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.c = params[&ParameterLabel::c];
        self.m = params[&ParameterLabel::m];
        self.d = params[&ParameterLabel::d];
        self.n = params[&ParameterLabel::n];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Burchill::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Nasgro  equation
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
/// Ref: AFGROW users guide and technical manual
/// James A. Harter
/// AFRL-VA-WP-TR-1999-3016
/// Feb 1999
impl Nasgro {
    pub const NAME: &'static str = "nasgro";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::smax_on_sigma0,
        ParameterLabel::alpha,
        ParameterLabel::k_crit,
        ParameterLabel::deltak0,
        ParameterLabel::cth,
        ParameterLabel::p,
        ParameterLabel::q,
        ParameterLabel::c,
        ParameterLabel::n,
        ParameterLabel::a_intr,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Nasgro {
        Nasgro {
            smax_on_sigma0: params[&ParameterLabel::smax_on_sigma0],
            alpha: params[&ParameterLabel::alpha],
            k_crit: params[&ParameterLabel::k_crit],
            deltak0: params[&ParameterLabel::deltak0],
            cth: params[&ParameterLabel::cth],
            p: params[&ParameterLabel::p],
            q: params[&ParameterLabel::q],
            c: params[&ParameterLabel::c],
            n: params[&ParameterLabel::n],
            a_intr: params[&ParameterLabel::a_intr],
            options,
        }
    }
}

// According to the AFGROW manual, the variable
// Smax_on_simga is the ratio of the maximum applied (far field) stress
// to the flow stress (of the material).  This does not make a lot of
// sense as they really have nothing to do with each other.
impl DaDn for Nasgro {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        _e = writeln!(
            f,
            "#  da/dN ({units}) = {c:e} * (((1 - f) / (1 - R)) ΔK ^ {n}) * G  [Nasgro99]",
            units = Nasgro::UNITS,
            c = self.c,
            n = self.n
        );
        writeln!(
            f,
            "#  G = (1 - (ΔK_th / ΔK)) ^ {p} /  (1 - (Kmax / Kcrit)) ^ {q}   ",
            p = self.p,
            q = self.q
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let r = kmin / kmax;
        let delta_k = kmax - kmin;

        // closure constants
        let a0 = (0.825 - 0.34 * self.alpha + 0.05 * self.alpha.powi(2))
            * (FRAC_PI_2 * self.smax_on_sigma0)
                .cos()
                .powf(1.0 / self.alpha);
        let a1 = (0.415 - 0.071 * self.alpha) * self.smax_on_sigma0;
        let a3 = 2.0 * a0 + a1 - 1.0;
        let a2 = 1.0 - a0 - a1 - a3;

        // closure level f
        let f = if r >= 0.0 {
            (a0 + a1 * r + a2 * r.powi(2) + a3 * r.powi(3)).max(r)
        } else if (-2.0..0.0).contains(&r) {
            a0 + a1 * r
        } else {
            a0 - 2.0 * a1
        };

        info!("Nasgro: {} {} {} {} {}", a0, a1, a3, a2, f);
        // a_int = intrinsic crack size, typically 38.1e-6 (m)
        let deltak_th = self.deltak0 * (state.a / (state.a + self.a_intr)).sqrt()
            / ((1.0 - f) / ((1.0 - a0) * (1.0 - r))).powf(1.0 + self.cth * r);

        info!("nasgro: deltak_th {}", deltak_th);
        let num = (1.0 - (deltak_th.min(delta_k - 1e-6) / delta_k)).powf(self.p);
        let denom = (1.0 - (kmax / self.k_crit)).powf(self.q);
        let dadn = self.c * (((1.0 - f) / (1.0 - r)) * delta_k).powf(self.n) * num / denom;
        info!("nasgro: dadn {} {} {}", dadn, num, denom);
        dadn
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.smax_on_sigma0 = params[&ParameterLabel::smax_on_sigma0];
        self.alpha = params[&ParameterLabel::alpha];
        self.k_crit = params[&ParameterLabel::k_crit];
        self.deltak0 = params[&ParameterLabel::deltak0];
        self.cth = params[&ParameterLabel::cth];
        self.p = params[&ParameterLabel::p];
        self.q = params[&ParameterLabel::q];
        self.c = params[&ParameterLabel::c];
        self.n = params[&ParameterLabel::n];
        self.a_intr = params[&ParameterLabel::a_intr];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Nasgro::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Forman equation
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
impl Forman {
    pub const NAME: &'static str = "forman";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::c,
        ParameterLabel::n,
        ParameterLabel::kf,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Forman {

        Forman {
            c: params[&ParameterLabel::c],
            n: params[&ParameterLabel::n],
            kf: params[&ParameterLabel::kf],
            options,
        }
    }
}

impl DaDn for Forman {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * ΔK ^ {n} / [(1 - R) {kf} - ΔK]  [Forman67]",
            units = Forman::UNITS,
            c = self.c,
            n = self.n,
            kf = self.kf
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        if delta_k > self.options.deltak_th {
            self.c * delta_k.powf(self.n) / ((1.0 - r) * self.kf - delta_k)
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.c = params[&ParameterLabel::c];
        self.n = params[&ParameterLabel::n];
        self.kf = params[&ParameterLabel::kf];
        if let Some(value) = deltak_th {
            self.options.deltak_th = value;
        }
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Forman::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Forman equation - Modified to behave the same as in AFGROW
/// da/dN = c (\Delta K)^n / ((1-R) Kf - \Delta K)
/// parameters = [c, n, k]
impl FormanAFGROW {
    pub const NAME: &'static str = "forman-ag";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::c,
        ParameterLabel::n,
        ParameterLabel::kf,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Self {
        Self {
            c: params[&ParameterLabel::c],
            n: params[&ParameterLabel::n],
            kf: params[&ParameterLabel::kf],
            options,
        }
    }
}

impl DaDn for FormanAFGROW {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * ΔK ^ {n} / [(1 - R) {kf} - ΔK]  [Forman67]",
            units = FormanAFGROW::UNITS,
            c = self.c,
            n = self.n,
            kf = self.kf
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        let result = self.c * delta_k.powf(self.n) / ((1.0 - r) * self.kf - delta_k);
        
        let dadn_th = self.c * self.options.deltak_th.powf(self.n) / (self.kf - self.options.deltak_th);

        if result <= dadn_th {0.0} else {result}
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.c = params[&ParameterLabel::c];
        self.n = params[&ParameterLabel::n];
        self.kf = params[&ParameterLabel::kf];
        if let Some(value) = deltak_th {
            self.options.deltak_th = value;
        }
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        FormanAFGROW::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Paris equation
/// da/dN = c (\Delta K)^m
/// parameters = [c, m]
impl Paris {
    pub const NAME: &'static str = "paris";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[ParameterLabel::c, ParameterLabel::m];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Paris {
        Paris {
            c: params[&ParameterLabel::c],
            m: params[&ParameterLabel::m],
            options,
        }
    }
}

impl DaDn for Paris {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({}) = {:e} * ΔK ^ {} [Paris63]",
            Paris::UNITS,
            self.c,
            self.m
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = kmax.max(0.0) - kmin.max(0.0);

        self.c * delta_k.powf(self.m)
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.c = params[&ParameterLabel::c];
        self.m = params[&ParameterLabel::m];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Paris::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Hartman-Schijve equation
/// c [(\Delta K - kth) / \sqrt{1 - (kmax/a)}]^n
/// parameters = [c, n, kth, a]
impl Hartman {
    pub const NAME: &'static str = "hartman";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::d,
        ParameterLabel::deltak_th,
        ParameterLabel::a,
        ParameterLabel::alpha,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> Hartman {
        options.deltak_th = params[&ParameterLabel::deltak_th];

        Hartman {
            d: params[&ParameterLabel::d],
            a: params[&ParameterLabel::a],
            alpha: params[&ParameterLabel::alpha],
            options,
        }
    }
}

impl DaDn for Hartman {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(f, "#  da/dN ({units}) =  {d:e} [(ΔK - {deltak_th}) / sqrt (1 - kmax/{a})]^{alpha} [Hartman70]",
               units=Hartman::UNITS, d=self.d, deltak_th=self.options.deltak_th, a=self.a, alpha=self.alpha)
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = kmax - kmin.max(0.0);
        let r = kmin.max(0.0) / kmax;

        let kmax = delta_k / (1.0 - r);

        if kmax > self.a {
            warn!("***Warning: A kmax of {} is > {} and therefore cannot be square-rooted in the Hartman-Schijve equation", kmax, self.a);
        }

        if delta_k > self.options.deltak_th {
            self.d
                * ((delta_k - self.options.deltak_th).max(0.0) / (1.0 - (kmax / self.a)).max(0.0).sqrt())
                    .powf(self.alpha)
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.d = params[&ParameterLabel::d];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        self.a = params[&ParameterLabel::a];
        self.alpha = params[&ParameterLabel::alpha];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Hartman::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Walker equation
///
/// da/dn = c * (\Delta K * (1 - R)^(m - 1))^n
/// parameters = [C, n, m]
impl Walker {
    pub const NAME: &'static str = "walker";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] =
        &[ParameterLabel::c, ParameterLabel::m, ParameterLabel::n];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, options: Options) -> Walker {
        Walker {
            c: params[&ParameterLabel::c],
            m: params[&ParameterLabel::m],
            n: params[&ParameterLabel::n],
            options,
        }
    }
}

/// In this implementation R is not limited to a minimium of 0.0 .
impl DaDn for Walker {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {c:e} * [ΔK / (1-R)^(1 - {m})] ^ {n} [Walker70]",
            units = Walker::UNITS,
            c = self.c,
            m = self.m,
            n = self.n
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = (kmax - kmin.max(0.0)).max(0.0);
        let r = (kmin / kmax).clamp(0.0, 0.99);

        self.c * (delta_k / (1.0 - r).powf(1.0 - self.m)).powf(self.n)
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.c = params[&ParameterLabel::c];
        self.m = params[&ParameterLabel::m];
        self.n = params[&ParameterLabel::n];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        Walker::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl Tabular {
    pub const NAME: &'static str = "file";

    pub fn new(file_name: &str, options: Options) -> Self {
        Self {
            table: make_table_from_file(file_name, &[]),
            options,
        }
    }
}

impl DaDn for Tabular {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let _ = write!(
            f,
            "  Spline interpolated tabular lookup
#       da/dn:"
        );

        for col in &self.table.columns {
            let _ = write!(f, " {:12}", col);
        }
        let _ = writeln!(f);
        for i in 0..self.table.row.len() {
            let _ = write!(f, "#  {:10.3e}: ", (10.0f64).powf(self.table.row[i]));
            for j in 0..self.table.values.len() {
                let _ = write!(f, "  {:10} ", self.table.values[j][i]);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0;
        }

        info!("kmin {} kmax {}", kmin, kmax);
        let rmax = self.table.columns[self.table.columns.len() - 1].min(self.options.rmax);
        let rmin = self.table.columns[0].max(self.options.rmin);
        let (r, delta_k) = get_r_deltak(kmin, kmax, rmin, rmax, self.options.kneg);

        info!("delta k {dk} r limited {r}", dk = delta_k, r = r);
        let interp = self.table.interp(delta_k, r);
        info!("Interp {}", interp);
        let value = (10.0f64).powf(interp);
        info!("value {}", value);
        value
    }

    fn update_parameters(&mut self, _params: &BTreeMap<ParameterLabel, f64>, _deltak_th: Option<f64>) {
        // Do nothing. Future implementations may optimise tabular data
    }

    fn get_name(&self) -> &str {
        Tabular::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl LinLowerTh5 {
    pub const NAME: &'static str = "lin_lower_th_5";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::deltak_th,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> LinLowerTh5 {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        LinLowerTh5 {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            options,
        }
    }
}

impl DaDn for LinLowerTh5 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a} (1 - R)^{c}",
            units = LinLowerTh5::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            c = self.c,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a) * (1.0 - r).powf(self.c)
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        LinLowerTh5::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl LinLowerTh5A {
    pub const NAME: &'static str = "lin_lower_th_5a";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::d,
        ParameterLabel::deltak_th,
        ParameterLabel::k_ut,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> LinLowerTh5A {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        LinLowerTh5A {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            d: params[&ParameterLabel::d],
            k_ut: params[&ParameterLabel::k_ut],
            options,
        }
    }
}

impl DaDn for LinLowerTh5A {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a} (1 - R)^{c} (1 - (Kmax / {k_ut}))^{d}",
            units = LinLowerTh5A::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            c = self.c,
            k_ut = self.k_ut,
            d = self.d,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a) * 
                (1.0 - r).powf(self.c) * (1.0 - (kmax / self.k_ut )).powf(self.d)
        } else {
            0.0
        }

    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.d = params[&ParameterLabel::d];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        self.k_ut = params[&ParameterLabel::k_ut];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        LinLowerTh5A::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl LinLowerTh5B {
    pub const NAME: &'static str = "lin_lower_th_5b";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::d,
        ParameterLabel::c1,
        ParameterLabel::c2,
        ParameterLabel::deltak_th,
        ParameterLabel::k_ut,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> LinLowerTh5B {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        LinLowerTh5B {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c1: params[&ParameterLabel::c1],
            c2: params[&ParameterLabel::c2],
            d: params[&ParameterLabel::d],
            k_ut: params[&ParameterLabel::k_ut],
            options,
        }
    }
}

impl DaDn for LinLowerTh5B {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a} (1 - R)^(c1: {c1}, c2: {c2}) (1 - (Kmax / {k_ut}))^{d}",
            units = LinLowerTh5B::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            c1 = self.c1,
            c2 = self.c2,
            k_ut = self.k_ut,
            d = self.d,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        let c = if r >= 0.0 {self.c1} else {self.c2};

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a) * 
                (1.0 - r).powf(c) * (1.0 - (kmax / self.k_ut )).powf(self.d)
        } else {
            0.0
        }

    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c1 = params[&ParameterLabel::c1];
        self.c2 = params[&ParameterLabel::c2];
        self.d = params[&ParameterLabel::d];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        self.k_ut = params[&ParameterLabel::k_ut];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        LinLowerTh5B::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl LinLowerTh6 {
    pub const NAME: &'static str = "lin_lower_th_6";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::d,
        ParameterLabel::deltak_th,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> LinLowerTh6 {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        LinLowerTh6 {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            d: params[&ParameterLabel::d],
            options,
        }
    }
}

impl DaDn for LinLowerTh6 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a} (1 - R)^({c} + {d} ΔK)",
            units = LinLowerTh6::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            c = self.c,
            d = self.d,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a) *
                (1.0 - r).powf(self.c + self.d * delta_k) 
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.d = params[&ParameterLabel::d];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        LinLowerTh6::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl LinLowerTh6A {
    pub const NAME: &'static str = "lin_lower_th_6a";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::d,
        ParameterLabel::e,
        ParameterLabel::deltak_th,
        ParameterLabel::k_ut,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> LinLowerTh6A {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        LinLowerTh6A {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            d: params[&ParameterLabel::d],
            e: params[&ParameterLabel::e],
            k_ut: params[&ParameterLabel::k_ut],
            options,
        }
    }
}

impl DaDn for LinLowerTh6A {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a} (1 - R)^({c} + {d} ΔK) (1 - (ΔK / {k_ut}))^{e}",
            units = LinLowerTh6A::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            c = self.c,
            d = self.d,
            k_ut = self.k_ut,
            e = self.e,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let (r, delta_k) = get_r_deltak(kmin, kmax, self.options.rmin, self.options.rmax, self.options.kneg);

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a) 
                * (1.0 - r).powf(self.c + self.d * delta_k) 
                * (1.0 - (delta_k / self.k_ut )).powf(self.e)
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.d = params[&ParameterLabel::d];
        self.e = params[&ParameterLabel::e];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        self.k_ut = params[&ParameterLabel::k_ut];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        LinLowerTh6A::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl NoRLowerTh1 {
    pub const NAME: &'static str = "no_r_low_th1";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::deltak_th,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> NoRLowerTh1 {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        NoRLowerTh1 {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            options,
        }
    }
}

impl DaDn for NoRLowerTh1 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a}",
            units = NoRLowerTh1::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = kmax - kmin;

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a)
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        NoRLowerTh1::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl NoRLowerTh2 {
    pub const NAME: &'static str = "no_r_low_th2";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::deltak_th,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> NoRLowerTh2 {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        NoRLowerTh2 {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            options,
        }
    }
}

impl DaDn for NoRLowerTh2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} ΔK^{a} (1 - ({delta_kth} / ΔK)^{c})",
            units = NoRLowerTh2::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            c = self.c,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = kmax - kmin;

        if delta_k > self.options.deltak_th {
            self.b * delta_k.powf(self.a) * (1.0 - (self.options.deltak_th / delta_k).powf(self.c)) 
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        NoRLowerTh2::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl NoRLowerUpperTh1 {
    pub const NAME: &'static str = "no_r_low_up_th1";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::deltak_th,
        ParameterLabel::k_ut,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> NoRLowerUpperTh1 {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        NoRLowerUpperTh1 {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            k_ut: params[&ParameterLabel::k_ut],
            options,
        }
    }
}

impl DaDn for NoRLowerUpperTh1 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} (ΔK - {delta_kth})^{a} * (1 - (Kmax / {k_ut}))^{c}",
            units = NoRLowerUpperTh1::UNITS,
            b = self.b,
            delta_kth = self.options.deltak_th,
            a = self.a,
            k_ut = self.k_ut,
            c = self.c,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        let delta_k = kmax - kmin;

        if delta_k > self.options.deltak_th {
            self.b * (delta_k - self.options.deltak_th).powf(self.a) * (1.0 - (kmax / self.k_ut)).powf(self.c)
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        self.k_ut = params[&ParameterLabel::k_ut];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        NoRLowerUpperTh1::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}

impl NoRLowerUpperTh2 {
    pub const NAME: &'static str = "no_r_low_up_th2";
    pub const UNITS: &'static str = "m";

    const PARAMETER_LABELS: &'static [ParameterLabel] = &[
        ParameterLabel::a,
        ParameterLabel::b,
        ParameterLabel::c,
        ParameterLabel::d,
        ParameterLabel::deltak_th,
        ParameterLabel::k_ut,
    ];

    pub fn new(params: &BTreeMap<ParameterLabel, f64>, mut options: Options) -> NoRLowerUpperTh2 {
        // Override the global deltak_th
        options.deltak_th = params[&ParameterLabel::deltak_th];

        NoRLowerUpperTh2 {
            a: params[&ParameterLabel::a],
            b: params[&ParameterLabel::b],
            c: params[&ParameterLabel::c],
            d: params[&ParameterLabel::d],
            k_ut: params[&ParameterLabel::k_ut],
            options,
        }
    }
}

impl DaDn for NoRLowerUpperTh2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut _e = writeln!(f, "{:?}", self);
        write!(
            f,
            "#  da/dN ({units}) = {b:e} ΔK^{a} * (1 - ({deltak_th} / ΔK)^{c}) * (1 / (1 - (Kmax / {k_ut})^{d}))",
            units = NoRLowerUpperTh2::UNITS,
            b = self.b,
            deltak_th = self.options.deltak_th,
            a = self.a,
            k_ut = self.k_ut,
            c = self.c,
            d = self.d,
        )
    }

    fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
        if kmax <= self.options.kmax_th {
            return 0.0
        }

        let delta_k = kmax - kmin;

        if delta_k > self.options.deltak_th {
            self.b * delta_k.powf(self.a)
            * (1.0 - (self.options.deltak_th / delta_k).powf(self.c))
            * (1.0 / (1.0 - (kmax / self.k_ut).powf(self.d)))
        } else {
            0.0
        }
    }

    fn update_parameters(&mut self, params: &BTreeMap<ParameterLabel, f64>, deltak_th: Option<f64>) {
        self.a = params[&ParameterLabel::a];
        self.b = params[&ParameterLabel::b];
        self.c = params[&ParameterLabel::c];
        self.d = params[&ParameterLabel::d];
        self.options.deltak_th = params[&ParameterLabel::deltak_th];
        self.k_ut = params[&ParameterLabel::k_ut];
        update_kmax_th(&mut self.options, deltak_th);
    }

    fn get_name(&self) -> &str {
        NoRLowerUpperTh2::NAME
    }

    fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
        Box::new(self.clone())
    }
}



/// Create a function closure for a dadn model.
///
/// The dadn model name is specified by two parts separated by a colon
/// e.g.  model:material e.g. forman:2024t3-sheet.
/// The data component must be one of the pre-existing models. But if
/// parameters are supplied these will be installed into the dadn
/// model and used.
pub fn make_model(
    model_name: &str,
    params: &BTreeMap<ParameterLabel, f64>,
    options: Options,
) -> Result<Box<dyn DaDn + Send + Sync>, String> {
    let components: Vec<&str> = model_name.split(':').collect();
    let dadn = components[0];
    // This is the material properties string, but here we only use it for the tabular dadn
    let file_name = *components.get(1).unwrap_or(&"");

    // Verify parameters and confirm a valid equation has been given
    let parameters_verified = match dadn {
        White::NAME => verify_parameters(params, White::PARAMETER_LABELS),
        Forman::NAME => verify_parameters(params, Forman::PARAMETER_LABELS),
        FormanAFGROW::NAME => verify_parameters(params, FormanAFGROW::PARAMETER_LABELS),
        Paris::NAME => verify_parameters(params, Paris::PARAMETER_LABELS),
        Hartman::NAME => verify_parameters(params, Hartman::PARAMETER_LABELS),
        Walker::NAME => verify_parameters(params, Walker::PARAMETER_LABELS),
        Burchill::NAME => verify_parameters(params, Burchill::PARAMETER_LABELS),
        Kujawski::NAME => verify_parameters(params, Kujawski::PARAMETER_LABELS),
        Nasgro::NAME => verify_parameters(params, Nasgro::PARAMETER_LABELS),
        LinLowerTh5::NAME => verify_parameters(params, LinLowerTh5::PARAMETER_LABELS),
        LinLowerTh5A::NAME => verify_parameters(params, LinLowerTh5A::PARAMETER_LABELS),
        LinLowerTh5B::NAME => verify_parameters(params, LinLowerTh5B::PARAMETER_LABELS),
        LinLowerTh6::NAME => verify_parameters(params, LinLowerTh6::PARAMETER_LABELS),
        LinLowerTh6A::NAME => verify_parameters(params, LinLowerTh6A::PARAMETER_LABELS),
        NoRLowerTh1::NAME => verify_parameters(params, NoRLowerTh1::PARAMETER_LABELS),
        NoRLowerTh2::NAME => verify_parameters(params, NoRLowerTh2::PARAMETER_LABELS),
        NoRLowerUpperTh1::NAME => verify_parameters(params, NoRLowerUpperTh1::PARAMETER_LABELS),
        NoRLowerUpperTh2::NAME => verify_parameters(params, NoRLowerUpperTh2::PARAMETER_LABELS),
        Tabular::NAME => true,
        _ => {
            return Err(format!("Unknown dadn equation: {:?}", model_name));
        }
    };

    if !parameters_verified {
        return Err(format!(
            "Incorrect parameters given for {}: {:?}",
            dadn, params
        ));
    }

    // Construct the model
    match dadn {
        White::NAME => Ok(Box::new(White::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Forman::NAME => Ok(Box::new(Forman::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        FormanAFGROW::NAME => Ok(Box::new(FormanAFGROW::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Paris::NAME => Ok(Box::new(Paris::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Hartman::NAME => Ok(Box::new(Hartman::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Walker::NAME => Ok(Box::new(Walker::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Burchill::NAME => Ok(Box::new(Burchill::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Kujawski::NAME => Ok(Box::new(Kujawski::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Nasgro::NAME => Ok(Box::new(Nasgro::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        LinLowerTh5::NAME => {Ok(Box::new(LinLowerTh5::new(params, options)) as Box<dyn DaDn + Send + Sync>)},
        LinLowerTh5A::NAME => {Ok(Box::new(LinLowerTh5A::new(params, options)) as Box<dyn DaDn + Send + Sync>)},
        LinLowerTh5B::NAME => {Ok(Box::new(LinLowerTh5B::new(params, options)) as Box<dyn DaDn + Send + Sync>)},
        LinLowerTh6::NAME => {Ok(Box::new(LinLowerTh6::new(params, options)) as Box<dyn DaDn + Send + Sync>)},
        LinLowerTh6A::NAME => {Ok(Box::new(LinLowerTh6A::new(params, options)) as Box<dyn DaDn + Send + Sync>)},
        NoRLowerTh1::NAME => Ok(Box::new(NoRLowerTh1::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        NoRLowerTh2::NAME => Ok(Box::new(NoRLowerTh2::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        NoRLowerUpperTh1::NAME => Ok(Box::new(NoRLowerUpperTh1::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        NoRLowerUpperTh2::NAME => Ok(Box::new(NoRLowerUpperTh2::new(params, options)) as Box<dyn DaDn + Send + Sync>),
        Tabular::NAME => Ok(Box::new(Tabular::new(file_name, options)) as Box<dyn DaDn + Send + Sync>),
        _ => Err(format!("Cannot create dadn equation: {:?}", model_name)),
    }
}

/// Build a dadn model from a file
///
/// If params is defined, then the table values will be updated to use them
pub fn make_table_from_file(file_name: &str, params: &[f64]) -> table::Table {
    let table = table::Table::read_file(file_name, true);
    // put the params back into the table
    let mut log_table = table::Table::new(
        table.columns,
        table.row.iter().map(|x| x.log10()).collect::<Vec<f64>>(),
        table.values,
        true,
    );

    if !params.is_empty() {
        if params.len() == log_table.values.iter().fold(0, |sum, x| sum + x.len()) {
            log_table.update(params);
        } else {
            error!(
                "Error: the number of parameters {} != the no. of table values {}",
                params.len(),
                log_table.values.len()
            );
            process::exit(1);
        }
    }

    log_table
}

/// Verify that the parameters given on the command line are complete and correct
pub fn verify_parameters(
    given: &BTreeMap<ParameterLabel, f64>,
    expected: &'static [ParameterLabel],
) -> bool {
    // Check for length first
    if given.len() != expected.len() {
        return false;
    }

    // If at any point, an expected key is not found, fail
    for key in expected {
        if !given.contains_key(key) {
            return false;
        }
    }

    true
}

pub fn relabel_parameters(data: &[f64], model_name: &str) -> Result<BTreeMap<ParameterLabel, f64>, String> {
    let components: Vec<&str> = model_name.split(':').collect();
    let model_name = components[0];

    let mut parameter_labels = match model_name {
        White::NAME => White::PARAMETER_LABELS.to_owned(),
        Forman::NAME => Forman::PARAMETER_LABELS.to_owned(),
        FormanAFGROW::NAME => FormanAFGROW::PARAMETER_LABELS.to_owned(),
        Paris::NAME => Paris::PARAMETER_LABELS.to_owned(),
        Hartman::NAME => Hartman::PARAMETER_LABELS.to_owned(),
        Walker::NAME => Walker::PARAMETER_LABELS.to_owned(),
        Burchill::NAME => Burchill::PARAMETER_LABELS.to_owned(),
        Kujawski::NAME => Kujawski::PARAMETER_LABELS.to_owned(),
        Nasgro::NAME => Nasgro::PARAMETER_LABELS.to_owned(),
        LinLowerTh5::NAME => LinLowerTh5::PARAMETER_LABELS.to_owned(),
        LinLowerTh5A::NAME => LinLowerTh5A::PARAMETER_LABELS.to_owned(),
        LinLowerTh5B::NAME => LinLowerTh5B::PARAMETER_LABELS.to_owned(),
        LinLowerTh6::NAME => LinLowerTh6::PARAMETER_LABELS.to_owned(),
        LinLowerTh6A::NAME => LinLowerTh6A::PARAMETER_LABELS.to_owned(),
        NoRLowerTh1::NAME => NoRLowerTh1::PARAMETER_LABELS.to_owned(),
        NoRLowerTh2::NAME => NoRLowerTh2::PARAMETER_LABELS.to_owned(), 
        NoRLowerUpperTh1::NAME => NoRLowerUpperTh1::PARAMETER_LABELS.to_owned(),
        NoRLowerUpperTh2::NAME => NoRLowerUpperTh2::PARAMETER_LABELS.to_owned(),
        Tabular::NAME => return Ok(BTreeMap::new()),    // This is not ideal
        _ => {
            return Err(format!("Unknown dadn equation: {:?}", model_name));
        }
    };

    if data.len() != parameter_labels.len() {
        return Err(format!("Number of values does not match number of labels: {:?} {:?}", data, parameter_labels));
    }

    parameter_labels.sort();

    let mut remapped_params = BTreeMap::new();
    for (i, name) in parameter_labels.iter().enumerate() {
        remapped_params.insert(name.clone(), data[i]);
    }

    Ok(remapped_params)
}

pub mod testing_objects {
    use super::*;
    
    /// A test dadn which returns a constant value
    pub struct Constant {
        value: f64
    }

    /// A test dadn which will return kmax - kmin 
    pub struct Subtractive {}

    impl Constant {
        pub const NAME: &'static str = "constant";

        pub fn new(value: f64) -> Self {
            Self {
                value,
            }
        }
    }

    impl Subtractive {
        pub const NAME: &'static str = "subtractive";
    }

    impl DaDn for Constant {
        fn fmt(&self, _f: &mut fmt::Formatter) -> fmt::Result {
            todo!()
        }

        fn dadn(&self, _kmin: f64, _kmax: f64, _state: CrackState) -> f64 {
            self.value
        }

        fn update_parameters(&mut self, _params: &BTreeMap<ParameterLabel, f64>, _deltak_th: Option<f64>) {
            todo!()
        }

        fn get_name(&self) -> &str {
            Constant::NAME
        }

        fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
            todo!()
        }
    }

    impl DaDn for Subtractive {
        fn fmt(&self, _f: &mut fmt::Formatter) -> fmt::Result {
            todo!()
        }

        fn dadn(&self, kmin: f64, kmax: f64, _state: CrackState) -> f64 {
            kmax - kmin
        }

        fn update_parameters(&mut self, _params: &BTreeMap<ParameterLabel, f64>, _deltak_th: Option<f64>) {
            todo!()
        }

        fn get_name(&self) -> &str {
            Constant::NAME
        }

        fn inner_clone(&self) -> Box<dyn DaDn + Send + Sync> {
            todo!()
        }
    }

    /// Get a boxed dadn which will return a constant value from dadn()
    /// 
    /// Currently only viable for testing in non-optimisation scenarios
    pub fn make_constant(value: f64) -> Box<dyn DaDn + Send + Sync> {
        Box::new(Constant::new(value)) as Box<dyn DaDn + Send + Sync>
    }

    /// Get a boxed dadn which will alwys return kmax - kmin from dadn()
    /// 
    /// Currently only viable for testing in non-optimisation scenarios
    pub fn make_subtractive() -> Box<dyn DaDn + Send + Sync> {
        Box::new(Subtractive{}) as Box<dyn DaDn + Send + Sync>
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material;
    use table;

    #[test]
    fn white_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("white:barter14-aa7050t7451").unwrap().params;
        let options = Options::default();
        let equation = White::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.1494279345488723e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn forman_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("forman:default").unwrap().params;
        let options = Options::default();
        let equation = Forman::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.6949152542372881e-12).abs() <= f64::EPSILON)
    }

    #[test]
    fn forman_afgrow_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("forman-ag:default").unwrap().params;
        let options = Options::default();
        let equation = FormanAFGROW::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.6949152542372881e-12).abs() <= f64::EPSILON)
    }

    #[test]
    fn paris_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("paris:default").unwrap().params;
        let options = Options::default();
        let equation = Paris::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1e-9).abs() <= f64::EPSILON)
    }

    #[test]
    fn hartman_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("hartman:default").unwrap().params;
        let options = Options::default();
        let equation = Hartman::new(&params, options);

        assert!(equation.dadn(0.0, 1.0, CrackState { a: 0.0 }).abs() <= f64::EPSILON)
    }

    #[test]
    fn walker_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("walker:default").unwrap().params;
        let options = Options::default();
        let equation = Walker::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn burchill_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("burchill:default").unwrap().params;
        let options = Options::default();
        let equation = Burchill::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn kujawski_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("kujawski:default").unwrap().params;
        let options = Options::default();
        let equation = Kujawski::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn nasgro_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("nasgro:default").unwrap().params;
        let options = Options::default();
        let equation = Nasgro::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 2.4406636167321105e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn lin_lower_th_5b_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("lin_lower_th_5b:default").unwrap().params;
        let options = Options::default();
        let equation = LinLowerTh5B::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.1142377317870345e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn no_r_low_th1_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("no_r_low_th1:aa7075t7351").unwrap().params;
        let options = Options::default();
        let equation = NoRLowerTh1::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 7.696890397044252e-11).abs() <= f64::EPSILON)
    }

    #[test]
    fn no_r_low_th2_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("no_r_low_th2:aa7075t7351").unwrap().params;
        let options = Options::default();
        let equation = NoRLowerTh2::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.1832687996044777e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn no_r_low_up_th1_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("no_r_low_up_th1:aa7075t7351").unwrap().params;
        let options = Options::default();
        let equation = NoRLowerUpperTh1::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.1606733952853102e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn no_r_low_up_th2_dadn_returns_correct_value_for_material() {
        let params = &material::get_dadn("no_r_low_up_th2:aa7075t7351").unwrap().params;
        let options = Options::default();
        let equation = NoRLowerUpperTh2::new(&params, options);

        assert!((equation.dadn(0.0, 1.0, CrackState { a: 0.0 }) - 1.228206125308067e-10).abs() <= f64::EPSILON)
    }

    #[test]
    fn forman_dadn_returns_zero_when_kmax_less_than_or_equal_to_zero() {
        let params = &material::get_dadn("forman:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = Forman::new(&params, options);

        assert!(equation.dadn(0.0, -1.0, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn forman_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_false() {
        let params = &material::get_dadn("forman:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = Forman::new(&params, options);
        println!("{:e}", equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }));

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 1.1428571428571429e-10).abs() <= f64::EPSILON);
    }

    #[test]
    fn forman_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_true() {
        let params = &material::get_dadn("forman:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: true,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = Forman::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 2.5714285714285719e-10).abs() <= f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6a_new_correctly_allocates_variables() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.1),
            (ParameterLabel::b, 2.2),
            (ParameterLabel::c, 3.3),
            (ParameterLabel::d, 4.4),
            (ParameterLabel::e, 5.5),
            (ParameterLabel::deltak_th, 6.6),
            (ParameterLabel::k_ut, 7.7),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&params, options);
        assert!(equation.a - 1.1 < f64::EPSILON);
        assert!(equation.b - 2.2 < f64::EPSILON);
        assert!(equation.c - 3.3 < f64::EPSILON);
        assert!(equation.d - 4.4 < f64::EPSILON);
        assert!(equation.e - 5.5 < f64::EPSILON);
        assert!(equation.options.deltak_th - 6.6 < f64::EPSILON);
        assert!(equation.k_ut - 7.7 < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6a_new_correctly_allocates_default_variables() {
        let material = material::get_dadn("lin_lower_th_6a:default").unwrap();
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&material.params, options);
        assert!((equation.a - 1.9544).abs() < f64::EPSILON);
        assert!((equation.b - 7.3345e-10).abs() < f64::EPSILON);
        assert!((equation.c - -0.2724).abs() < f64::EPSILON);
        assert!((equation.d - 0.0589).abs() < f64::EPSILON);
        assert!((equation.e - -6165.6640).abs() < f64::EPSILON);
        assert!((equation.options.deltak_th - 0.6709).abs() < f64::EPSILON);
        assert!((equation.k_ut - 61964.1148).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6a_dadn_returns_zero_when_kmax_less_than_or_equal_to_zero() {
        let params = &material::get_dadn("lin_lower_th_6a:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };
        
        let equation = LinLowerTh6A::new(&params, options);

        assert!(equation.dadn(0.0, -1.0, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6a_dadn_returns_zero_when_delta_k_less_than_or_equal_to_delta_kth() {
        let params = &material::get_dadn("lin_lower_th_6a:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&params, options);

        // By halving the supplied deltak_th, we can ensure a delta_k <= deltak_th
        let kmax = equation.options.deltak_th / 2.0;

        assert!(equation.dadn(0.0, kmax, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6a_dadn_returns_correct_value_when_r_less_than_or_equal_to_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.594),
            (ParameterLabel::b, 1.078e-9),
            (ParameterLabel::c, 0.2418),
            (ParameterLabel::d, -0.06825),
            (ParameterLabel::e, -9129.0),
            (ParameterLabel::deltak_th, 0.6903),
            (ParameterLabel::k_ut, 78496.0),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&params, options);

        assert!((equation.dadn(6.0, 10.0, CrackState { a: 0.0 }) - 1.190185e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_6a_dadn_returns_correct_value_when_r_greater_than_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.594),
            (ParameterLabel::b, 1.078e-9),
            (ParameterLabel::c, 0.2418),
            (ParameterLabel::d, -0.06825),
            (ParameterLabel::e, -9129.0),
            (ParameterLabel::deltak_th, 0.6903),
            (ParameterLabel::k_ut, 78496.0),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&params, options);

        assert!((equation.dadn(18.0, 20.0, CrackState { a: 0.0 }) - 1.765251e-9).abs() < 1e-14);
    }

    #[test]
    fn lin_lower_th_6a_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_false() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.594),
            (ParameterLabel::b, 1.078e-9),
            (ParameterLabel::c, 0.2418),
            (ParameterLabel::d, -0.06825),
            (ParameterLabel::e, -9129.0),
            (ParameterLabel::deltak_th, 0.6903),
            (ParameterLabel::k_ut, 78496.0),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 1.142102e-8).abs() < 1e-14);
    }

    #[test]
    fn lin_lower_th_6a_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_true() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.594),
            (ParameterLabel::b, 1.078e-9),
            (ParameterLabel::c, 0.2418),
            (ParameterLabel::d, -0.06825),
            (ParameterLabel::e, -9129.0),
            (ParameterLabel::deltak_th, 0.6903),
            (ParameterLabel::k_ut, 78496.0),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: true,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6A::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 2.896738e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_6_new_correctly_allocates_variables() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.1),
            (ParameterLabel::b, 2.2),
            (ParameterLabel::c, 3.3),
            (ParameterLabel::d, 4.4),
            (ParameterLabel::deltak_th, 5.5),
        ]);
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);
        assert!(equation.a - 1.1 < f64::EPSILON);
        assert!(equation.b - 2.2 < f64::EPSILON);
        assert!(equation.c - 3.3 < f64::EPSILON);
        assert!(equation.d - 4.4 < f64::EPSILON);
        assert!(equation.options.deltak_th - 5.5 < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6_new_correctly_allocates_default_variables() {
        let material = material::get_dadn("lin_lower_th_6:default").unwrap();
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&material.params, options);
        assert!((equation.a - 2.4299).abs() < f64::EPSILON);
        assert!((equation.b - 6.8283e-10).abs() < f64::EPSILON);
        assert!((equation.c - -0.8636).abs() < f64::EPSILON);
        assert!((equation.d - 0.1065).abs() < f64::EPSILON);
        assert!((equation.options.deltak_th - 0.6402).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6_dadn_returns_zero_when_kmax_less_than_or_equal_to_zero() {
        let params = &material::get_dadn("lin_lower_th_6:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);

        assert!(equation.dadn(0.0, -1.0, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6_dadn_returns_zero_when_delta_k_less_than_or_equal_to_delta_kth() {
        let params = &material::get_dadn("lin_lower_th_6:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);

        // By halving the supplied deltak_th, we can ensure a delta_k <= deltak_th
        let kmax = equation.options.deltak_th / 2.0;

        assert!(equation.dadn(0.0, kmax, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_6_dadn_returns_correct_value_when_r_less_than_or_equal_to_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.124),
            (ParameterLabel::b, 8.788e-10),
            (ParameterLabel::c, -0.3688),
            (ParameterLabel::d, -0.1321),
            (ParameterLabel::deltak_th, 0.6466),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);

        assert!((equation.dadn(6.0, 10.0, CrackState { a: 0.0 }) - 2.612479e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_6_dadn_returns_correct_value_when_r_greater_than_rmax() {
        //let params = &material::get_dadn("lin_lower_th_6:default").unwrap().params;
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.124),
            (ParameterLabel::b, 8.788e-10),
            (ParameterLabel::c, -0.3688),
            (ParameterLabel::d, -0.1321),
            (ParameterLabel::deltak_th, 0.6466),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);

        assert!((equation.dadn(18.0, 20.0, CrackState { a: 0.0 }) - 4.629001e-9).abs() < 1e-14);
    }

    #[test]
    fn lin_lower_th_6_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_false() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.124),
            (ParameterLabel::b, 8.788e-10),
            (ParameterLabel::c, -0.3688),
            (ParameterLabel::d, -0.1321),
            (ParameterLabel::deltak_th, 0.6466),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 7.980533e-9).abs() < 1e-14);
    }

    #[test]
    fn lin_lower_th_6_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_true() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.124),
            (ParameterLabel::b, 8.788e-10),
            (ParameterLabel::c, -0.3688),
            (ParameterLabel::d, -0.1321),
            (ParameterLabel::deltak_th, 0.6466),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: true,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh6::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 1.936365e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_5a_new_correctly_allocates_variables() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.1),
            (ParameterLabel::b, 2.2),
            (ParameterLabel::c, 3.3),
            (ParameterLabel::d, 4.4),
            (ParameterLabel::deltak_th, 5.5),
            (ParameterLabel::k_ut, 6.6),
        ]);
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&params, options);
        assert!(equation.a - 1.1 < f64::EPSILON);
        assert!(equation.b - 2.2 < f64::EPSILON);
        assert!(equation.c - 3.3 < f64::EPSILON);
        assert!(equation.d - 4.4 < f64::EPSILON);
        assert!(equation.options.deltak_th - 5.5 < f64::EPSILON);
        assert!(equation.k_ut - 6.6 < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5a_new_correctly_allocates_default_variables() {
        let material = material::get_dadn("lin_lower_th_5a:default").unwrap();
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&material.params, options);
        assert!((equation.a - 1.8476).abs() < f64::EPSILON);
        assert!((equation.b - 7.7723e-10).abs() < f64::EPSILON);
        assert!((equation.c - -0.0403).abs() < f64::EPSILON);
        assert!((equation.d - -998.5654).abs() < f64::EPSILON);
        assert!((equation.options.deltak_th - 0.6717).abs() < f64::EPSILON);
        assert!((equation.k_ut - 8645.19).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5a_dadn_returns_zero_when_kmax_less_than_or_equal_to_zero() {
        let params = &material::get_dadn("lin_lower_th_5a:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&params, options);

        assert!(equation.dadn(0.0, -1.0, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5a_dadn_returns_zero_when_delta_k_less_than_or_equal_to_delta_kth() {
        let params = &material::get_dadn("lin_lower_th_5a:default").unwrap().params;
        let options = Options {
            rmin: -1.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };
    
        let equation = LinLowerTh5A::new(&params, options);

        // By halving the supplied deltak_th, we can ensure a delta_k <= deltak_th
        let kmax = equation.options.deltak_th / 2.0;

        assert!(equation.dadn(0.0, kmax, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5a_dadn_returns_correct_value_when_r_less_than_or_equal_to_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.6220),
            (ParameterLabel::b, 1.0348e-9),
            (ParameterLabel::c, 0.1318),
            (ParameterLabel::d, -832.65),
            (ParameterLabel::deltak_th, 0.6914),
            (ParameterLabel::k_ut, 6761.66),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&params, options);

        assert!((equation.dadn(6.0, 10.0, CrackState { a: 0.0 }) - 2.190127e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_5a_dadn_returns_correct_value_when_r_greater_than_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.6220),
            (ParameterLabel::b, 1.0348e-9),
            (ParameterLabel::c, 0.1318),
            (ParameterLabel::d, -832.65),
            (ParameterLabel::deltak_th, 0.6914),
            (ParameterLabel::k_ut, 6761.66),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&params, options);

        assert!((equation.dadn(18.0, 20.0, CrackState { a: 0.0 }) - 1.525398e-8).abs() < 1e-14);
    }

    #[test]
    fn lin_lower_th_5a_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_false() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.6220),
            (ParameterLabel::b, 1.0348e-9),
            (ParameterLabel::c, 0.1318),
            (ParameterLabel::d, -832.65),
            (ParameterLabel::deltak_th, 0.6914),
            (ParameterLabel::k_ut, 6761.66),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 1.244264e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_5a_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_true() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.6220),
            (ParameterLabel::b, 1.0348e-9),
            (ParameterLabel::c, 0.1318),
            (ParameterLabel::d, -832.65),
            (ParameterLabel::deltak_th, 0.6914),
            (ParameterLabel::k_ut, 6761.66),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: true,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5A::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 2.678966e-8).abs() < 1e-13);
    }

    #[test]
    fn lin_lower_th_5_new_correctly_allocates_variables() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 1.1),
            (ParameterLabel::b, 2.2),
            (ParameterLabel::c, 3.3),
            (ParameterLabel::deltak_th, 4.4),
        ]);
        let options = Options {
            rmin: 0.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);
        assert!(equation.a - 1.1 < f64::EPSILON);
        assert!(equation.b - 2.2 < f64::EPSILON);
        assert!(equation.c - 3.3 < f64::EPSILON);
        assert!(equation.options.deltak_th - 4.4 < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5_new_correctly_allocates_default_variables() {
        let material = material::get_dadn("lin_lower_th_5:default").unwrap();
        let options = Options {
            rmin: 0.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&material.params, options);
        assert!((equation.a - 2.3805).abs() < f64::EPSILON);
        assert!((equation.b - 7.3607e-10).abs() < f64::EPSILON);
        assert!((equation.c - -0.5799).abs() < f64::EPSILON);
        assert!((equation.options.deltak_th - 0.6238).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5_dadn_returns_zero_when_kmax_less_than_or_equal_to_zero() {
        let params = &material::get_dadn("lin_lower_th_5:default").unwrap().params;
        let options = Options {
            rmin: 0.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);

        assert!(equation.dadn(0.0, -1.0, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5_dadn_returns_zero_when_delta_k_less_than_or_equal_to_delta_kth() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.2317),
            (ParameterLabel::b, 8.0567e-10),
            (ParameterLabel::c, -0.6636),
            (ParameterLabel::deltak_th, 0.6506),
        ]);
        
        let options = Options {
            rmin: 0.0,
            rmax: 1.0,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);

        // By halving the supplied deltak_th, we can ensure a delta_k <= deltak_th
        let kmax = equation.options.deltak_th / 2.0;

        assert!(equation.dadn(0.0, kmax, CrackState { a: 0.0 }).abs() < f64::EPSILON);
    }

    #[test]
    fn lin_lower_th_5_dadn_returns_correct_value_when_r_less_than_or_equal_to_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.2317),
            (ParameterLabel::b, 8.0567e-10),
            (ParameterLabel::c, -0.6636),
            (ParameterLabel::deltak_th, 0.6506),
        ]);

        let options = Options {
            rmin: 0.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);

        assert!((equation.dadn(6.0, 10.0, CrackState { a: 0.0 }) - 2.196843e-8).abs() < 10e-11);
    }

    #[test]
    fn lin_lower_th_5_dadn_returns_correct_value_when_r_greater_than_rmax() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.2317),
            (ParameterLabel::b, 8.0567e-10),
            (ParameterLabel::c, -0.6636),
            (ParameterLabel::deltak_th, 0.6506),
        ]);
        
        let options = Options {
            rmin: 0.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);

        assert!((equation.dadn(18.0, 20.0, CrackState { a: 0.0 }) - 4.575396e-9).abs() < 10e-12);
    }

    #[test]
    fn lin_lower_th_5_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_false() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.2317),
            (ParameterLabel::b, 8.0567e-10),
            (ParameterLabel::c, -0.6636),
            (ParameterLabel::deltak_th, 0.6506),
        ]);

        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 9.138437e-9).abs() < 10e-12);
    }

    #[test]
    fn lin_lower_th_5_dadn_returns_correct_value_when_kmin_less_than_zero_and_kneg_is_true() {
        let params = BTreeMap::from([
            (ParameterLabel::a, 2.2317),
            (ParameterLabel::b, 8.0567e-10),
            (ParameterLabel::c, -0.6636),
            (ParameterLabel::deltak_th, 0.6506),
        ]);
        
        let options = Options {
            rmin: -1.0,
            rmax: 0.8,
            kneg: true,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let equation = LinLowerTh5::new(&params, options);

        assert!((equation.dadn(-2.0, 4.0, CrackState { a: 0.0 }) - 2.598134e-8).abs() < 10e-11);
    }

    #[test]
    fn verify_parameters_returns_true_when_expected_matches_given() {
        let expected: &'static [ParameterLabel] =
            &[ParameterLabel::a, ParameterLabel::b, ParameterLabel::c];

        let given = BTreeMap::from([
            (ParameterLabel::a, 1.1),
            (ParameterLabel::b, 2.2),
            (ParameterLabel::c, 3.3),
        ]);

        assert!(verify_parameters(&given, expected));
    }

    #[test]
    fn verify_parameters_returns_false_when_expected_does_not_match_given() {
        let expected: &'static [ParameterLabel] = &[ParameterLabel::a];

        let given = BTreeMap::from([(ParameterLabel::b, 1.1)]);

        assert!(!verify_parameters(&given, expected));
    }

    #[test]
    fn verify_parameters_returns_false_when_different_number_of_parameters() {
        let expected: &'static [ParameterLabel] = &[ParameterLabel::a, ParameterLabel::b];

        let given = BTreeMap::from([
            (ParameterLabel::a, 1.1),
            (ParameterLabel::b, 2.2),
            (ParameterLabel::c, 3.3),
        ]);

        assert!(!verify_parameters(&given, expected));
    }

    #[test]
    fn check_kujawski() {
        let params = match material::get_dadn("kujawski:default") {
            Some(m) => &m.params,
            None => panic!(),
        };
        let options = Options::default();

        let kujawski_eqn = Kujawski::new(&params, options.to_owned());

        // 6 ksi sqrt(in) gives dadn = 1.25e-6 in
        let da = kujawski_eqn.dadn(0.0, 6.6, CrackState { a: 0.0 });
        println!("kujawski da/dn {}", da);
        assert!((da - 2.87496e-8).abs() < 1.0e-10);

        let params = BTreeMap::from([
            (ParameterLabel::c, 1e-10),
            (ParameterLabel::m, 3.0),
            (ParameterLabel::alpha, 0.25),
        ]);

        let kujawski_eqn = Kujawski::new(&params, options);
        let da = kujawski_eqn.dadn(0.0, 6.6, CrackState { a: 0.0 });
        println!("kujawski da/dn {}", da);
        assert!((da - 2.87496e-8).abs() < 1.0e-10);

        let da = kujawski_eqn.dadn(3.0, 6.6, CrackState { a: 0.0 });
        println!("kujawski da/dn {}", da);
        assert!((da - 0.735_086_7e-8).abs() < 1.0e-10);
    }

    #[test]
    fn check_tabular() {
        let table = table::Table::new(
            vec![0.0, 1.0],
            // already logged.
            vec![-8.0, -7.0, -6.0, -5.0],
            // these are the columns
            vec![vec![10.0, 15.0, 20.0, 25.0], vec![5.0, 7.5, 10.0, 12.5]],
            true,
        );

        let options = Options {
            rmin: 0.0,
            rmax: 0.99,
            kneg: false,
            deltak_th: 0.0,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let tabular = Tabular {
            table,
            options,
        };

        assert!((tabular.dadn(0.0, 10.0, CrackState { a: 0.0 }) - 1e-8).abs() < std::f64::EPSILON);
        println!("dadn: min -5.0, max 10.0");
        assert!((tabular.dadn(-5.0, 10.0, CrackState { a: 0.0 }) - 1e-8).abs() < std::f64::EPSILON);

        println!("dadn: min 0.0, max 12.5");
        assert!(
            (tabular.dadn(0.0, 12.5, CrackState { a: 0.0 }) - 3.162_277_660_168_379e-08).abs()
                < std::f64::EPSILON
        );
        assert!(
            (tabular.dadn(12.5, 25.0, CrackState { a: 0.0 }) - 0.000_000_562_341_325_190_349).abs()
                < std::f64::EPSILON
        );
    }

    #[test]
    fn check_walker_rratio() {
        // check that changing the r exponent makes no difference when the r ratio is 0.0
        // this makes sure we know which exponent is which
        let params = BTreeMap::from([
            (ParameterLabel::c, 1e-8),
            (ParameterLabel::m, 0.5),
            (ParameterLabel::n, 2.0),
        ]);
        let options = Options::default();

        let w1 = Walker::new(&params, options.to_owned());
        let a1 = w1.dadn(0.0, 10.0, CrackState { a: 0.1 });

        let params = BTreeMap::from([
            (ParameterLabel::c, 1e-8),
            (ParameterLabel::m, 2.5),
            (ParameterLabel::n, 2.0),
        ]);
        let w2 = Walker::new(&params, options);
        let a2 = w2.dadn(0.0, 10.0, CrackState { a: 0.1 });

        assert!((a1 - a2).abs() < 1e-8);
    }

    #[test]
    fn check_nasgro() {
        let params = match material::get_dadn("nasgro:default") {
            Some(m) => &m.params,
            None => panic!(),
        };
        let options = Options::default();

        let nasgro_eqn = Nasgro::new(&params, options);

        // 6 ksi sqrt(in) gives dadn = 1.25e-6 in
        let da = nasgro_eqn.dadn(0.0, 6.6, CrackState { a: 0.0 });
        println!("nasgro da/dn {}", da);
        assert!(((da - 3.2668e-8).abs() / da) < 1.0e-3);
    }

    //#[test]
    fn _check_walker_closed_vs_tabular() {
        let options = Options::default();
        let params = BTreeMap::from([
            (ParameterLabel::c, 1.00e-10),
            (ParameterLabel::m, 0.5),
            (ParameterLabel::n, 3.0),
        ]);

        let walker = Walker::new(&params, options.to_owned());
        let tabular = Tabular::new("walker_default.txt", options);

        // R = 0.0,  dk = 7.937005259841,   da = 5e-8
        // R = 0.5,  dk = 5.61231024154686, da = 5e-8
        // R = 0.45, dk = 3.44229440855592, da = 1e-8

        let r = 0.45;
        let dk = 58.8624064001031;
        let kmax = dk / (1.0 - r);
        let kmin = r * kmax;
        let da = walker.dadn(kmin, kmax, CrackState { a: 0.0 });
        let tab_da = tabular.dadn(kmin, kmax, CrackState { a: 0.0 });

        println!("kmax {} kmin {} diff {}, da = {:.3e}", kmax, kmin, kmax - kmin, da);
        println!("da {:.3e} tab_da {:.3e} diff {:.3}%", da, tab_da, (tab_da - da) / da);
    }

    //#[test]
    fn _compare_lin_lower_th_5b_with_tabular() {
        // a=1.99576,b=4.04968e-10,d=-3793.51,c1=-0.41677,c2=-1.60455,deltak_th=8.4017e-05,k_ut=87910.2 --rmax 0.9 --rmin -1

        let options = Options {
            rmin: -1.0,
            rmax: 0.9,
            kneg: true,
            deltak_th: 8.4017e-05,
            kmax_th: 0.0,
            kmax_th_follows_deltak_th: true,
        };

        let params = BTreeMap::from([
            (ParameterLabel::a, 1.99576),
            (ParameterLabel::b, 4.04968e-10),
            (ParameterLabel::d, -3793.51),
            (ParameterLabel::c1, -0.41677),
            (ParameterLabel::c2, -1.60455),
            (ParameterLabel::deltak_th, options.deltak_th),
            (ParameterLabel::k_ut, 87910.2),
        ]);
    
        let linlower5b = LinLowerTh5B::new(&params, options.to_owned());
        let tabular = Tabular::new("dadN_DK_gen_knegT_widebound_PS2_new.txt", options);
        
        // r = -1.0, dk = 100.4119634, da = 1e-4
        // r = -1.0, dk = 2.596331569, da = 1e-9
        // r = -0.5, dk = 2.083204455, da = 1e-9
        // r = 0.0, dk = 1.522113689, da = 1e-9

        let r = -1.0;
        let dk = 100.4119634;
        let kmax = dk / (1.0 - r);
        let kmin = r * kmax;

        let da = linlower5b.dadn(kmin, kmax, CrackState { a: 0.0 });
        let tab_da = tabular.dadn(kmin, kmax, CrackState { a: 0.0 });

        println!("kmax {} kmin {} diff {}, da = {:.3e}", kmax, kmin, kmax - kmin, da);
        println!("da {:.3e} tab_da {:.3e} diff {:.3}%", da, tab_da, (tab_da - da) / da);

        println!("1/10 = {}", linlower5b.dadn(1.0, 10.0, CrackState { a: 0.0 }));
        println!("10/100 = {}", linlower5b.dadn(10.0, 100.0, CrackState { a: 0.0 }));
    }
}
