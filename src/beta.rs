//! A variety of beta solutions for cracks in different shaped coupons.
//!
//! Beta functions relate the far field loading to the crack tip so
//! they are a function of the geometry and of the size of the crack
//! in the component. In all the functions here we use
//! non-dimensionalised factors. All the betas can be called with the
//! same set of non-dimensionalised parameters about the crack size
//! but each particular functino will only use those vaiables that are
//! relavent to it.
//!
//! The naming convention is (crack shape, geometry, loading) - (author, publication year)
//! The letter for crack is not included in the abbreviation since it seems redundant.
//! No need to say surface for semi-elliptical or quarter cracks.

// cargo test -- --nocapture
#![allow(clippy::unreadable_literal)]
use crate::{io, table};
use std::process;
use std::f64::consts::{FRAC_PI_2, PI};
use crate::grow::{self, CrackGeometry, CrackFront};
use std::fmt;

use log::{error, debug};

const PHI_A: f64 = 0.0;
const PHI_C: f64 = FRAC_PI_2;

#[derive(Debug, Clone)]
pub enum DirectionOfInterest {
    A,
    C,
}

impl fmt::Display for DirectionOfInterest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

/// Result of a beta function
#[derive(Debug, Clone)]
pub struct BetaResult {
    pub a: Option<f64>,
    pub c: Option<f64>,
}

impl BetaResult {
    /// Extract the beta for the given direction of interest
    pub fn beta_of_interest(&self, direction_of_interest: &DirectionOfInterest) -> f64 {
        match direction_of_interest {
            DirectionOfInterest::A => self.a.unwrap(),
            DirectionOfInterest::C => self.c.unwrap(),
        }
    }
}


// This trait is constructed for betas to allow extra information to
// be supplied with the betas such as required for a tabular beta or
// the coupon beta which uses the dimensions of the coupon.
pub trait Beta {
    /// Construct an initial geometry
    fn geometry(&self, a: f64, c: f64) -> CrackGeometry;

    /// Calculate a beta factor for some number of crack fronts
    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)>;

    fn area(&self, _a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        1.0
    }

    fn direction_of_interest(&self) -> &DirectionOfInterest;

    fn set_b(&mut self, _b: f64) {}  // Default implementation does nothing

    fn set_d(&mut self, _d: f64) {}

    /// Convert a beta function into a beta table    
    fn as_table(&mut self) -> TableBeta {
        let a_on_cs = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let ratios = (0..=99).map(|x| x as f64 * 0.01).collect::<Vec<_>>();
        
        let direction = self.direction_of_interest().to_owned();
        let mut geometry = self.geometry(0.0, 0.0);
        let mut beta_result;
        let mut values = Vec::new();

        for a_on_c in &a_on_cs {
            let mut column = Vec::new();
            for ratio in &ratios {
                match direction {
                    DirectionOfInterest::A => {
                        geometry.set_length_a(*a_on_c);
                        geometry.set_length_c(1.0);
                        self.set_d(a_on_c / ratio);
                    } 
                    DirectionOfInterest::C => {
                        geometry.set_length_a(1.0);
                        geometry.set_length_c(*a_on_c);
                        self.set_b(a_on_c / ratio);
                    }
                };

                beta_result = self.beta(&geometry);

                match beta_result {
                    Ok(result) => {
                        column.push(result.beta_of_interest(&direction));
                    }
                    Err(_) => todo!(),
                };
            }
            values.push(column);
        }

        let table = table::Table::new(a_on_cs, ratios, values, false);
        TableBeta { table, direction_of_interest: self.direction_of_interest().to_owned(), b: 1.0, d: 1.0, ratio: 1.0 }
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync>;
}

impl Clone for Box<dyn Beta + Send + Sync> {
    fn clone(&self) -> Self {
        self.inner_clone()
    }
}

impl fmt::Display for TableBeta {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let ratio = match self.direction_of_interest() {
            DirectionOfInterest::A => "a/d",
            DirectionOfInterest::C => "c/b"
        };
        let _ = writeln!(f, "# {}{:9}", ratio, " \\ a/c");
        let _ = write!(f, "            ");
        for a_on_c in &self.table.columns {
            let _ = write!(f, "{:6.3} ", a_on_c);
        }
        let _ = writeln!(f);

        for i in 0..self.table.row.len() {
            let _ = write!(f, "{:9.3}   ", self.table.row[i]);
            for j in 0..self.table.columns.len() {
                let value = self.table.values[j][i];
                let _ = write!(f, "{:6.3} ", value);
            }
            let _ = writeln!(f);
        }
        write!(f, "")
    }
}

/// Construct a beta from its name.  If the beta it is not found
/// it will trigger an error! call.
pub fn get_beta_fn(beta_name: &str, component: &grow::Component, initial_geometry: CrackGeometry, doi: Option<DirectionOfInterest>) -> Box<dyn Beta + Send + Sync> {
    if beta_name.starts_with("file-") {
        return get_beta_from_file(beta_name, component, initial_geometry.ratio.unwrap());
    }

    match beta_name {
        QuarterBroek86::NAME => Box::new(QuarterBroek86 {}),
        SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::NAME =>
            Box::new(SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::new(
                component.sideways,
                component.forward,
                initial_geometry.a.unwrap().angle,
                initial_geometry.c.unwrap().angle,
                doi
            )),
        SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::NAME =>
            Box::new(SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::new(
                initial_geometry.a.unwrap().angle,
                initial_geometry.c.unwrap().angle,
                doi,
            )),
        // QuarterCircularCornerCrackFinitePlateTensionMurakami87::NAME =>
        //    Box::new(QuarterCircularCornerCrackFinitePlateTensionMurakami87 { d: component.forward }),
        QuarterEllipticalCornerCrackFinitePlateTensionNewman84::NAME =>
            Box::new(QuarterEllipticalCornerCrackFinitePlateTensionNewman84::new(
                component.sideways,
                component.forward,
                initial_geometry.a.unwrap().angle,
                initial_geometry.c.unwrap().angle,
                doi,
            )),
        EllipticalEmbeddedCrackFinitePlateTensionNewman84::NAME =>
            Box::new(EllipticalEmbeddedCrackFinitePlateTensionNewman84::new(
                component.sideways,
                component.forward,
                initial_geometry.a.unwrap().angle,
                initial_geometry.c.unwrap().angle,
                doi,
            )),
        SingleSidedEdgeCrackTensionTada73::NAME =>
            Box::new(SingleSidedEdgeCrackTensionTada73 { b: component.sideways }),
        DoubleSidedEdgeCrackTensionTada73::NAME =>
            Box::new(DoubleSidedEdgeCrackTensionTada73 { b: component.sideways }),
        DoubleSidedEdgeCrackTensionTada73_2::NAME =>
            Box::new(DoubleSidedEdgeCrackTensionTada73_2 { b: component.sideways }),
        CompactCoupon::NAME =>
            Box::new(CompactCoupon{ b: component.sideways, d: component.forward}),
        CentreCrackTensionFedderson66::NAME =>
            Box::new(CentreCrackTensionFedderson66 { b: component.sideways }),
        CentreCrackTensionKoiter65::NAME =>
            Box::new(CentreCrackTensionKoiter65 { b: component.sideways }),
        // CornerCrackConstrainedTensionMcdonald07::NAME =>
        //     Box::new(CornerCrackConstrainedTensionMcdonald07::new(component.sideways)),
        SingleSidedThroughCrackCircularHoleTensionBowie56::NAME =>
            Box::new(SingleSidedThroughCrackCircularHoleTensionBowie56 { r: component.radius }),
        DoubleSidedThroughCrackCircularHoleTensionBowie56::NAME =>
            Box::new(DoubleSidedThroughCrackCircularHoleTensionBowie56 { r: component.radius }),
        // DoubleSidedCornerCrackHoleTensionNewman81::NAME =>
        //    Box::new(DoubleSidedCornerCrackHoleTensionNewman81 { d: component.forward, b: component.sideways, r: component.radius }),
        SemiEllipticalSurfaceCrackRoundBarBendingMurakami87::NAME =>
            Box::new(SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 { d: component.forward, ratio: initial_geometry.ratio.unwrap() }),
        SemiEllipticalSurfaceCrackRoundBarTensionMurakami87::NAME =>
            Box::new(SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 { d: component.forward, ratio: initial_geometry.ratio.unwrap() }),
        SemiEllipticalSurfaceCrackRoundBarBendingMurakami86::NAME => 
            Box::new(SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 { d: component.forward, ratio: initial_geometry.ratio.unwrap() }),
        SemiEllipticalSurfaceCrackRoundBarTensionMurakami86::NAME =>
            Box::new(SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 { d: component.forward, ratio: initial_geometry.ratio.unwrap() }),
        EdgeCrackStripBendingMurakami87::NAME =>
            Box::new(EdgeCrackStripBendingMurakami87 { b: component.sideways }),
        EdgeCrackStripTensionMurakami87::NAME => 
            Box::new(EdgeCrackStripTensionMurakami87 { b: component.sideways }),
        SemiCircularSurfaceCrackRoundBarTensionForman86::NAME =>
            Box::new(SemiCircularSurfaceCrackRoundBarTensionForman86 { d: component.forward }),
        SemiCircularSurfaceCrackRoundBarBendingForman86::NAME =>
            Box::new(SemiCircularSurfaceCrackRoundBarBendingForman86 { d: component.forward }),
        _ => {
            error!("Error: Unknown beta '{}'", beta_name);
            process::exit(1);
        }
    }
}

fn get_beta_from_file(beta_name: &str, component: &grow::Component, initial_ratio: f64) -> Box<dyn Beta + Send + Sync> {
    let (beta_name, filename) = match beta_name.split_once(':') {
        Some(result) => result,
        None => {
            error!("Error: Unknown beta '{}'", beta_name);
            process::exit(1);
        }
    };

    match beta_name {
        TableBeta::NAME => {
            let table = table::Table::read_file(filename, false);
            Box::new(TableBeta::new(component, table, initial_ratio))
        },
        TableBeta1D::NAME => {
            Box::new(TableBeta1D::new(component.sideways, io::read_table(filename)))
        },
        _ => {
            error!("Error: Unknown beta '{}'", beta_name);
            process::exit(1);
        }
    }
}

pub struct BetaCite<'a> {
    pub name: String,
    pub summary: &'a str,
    pub cite: &'a str,
    pub args: &'a str,
    pub direction_of_interest: String,
}

/// Return a Vec of all the beta functions that are available.
/// Note: all semi-elliptical and quarter cracks are surface cracks.
/// The component geometry is required only for the compact tension beta.
pub fn get_all_betas() -> Vec<BetaCite<'static>> {
    
    vec![
        BetaCite {
            name: QuarterBroek86::NAME.to_owned(),
            summary: "quarter circular crack in an infinite plate in tension",
            cite: "[Broek 86]",
            args: "N/A",
            direction_of_interest: QuarterBroek86::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::NAME.to_owned(),
            summary: "semi-elliptical surface crack in a finite plate in tension",
            cite: "[Newman 86]",
            args: "f, s, aa, ca, as, cs",
            direction_of_interest: SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::NAME.to_owned(),
            summary: "semi-elliptical surface crack in an infinite plate in tension",
            cite: "[Anderson 05]",
            args: "aa, ca, as, cs",
            direction_of_interest: SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: QuarterEllipticalCornerCrackFinitePlateTensionNewman84::NAME.to_owned(),
            summary: "quarter elliptical corner crack in a finite plate in tension",
            cite: "[Newman 86]",
            args: "f, s, aa, ca, as, cs",
            direction_of_interest: QuarterEllipticalCornerCrackFinitePlateTensionNewman84::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: EllipticalEmbeddedCrackFinitePlateTensionNewman84::NAME.to_owned(),
            summary: "elliptical crack in a finite plate in tension",
            cite: "[Newman 86]",
            args: "f, s, aa, ca, as, cs",
            direction_of_interest: EllipticalEmbeddedCrackFinitePlateTensionNewman84::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SingleSidedEdgeCrackTensionTada73::NAME.to_owned(),
            summary: "single sided edge crack in a plate in tension",
            cite: "[Tada 73]",
            args: "s, cs",
            direction_of_interest: SingleSidedEdgeCrackTensionTada73::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: DoubleSidedEdgeCrackTensionTada73::NAME.to_owned(),
            summary: "double sided edge crack in a plate in tension",
            cite: "[Tada 73]",
            args: "s, cs",
            direction_of_interest: DoubleSidedEdgeCrackTensionTada73::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: DoubleSidedEdgeCrackTensionTada73_2::NAME.to_owned(),
            summary: "double sided edge crack in a plate in tension",
            cite: "[Tada 73]",
            args: "s, cs",
            direction_of_interest: DoubleSidedEdgeCrackTensionTada73_2::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: CompactCoupon::NAME.to_owned(),
            summary: "compact specimen in tension (scale is in load units not stress units) ",
            cite: "[ASTM 97]",
            args: "f, s, cs",
            direction_of_interest: CompactCoupon::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: CentreCrackTensionFedderson66::NAME.to_owned(),
            summary: "centre cracked plate in tension",
            cite: "[Fedderson 66]",
            args: "s, cs",
            direction_of_interest: CentreCrackTensionFedderson66::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: CentreCrackTensionKoiter65::NAME.to_owned(),
            summary: "centre cracked plate in tension",
            cite: "[Tada 73]",
            args: "s, cs",
            direction_of_interest: CentreCrackTensionKoiter65::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SingleSidedThroughCrackCircularHoleTensionBowie56::NAME.to_owned(),
            summary: "single sided through-thickness crack from a circular hole in tension",
            cite: "[Bowie 56]",
            args: "s, cs, r",
            direction_of_interest: SingleSidedThroughCrackCircularHoleTensionBowie56::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: DoubleSidedThroughCrackCircularHoleTensionBowie56::NAME.to_owned(),
            summary: "double sided through-thickness crack from a circular hole in tension",
            cite: "[Bowie 56]",
            args: "s, cs, r",
            direction_of_interest: DoubleSidedThroughCrackCircularHoleTensionBowie56::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiEllipticalSurfaceCrackRoundBarBendingMurakami87::NAME.to_owned(),
            summary: "semi-elliptical surface crack in a round bar in bending",
            cite: "[Murakami 87]",
            args: "f, as, cs",
            direction_of_interest: SemiEllipticalSurfaceCrackRoundBarBendingMurakami87::DEFAULT_DOI.to_string(),
        },       
        BetaCite {
            name: SemiEllipticalSurfaceCrackRoundBarTensionMurakami87::NAME.to_owned(),
            summary: "semi-elliptical surface crack in a round bar in tension",
            cite: "[Murakami 87]",
            args: "f, as, cs",
            direction_of_interest: SemiEllipticalSurfaceCrackRoundBarTensionMurakami87::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiEllipticalSurfaceCrackRoundBarBendingMurakami86::NAME.to_owned(),
            summary: "semi-elliptical surface crack in a round bar in bending",
            cite: "[Murakami 86]",
            args: "f, as, cs",
            direction_of_interest: SemiEllipticalSurfaceCrackRoundBarBendingMurakami86::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiEllipticalSurfaceCrackRoundBarTensionMurakami86::NAME.to_owned(),
            summary: "semi-elliptical surface crack in a round bar in tension",
            cite: "[Murakami 86]",
            args: "f, as, cs",
            direction_of_interest: SemiEllipticalSurfaceCrackRoundBarTensionMurakami86::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: EdgeCrackStripBendingMurakami87::NAME.to_owned(),
            summary: "edge crack in a strip in bending",
            cite: "[Murakami 87]",
            args: "s, cs",
            direction_of_interest: EdgeCrackStripBendingMurakami87::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: EdgeCrackStripTensionMurakami87::NAME.to_owned(),
            summary: "edge crack in a strip in tension",
            cite: "[Murakami 87]",
            args: "s, cs",
            direction_of_interest: EdgeCrackStripTensionMurakami87::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiCircularSurfaceCrackRoundBarTensionForman86::NAME.to_owned(),
            summary: "semi-circular surface crack in a round bar in tension",
            cite: "[Forman 86]",
            args: "f, as",
            direction_of_interest: SemiCircularSurfaceCrackRoundBarTensionForman86::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: SemiCircularSurfaceCrackRoundBarBendingForman86::NAME.to_owned(),
            summary: "semi-circular surface crack in a round bar in bending",
            cite: "[Forman 86]",
            args: "f, as",
            direction_of_interest: SemiCircularSurfaceCrackRoundBarBendingForman86::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: TableBeta1D::NAME.to_owned() + ":FILE",
            summary: "read FILE for 1-dimensional beta values",
            cite: "",
            args: "s, cs",
            direction_of_interest: TableBeta1D::DEFAULT_DOI.to_string(),
        },
        BetaCite {
            name: TableBeta::NAME.to_owned() + ":FILE",
            summary: "read FILE for 2-dimensional beta values",
            cite: "",
            args: "f, as, cs",
            direction_of_interest: TableBeta::DEFAULT_DOI.to_string(),
        },
    ]
}


#[derive(Clone)]
pub struct TableBeta {
    table: table::Table,
    direction_of_interest: DirectionOfInterest,
    b: f64,
    d: f64,
    ratio: f64,
}

impl TableBeta {
    pub const NAME: &'static str = "file-2d";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;

    pub fn new(component: &grow::Component, table: table::Table, initial_ratio: f64) -> Self {
        Self {
            table,
            direction_of_interest: TableBeta::DEFAULT_DOI,
            b: component.sideways,
            d: component.forward,
            ratio: initial_ratio,
        }
    }
}

impl Beta for TableBeta {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &self.direction_of_interest
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: PHI_A
            }),
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: Some(self.ratio),
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn set_d(&mut self, d: f64) {
        self.d = d;
    }

    fn beta(
        &self,
        geometry: &CrackGeometry
    ) -> Result<BetaResult, (BetaResult, String)> {
        let a = geometry.a.as_ref().unwrap().length;
        let a_on_d = a / self.d;
        let a_on_c = self.ratio;

        Ok(BetaResult{
            a: Some(self.table.interp(a_on_d, a_on_c)),
            c: None
        })
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        // This may not fit all 2d betas
        PI * a_on_d.powi(2) / 4.0
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

#[derive(Clone)]
pub struct TableBeta1D {
    data: Vec<Vec<f64>>,
    b: f64,
}

impl TableBeta1D {
    pub const NAME: &'static str = "file-1d";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;

    pub fn new(b: f64, data: Vec<Vec<f64>>) -> TableBeta1D {
        // Data should already be sorted, but do it anyway to be certain
        let mut data = data;
        data.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());

        TableBeta1D {
            data,
            b,
        }
    }

    fn interpolate(&self, x: f64) -> Result<f64, (f64, String)> {
        if x.is_nan() {
            return Err((self.data.first().unwrap()[1], "Invalid number used in interpolation, using smallest available value".to_owned()));
        }

        if x > self.data.last().unwrap()[0] {
            return Err((self.data.last().unwrap()[1], "Beta value outside range, using largest available value".to_owned()));
        }

        let index = match self.data.binary_search_by(|probe| probe[0].partial_cmp(&x).unwrap()) {
            // Found an exact match, return the exact value
            Ok(index) => return Ok(self.data[index][1]),
            Err(index) => index,
        };

        if index == 0 {
            return Err((self.data[index][1], "Beta value outside range, using smallest available value".to_owned()));
        }

        let x0 = self.data[index - 1][0].ln();
        let x1 = self.data[index][0].ln();
        let y0 = self.data[index - 1][1];
        let y1 = self.data[index][1];

        Ok(y0 + (x.ln() - x0) * ((y1 - y0) / (x1 - x0)))
    }
}

impl Beta for TableBeta1D {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &TableBeta1D::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        // todo!();
        let c = geometry.c.as_ref().unwrap().length;
        let c_on_b = c / self.b;

        match self.interpolate(c_on_b) {
            Ok(value) => { 
                Ok(BetaResult {
                    a: None,
                    c: Some(value),
                })
            },
            Err((value, message)) => {
                Err((BetaResult {
                    a: None,
                    c: Some(value),
                },
                message))
            }
        }
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

#[derive(Clone)]
struct QuarterBroek86 {}

impl QuarterBroek86 {
    pub const NAME: &'static str = "qct-broek86";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for QuarterBroek86 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &QuarterBroek86::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: PHI_A
            }),
            c: None,
            ratio: None,
        }
    }

    fn beta(&self, _geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let embedded_beta = 2.0 / PI;
        let result = 1.12 * 1.12 * embedded_beta;

        Ok(BetaResult {
            a: Some(result),
            c: None
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Beta factor for a corner crack in a plate in tension.
/// Plate is of thickness theta = 45 degrees.
///
/// Ref. Stress intensity factors handbook. Vol 2,
/// by Y. Murakami,
/// Pergamon Press, Oxford, , 1987.
// #[derive(Clone)]
// struct QuarterCircularCornerCrackFinitePlateTensionMurakami87 {
//     d: f64,
// }

// impl QuarterCircularCornerCrackFinitePlateTensionMurakami87 {
//     pub const NAME: &'static str = "qcft-murakami87";
// }

// impl Beta for QuarterCircularCornerCrackFinitePlateTensionMurakami87 {
//     fn direction_of_interest(&self) -> &DirectionOfInterest {
//         &DirectionOfInterest::A
//     }

//     fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
//         CrackGeometry {
//             a: Some(CrackFront {
//                 length: a,
//                 angle: PHI_A
//             }),
//             c: None,
//             ratio: None,
//         }
//     }

//     fn set_d(&mut self, d: f64) {
//         self.d = d;
//     }

//     fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
//         let a = geometry.a.as_ref().unwrap().length;
//         let a_on_d = a / self.d;

//         let result = ((1.05 + (-0.44 + 1.06 / 1.3) * a_on_d.powi(2) - 0.25 * a_on_d.powi(4))
//         * (1.0 + (0.08 + 0.40 * a_on_d.powi(2)) * (1.0 - (PI / 4.0).sin()).powi(3))
//         * (1.0 + (0.08 + 0.15 * a_on_d.powi(2)) * (1.0 - (PI / 4.0).cos()).powi(3))
//         * ((PI / 4.0).cos().powi(2) + (PI / 4.0).sin().powi(2)).powf(0.25)) / 2.464f64.sqrt();

//         Ok(BetaResult {
//             a: Some(result),
//             c: None,
//         })
//     }

//     fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
//         PI * a_on_d.powi(2) / 4.0
//     }

//     fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
//         Box::new(self.clone())
//     }
// }

/// Semi-elliptical beta factor for an edge crack.
/// Where *a* is the depth of crack, *2c*  is the surface length of the crack,
/// *`a_on_c` = a/c* and `phi` are the angle to a point on the crack front from surface .
///
/// Ref. Fracture Mechanics p.48,
/// by T. L. Anderson,
/// Taylor and Francis Group 2005.
///
/// Calculates the beta factor for a semi-elliptical crack in infinite
/// plate phi is the angle of the position on the crack front.  This
/// has been re-defined phi from the book so that it is from the
/// centreline i.e. phi = 0 is now at the deepest point and PI/2 is at
/// the surface.
#[derive(Clone)]
struct SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {
    a_angle: f64,
    c_angle: f64,
    doi: DirectionOfInterest,
}

impl SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {
    pub const NAME: &'static str = "seit-anderson05";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
    
    pub fn new(a_angle: f64, c_angle: f64, doi: Option<DirectionOfInterest>) -> Self {
        let doi = match doi {
            Some(value) => value,
            None => SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::DEFAULT_DOI,
        };
        
        Self {
            a_angle,
            c_angle,
            doi,
        }
    }

    fn beta(a_on_c: f64, phi: f64) -> f64 {
        let newman_phi = FRAC_PI_2 - phi;

        // surface correction factor
        let lambda_s = (1.13 - 0.09 * a_on_c) * (1.0 + 0.1 * (1.0 - newman_phi.sin()).powi(2));

        // angle correction factor
        let f_phi = (newman_phi.sin().powi(2) + (a_on_c * newman_phi.cos()).powi(2)).powf(0.25);
        let q = 1.0 + 1.464 * a_on_c.powf(1.65);

        lambda_s * f_phi / q.sqrt()
    }
}

impl Beta for SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &self.doi
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {   
        CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: self.a_angle,
            }),
            c: Some(CrackFront {
                length: c,
                angle: self.c_angle,
            }),
            ratio: None,
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let CrackFront{length: a, angle: phi_a} = geometry.a.as_ref().unwrap();
        let CrackFront{length: c, angle: phi_c} = geometry.c.as_ref().unwrap();
        let a_on_c = a / c;

        if a_on_c > 1.0 {
            println!(
                "Error: semi_elliptical_infinite_anderson05: invalid c > a. a/c = {}",
                a_on_c
            );
            process::exit(1);
        }

        let beta_a = SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::beta(a_on_c, *phi_a);
        let beta_c = SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::beta(a_on_c, *phi_c) * a_on_c.sqrt();

        Ok(BetaResult {
            a: Some(beta_a),
            c: Some(beta_c)
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Beta factor for a semi-elliptical crack under tension.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// NASA Technical Memorandum 85793.
#[derive(Clone)]
struct SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {
    d: f64,
    b: f64,
    a_angle: f64,
    c_angle: f64,
    doi: DirectionOfInterest
}

impl SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {
    pub const NAME: &'static str = "seft-newman84";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;

    pub fn new(b: f64, d: f64, a_angle: f64, c_angle: f64, doi: Option<DirectionOfInterest>) -> Self {
        let doi = match doi {
            Some(value) => value,
            None => SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::DEFAULT_DOI,
        };
        
        Self {
            b,
            d,
            a_angle,
            c_angle,
            doi,
        }
    }

    fn beta(a_on_c: f64, a_on_d: f64, c_on_a: f64, c_on_b: f64, phi: f64) -> f64 {
        let newman_phi = FRAC_PI_2 - phi;

        // finite width correction
        // sec x = 1/cos x
        let f_w = (FRAC_PI_2 * c_on_b * a_on_d.sqrt()).cos().recip().sqrt();

        let f = if a_on_c < 1.0 {
            let m1 = 1.13 - 0.09 * a_on_c;
            let m2 = -0.54 + 0.89 / (0.2 + a_on_c);
            let m3 = 0.5 - 1.0 / (0.65 + a_on_c) + 14.0 * (1.0 - a_on_c).powf(24.0);
            let g = 1.0 + (0.1 + 0.35 * a_on_d.powi(2)) * (1.0 - newman_phi.sin()).powi(2);
            debug!("m1 {} m2 {} m3 {} g {} f_w {} f_phi {}", m1, m2, m3, g, f_w, f_phi(a_on_c, newman_phi));

            (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powi(4)) * g * f_w
                * f_phi(a_on_c, newman_phi)
        } else {
            let m1 = c_on_a.sqrt() * (1.0 + 0.04 * c_on_a);
            let m2 = 0.2 * c_on_a.powf(4.0);
            let m3 = -0.11 * c_on_a.powf(4.0);

            // check on this it seems unsymmetrical
            // corrected the angle
            let g = 1.0
                + (0.1 + 0.35 * c_on_a * a_on_d.powi(2))
                    * (1.0 - newman_phi.sin()).powi(2);

            (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powi(4)) * g * f_w
                * f_phi(a_on_c, newman_phi)
        };

        f / shape_factor_newman84(a_on_c).sqrt()
    }
}

impl Beta for SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &self.doi
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: self.a_angle,
            }),
            c: Some(CrackFront {
                length: c,
                angle: self.c_angle,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn set_d(&mut self, d: f64) {
        self.d = d;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let CrackFront{length: a, angle: phi_a} = geometry.a.as_ref().unwrap();
        let CrackFront{length: c, angle: phi_c} = geometry.c.as_ref().unwrap();
        let a_on_c = a / c;
        let a_on_d = a / self.d;
        let c_on_a = c / a;
        let c_on_b = c / self.b;

        let beta_a = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, c_on_a, c_on_b, *phi_a);
        let beta_c = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, c_on_a, c_on_b, *phi_c) * a_on_c.sqrt();
        
        Ok(BetaResult {
            a: Some(beta_a),
            c: Some(beta_c)
        })
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        // This is half the crack and should be compared to b*d for failure check due to symmetry
        PI / 4.0 * a_on_d * c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Shape correction factor.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// Nasa Technical Memorandum 85793.

fn shape_factor_newman84(a_on_c: f64) -> f64 {
    if a_on_c < 1.0 {
        1.0 + 1.464 * a_on_c.powf(1.65)
    } else {
        1.0 + 1.464 * a_on_c.recip().powf(1.65)
    }
}

/// Beta factor for a quarter elliptical crack in tension.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// Nasa Technical Memorandum 85793.
#[derive(Clone)]
struct QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {
    d: f64,
    b: f64,
    a_angle: f64,
    c_angle: f64,
    doi: DirectionOfInterest,
}

impl QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {
    pub const NAME: &'static str = "qeft-newman84";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;

    fn new(b: f64, d: f64, a_angle: f64, c_angle: f64, doi: Option<DirectionOfInterest>) -> Self {
        let doi = match doi {
            Some(value) => value,
            None => QuarterEllipticalCornerCrackFinitePlateTensionNewman84::DEFAULT_DOI,
        };
        
        Self {
            b,
            d,
            a_angle,
            c_angle,
            doi,
        }
    }

    fn beta(a_on_c: f64, a_on_d: f64, c_on_a: f64, c_on_b: f64, c_on_d: f64, phi: f64) -> f64 {
        let newman_phi = FRAC_PI_2 - phi;

        // finite width correction
        // let f_w = (PI * c_on_b * a_on_d.sqrt() /2.0).cos().sqrt().recip(); // murakami
        let lambda = c_on_b * a_on_d.sqrt();
        let f_w = 1.0 - 0.2 * lambda + 9.4 * lambda.powi(2) - 19.4 * lambda.powi(3)
            + 27.1 * lambda.powi(4);

        let f = if a_on_c <= 1.0 {
            let m1 = 1.08 - 0.03 * a_on_c;
            let m2 = -0.44 + 1.06 / (0.3 + a_on_c);
            let m3 = -0.5 + 0.25 * a_on_c + 14.8 * (1.0 - a_on_c).powf(15.0);
            let g1 = 1.0 + (0.08 + 0.4 * a_on_d.powi(2)) * (1.0 - newman_phi.sin()).powi(3);
            let g2 = 1.0 + (0.08 + 0.15 * a_on_d.powi(2)) * (1.0 - newman_phi.cos()).powi(3);

            // println!("a/c <= 1 m1 {} m2 {} m3 {} g1 {} g2 {}", m1, m2, m3, g1, g2 );
            (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g1 * g2 * f_w
                * f_phi(a_on_c, newman_phi)
        } else {
            let m1 = c_on_a.sqrt() * (1.08 - 0.03 * c_on_a);
            let m2 = 0.375 * c_on_a.powi(2);
            let m3 = -0.25 * c_on_a.powi(2);

            // the report says c/t instead of c/b
            //            let g1 = 1.0 + (0.08 + 0.4 * c_on_b.powi(2)) * (1.0 - (FRAC_PI_2 - phi).sin()).powi(3);
            //            let g2 = 1.0 + (0.08 + 0.15 * c_on_b.powi(2)) * (1.0 - (FRAC_PI_2 - phi).cos()).powi(3);
            let g1 = 1.0 + (0.08 + 0.4 * c_on_d.powi(2)) * (1.0 - newman_phi.sin()).powi(3);
            let g2 = 1.0 + (0.08 + 0.15 * c_on_d.powi(2)) * (1.0 - newman_phi.cos()).powi(3);

            (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g1 * g2 * f_w
                * f_phi(a_on_c, newman_phi)
        };

        // println!("a_on_c {} f_w {} f {} shape factor {} beta {}", a_on_c, f_w, f, shape_factor(a_on_c).sqrt(), f/shape_factor(a_on_c).sqrt());
        f / shape_factor_newman84(a_on_c).sqrt()
    }
}

impl Beta for QuarterEllipticalCornerCrackFinitePlateTensionNewman84 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &self.doi
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: self.a_angle,
            }),
            c: Some(CrackFront {
                length: c,
                angle: self.c_angle,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn set_d(&mut self, d: f64) {
        self.d = d;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let CrackFront{length: a, angle: phi_a} = geometry.a.as_ref().unwrap();
        let CrackFront{length: c, angle: phi_c} = geometry.c.as_ref().unwrap();
        let a_on_d = a / self.d;
        let a_on_c = a / c;
        let c_on_a = c / a;
        let c_on_b = c / self.b;
        let c_on_d = c / self.d;

        let beta_a = QuarterEllipticalCornerCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, c_on_a, c_on_b, c_on_d, *phi_a);
        let beta_c = QuarterEllipticalCornerCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, c_on_a, c_on_b, c_on_d, *phi_c) * a_on_c.sqrt();

        Ok(BetaResult {
            a: Some(beta_a),
            c: Some(beta_c)
        })
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        PI * a_on_d * c_on_b / 4.0
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Angular correction `f_phi`.
fn f_phi(a_on_c: f64, newman_phi: f64) -> f64 {
    if a_on_c <= 1.0 {
        (a_on_c.powi(2) * newman_phi.cos().powi(2) + newman_phi.sin().powi(2)).powf(0.25)
    } else {
        (a_on_c.recip().powi(2) * newman_phi.sin().powi(2) + newman_phi.cos().powi(2)).powf(0.25)
    }
}

/// Shape correction factor.
///
/// Ref. Analyses of surface cracks in finite plates under tension or bending loads,
/// by J. C. Newman, Jr., and I. S. Raju,
/// April 1984,
/// NASA Technical Memorandum 85793.
#[derive(Clone)]
struct EllipticalEmbeddedCrackFinitePlateTensionNewman84 {
    d: f64,
    b: f64,
    a_angle: f64,
    c_angle: f64,
    doi: DirectionOfInterest,
}

impl EllipticalEmbeddedCrackFinitePlateTensionNewman84 {
    pub const NAME: &'static str = "eft-newman84";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;

    fn new(b: f64, d: f64, a_angle: f64, c_angle: f64, doi: Option<DirectionOfInterest>) -> Self {
        let doi = match doi {
            Some(value) => value,
            None => EllipticalEmbeddedCrackFinitePlateTensionNewman84::DEFAULT_DOI,
        };
        
        Self {
            b,
            d,
            a_angle,
            c_angle,
            doi,
        }
    }

    fn beta(a_on_c: f64, a_on_d: f64, c_on_a: f64, c_on_b: f64, phi: f64) -> f64 {
        let newman_phi = FRAC_PI_2 - phi;
        // finite width correction
        let f_w = (FRAC_PI_2 * c_on_b * a_on_d.sqrt()).cos().recip().sqrt();

        let m2 = 0.05 / (0.11 + a_on_c.powf(1.5));
        let m3 = 0.29 / (0.23 + a_on_c.powf(1.5));
        let g = 1.0
            - (a_on_d.powi(4) * (2.6 - 2.0 * a_on_d).sqrt() / (1.0 + 4.0 * a_on_c))
                * newman_phi.cos().abs();

        let m1 = if a_on_c <= 1.0 { 1.0 } else { c_on_a.sqrt() };

        let f = (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g * f_w
            * f_phi(a_on_c, newman_phi);

        f / shape_factor_newman84(a_on_c).sqrt()
    }
}

impl Beta for EllipticalEmbeddedCrackFinitePlateTensionNewman84 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &self.doi
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: self.a_angle,
            }),
            c: Some(CrackFront {
                length: c,
                angle: self.c_angle,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn set_d(&mut self, d: f64) {
        self.d = d;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let CrackFront{length: a, angle: phi_a} = geometry.a.as_ref().unwrap();
        let CrackFront{length: c, angle: phi_c} = geometry.c.as_ref().unwrap();
        let a_on_c = a / c;
        let a_on_d = a / self.d;
        let c_on_a = c / a;
        let c_on_b = c / self.b;
        
        let beta_a = EllipticalEmbeddedCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, c_on_a, c_on_b, *phi_a);
        let beta_c = EllipticalEmbeddedCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, c_on_a, c_on_b, *phi_c) * a_on_c.sqrt();

        Ok(BetaResult {
            a: Some(beta_a),
            c: Some(beta_c)
        })
    }

    fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        PI / 4.0 * a_on_d * c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Double corner crack in a hole in tension.
///
/// Taken from Y. Murakami
/// Stress Intensity Factors handbook
/// Vol. 2. P. 716 Pergamon Press, 1987.
/// In turn taken from J. C. Newman Jr. and I. S. Raju
/// Stress Intensity Factor Equations for Cracks
/// in three-dimensional finite bodies, Nasa Technical Memorandum 83299 1981 p 1--49.
// #[derive(Clone)]
// struct DoubleSidedCornerCrackHoleTensionNewman81 {
//     d: f64,
//     b: f64,
//     r: f64,
// }

// impl DoubleSidedCornerCrackHoleTensionNewman81 {
//     pub const NAME: &'static str = "dccht-newman81";

//     fn beta(a_on_c: f64, a_on_d: f64, c_on_a: f64, c_on_b: f64, c_on_r: f64, phi: f64) -> f64 {
//         let m1 = if a_on_c <= 1.0 {
//             1.13 - 0.09 * a_on_c
//         } else {
//             c_on_a.sqrt() * (1.0 + 0.04 * c_on_a)
//         };
//         let m2 = if a_on_c <= 1.0 {
//             -0.54 + (0.89 / (0.2 + a_on_c))
//         } else {
//             0.2 * c_on_a.powi(4)
//         };

//         let m3 = if a_on_c <= 1.0 {
//             0.5 - 1.0 / (0.65 + a_on_c)
//         } else {
//             -0.11 * c_on_a.powi(4)
//         };

//         let f_w = (FRAC_PI_2 * c_on_b * a_on_d.sqrt()).cos().recip().sqrt();

//         let newman_phi = FRAC_PI_2 - phi;

//         // finite width correction
//         let g1 = if a_on_c < 1.0 {
//             1.0 + (0.1 + 0.35 * c_on_a * a_on_d.powi(2)) * (1.0 - phi.sin()).powi(2)
//         } else {
//             1.0 + (0.1 + 0.35 * a_on_d.powi(2)) * (1.0 - phi.sin()).powi(2)
//         };
//         // println!("c_on_r {}", c_on_r);
//         let lambda = 1.0 / (1.0 + c_on_r * (0.85 * phi).cos());
//         let g2 = (1.0 - 0.15 * lambda + 3.46 * lambda.powi(2) - 4.47 * lambda.powi(3)
//             + 3.52 * lambda.powi(4)) / (1.0 + 0.13 * lambda.powi(2));

//         let g3 = if a_on_c <= 1.0 {
//             (1.0 + 0.04 * a_on_c) * (1.0 + 0.1 * (1.0 - phi.cos()).powi(2))
//                 * (0.8 + 0.2 * a_on_d.powf(0.25))
//         } else {
//             (1.13 - 0.09 * c_on_a) * (1.0 + 0.1 * (1.0 - phi.cos()).powi(2))
//                 * (0.8 + 0.2 * a_on_d.powf(0.25))
//         };

//         let f = (m1 + m2 * a_on_d.powi(2) + m3 * a_on_d.powf(4.0)) * g1 * g2 * g3 * f_w
//             * f_phi(a_on_c, newman_phi);
//         debug!("m1 {}, m2 {}, m3 {}, lamda {}, g1 {}, g2 {}, g3 {}, f_w {}, f {}", m1, m2, m3, lambda, g1, g2, g3, f_w, f);
//         f
//     }
// }

// impl Beta for DoubleSidedCornerCrackHoleTensionNewman81 {
//     fn direction_of_interest(&self) -> &DirectionOfInterest {
//         // TODO: This is an assumption, it could be either
//         &DirectionOfInterest::A
//     }

//     fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
//         CrackGeometry {
//             a: Some(CrackFront {
//                 length: a,
//                 angle: PHI_A
//             }),
//             c: Some(CrackFront {
//                 length: c,
//                 angle: PHI_C,
//             }),
//             ratio: None,
//         }
//     }

//     fn set_b(&mut self, b: f64) {
//         self.b = b;
//     }

//     fn set_d(&mut self, d: f64) {
//         self.d = d;
//     }

//     fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
//         let CrackFront{length: a, angle: phi_a} = geometry.a.as_ref().unwrap();
//         let CrackFront{length: c, angle: phi_c} = geometry.c.as_ref().unwrap();
//         let a_on_c = a / c;
//         let a_on_d = a / self.d;
//         let c_on_a = c / a;
//         let c_on_b = c / self.b;
//         let c_on_r = c / self.r;

//         let beta_a = DoubleSidedCornerCrackHoleTensionNewman81::beta(a_on_c, a_on_d, c_on_a, c_on_b, c_on_r, *phi_a);
//         let beta_c = DoubleSidedCornerCrackHoleTensionNewman81::beta(a_on_c, a_on_d, c_on_a, c_on_b, c_on_r, *phi_c);

//         Ok(BetaResult {
//             a: Some(beta_a),
//             c: Some(beta_c)
//         })
//     }

//     fn area(&self, a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
//         PI / 4.0 * a_on_d * c_on_b
//     }

//     fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
//         Box::new(self.clone())
//     }
// }

/// Single edge notched tension (SENT).
///
/// Ref. H. Tada, P.C. Paris and G. R. Irwin
/// The stress analysis of cracks handbook
/// P. 2.11
/// Compiled from Tada 1973
/// Uses the third equation which has accuracy better than 0.5% for any a/d
#[derive(Clone)]
struct SingleSidedEdgeCrackTensionTada73 {
    b: f64,
}

impl SingleSidedEdgeCrackTensionTada73 {
    pub const NAME: &'static str = "sset-tada73";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for SingleSidedEdgeCrackTensionTada73 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SingleSidedEdgeCrackTensionTada73::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: Some(a / c),
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.c.as_ref().unwrap().length;
        let limited_c_on_b = (c / self.b).min(0.8);
        let num = if limited_c_on_b == 0.0 {
            1.0
        } else {
            ((1.0 / (FRAC_PI_2 * limited_c_on_b)) * (FRAC_PI_2 * limited_c_on_b).tan()).sqrt()
        };

        let denom = (PI * limited_c_on_b / 2.0).cos();
        let result = (num / denom) * (0.752 + 2.02 * limited_c_on_b
            + 0.37 * (1.0 - (FRAC_PI_2 * limited_c_on_b).sin()).powi(3));

        Ok(BetaResult {
            a: None,
            c: Some(result),
        })
    }

    fn area(&self, _a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Middle cracked tension (centre cracked tension 2nd eq)
///
/// Fedderson 1966
/// Taken from Damage Tolerant Design handbook from Afgrow documentation.
/// <http://www.afgrow.net/applications/DTDHandbook/Sections/page11_3.aspx#standard_center_cracked_tension_specimen>
/// Note that `a_on_d` divided by 2 to account for the different definition of the width/depth in the equation.
//TODO: Confirm that the above condition is still required now that we've changed from a_on_d to c_on_b
#[derive(Clone)]
struct CentreCrackTensionFedderson66 {
    b: f64,
}

impl CentreCrackTensionFedderson66 {
    pub const NAME: &'static str = "ct-fedderson66";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for CentreCrackTensionFedderson66 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &CentreCrackTensionFedderson66::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.get_length_c().unwrap();
        let c_on_b = c / self.b;
        let f_0 = (PI * c_on_b / 2.0).cos().powf(-0.5);

        Ok(BetaResult {
            a: None,
            c: Some(f_0)
        })
    }

    fn area(&self, _a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Centre cracked tension.
///
/// Modification of Koiter 1965 from Tada 73
/// taken from H. Tada, P.C. Paris and G. R. Irwin (p 2.2)
//TODO: Confirm this is still accurate since changing from a_on_d to c_on_b
#[derive(Clone)]
struct CentreCrackTensionKoiter65 {
    b: f64,
}

impl CentreCrackTensionKoiter65 {
    pub const NAME: &'static str = "ct-koiter65";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for CentreCrackTensionKoiter65 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &CentreCrackTensionKoiter65::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: Some(a / c),
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.c.as_ref().unwrap().length;
        let c_on_b = c / self.b;
        let denom = (1.0 - c_on_b).sqrt();
        let numerator = 1.0 - 0.5 * c_on_b + 0.370 * c_on_b.powi(2) - 0.044 * c_on_b.powi(3);
        let result = numerator / denom;

        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn area(&self, _a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Double edge notched tension - DENT
///
/// Two edge cracks each of length a in a plate of width 2d
/// by H. Tada, P.C. Paris and G. R. Irwin.
#[derive(Clone)]
struct DoubleSidedEdgeCrackTensionTada73 {
    b: f64,
}

impl DoubleSidedEdgeCrackTensionTada73 {
    pub const NAME: &'static str = "dset-tada73";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for DoubleSidedEdgeCrackTensionTada73 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &DoubleSidedEdgeCrackTensionTada73::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.c.as_ref().unwrap().length;
        let c_on_b = c / self.b;
        let result = (1.122 - 0.561 * c_on_b - 0.205 * c_on_b.powi(2) + 0.471 * c_on_b.powi(3)
        - 0.190 * c_on_b.powi(4)) / (1.0 - c_on_b).sqrt();
        
        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn area(&self, _a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        // b = half of plate width because of symmetry
        c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


#[derive(Clone)]
struct DoubleSidedEdgeCrackTensionTada73_2 {
    b: f64,
}

impl DoubleSidedEdgeCrackTensionTada73_2 {
    pub const NAME: &'static str = "dset-tada73_2";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for DoubleSidedEdgeCrackTensionTada73_2 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &DoubleSidedEdgeCrackTensionTada73_2::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn set_b(&mut self, b: f64) {
        self.b = b;
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.c.as_ref().unwrap().length;
        let c_on_b = c / self.b;
        let result = (1.0 + 0.122 * (FRAC_PI_2 * c_on_b).cos().powi(4)) *
            (((2.0 * self.b) / (PI * c)) * (FRAC_PI_2 * c_on_b).tan()).sqrt();
        
        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn area(&self, _a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        // b = half of plate width because of symmetry
        c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Compact specimen width W, thickness B
///
/// by H. Tada, P.C. Paris and G. R. Irwin.
#[derive(Clone)]
struct CompactCoupon {
    b: f64,
    d: f64,
}

impl CompactCoupon {
    pub const NAME: &'static str = "compact";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for CompactCoupon {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &CompactCoupon::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: Some(a/c),
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.get_length_c().unwrap();
        let c_on_b = c / self.b;
        // This beta factor is non-dimensionalised with load
        // hence we have to divide by (PI * c).sqrt() so the will be cancelled out in the later calculation
        // this expression is valid for c/b > 0.2
        let correction = (PI * c).sqrt();
        let front = (self.d * self.b.sqrt()).recip();
        let scale = (2.0 + c_on_b) / (1.0 - c_on_b).powf(3.0 / 2.0);
        let rest = 0.886 + 4.64 * c_on_b - 13.32 * c_on_b.powi(2) + 14.72 * c_on_b.powi(3)
            - 5.60 * c_on_b.powi(4);

        let result = front * scale * rest / correction;

        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}



/// Semi-elliptical surface crack in a round bar in bending
/// Murakami87 p. 657
/// This is obtained by Murakami using the results for an edge crack in a strip
/// for which he has equations for both tension and bending whereas the crack in a
/// round bar is available for tension.
#[derive(Clone)]
struct SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {
    d: f64,
    ratio: f64,
}

impl SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {
    pub const NAME: &'static str = "serbb-murakami87";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SemiEllipticalSurfaceCrackRoundBarBendingMurakami87::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: None,
            ratio: Some(self.ratio),
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let semi_tension = SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 { d: self.d, ratio: self.ratio };
        let st = semi_tension.beta(geometry)?;

        // Need to use the a-direction from the main beta as a c-direction below
        let geometry_swapped = CrackGeometry {
            a: None,
            c: geometry.a.clone(),
            ratio: Some(self.ratio),
        };

        // Note using forward direction in place of sideways
        let edge_bending = EdgeCrackStripBendingMurakami87 { b: self.d };
        let eb = edge_bending.beta(&geometry_swapped)?;

        let edge_tension = EdgeCrackStripTensionMurakami87 { b: self.d };
        let et = edge_tension.beta(&geometry_swapped)?;
        
        let result =
            st.beta_of_interest(semi_tension.direction_of_interest()) * 
            eb.beta_of_interest(edge_bending.direction_of_interest()) / 
            et.beta_of_interest(edge_tension.direction_of_interest());

        Ok(BetaResult {
            a: Some(result),
            c: None
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Stress Intensity Factor Equations for a Semi-Elliptical Surface Crack in a shaft under bending
/// Yukitaka Murakami and Hideto Tsuru 86
/// No. 87-0164B
/// This is obtained by Murakami using the results for an edge crack in a strip
/// for which he has equations for both tension and bending whereas the crack in a
/// round bar is available for tension.
#[derive(Clone)]
struct SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {
    d: f64,
    ratio: f64,
}

impl SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {
    pub const NAME: &'static str = "serbb-murakami86";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SemiEllipticalSurfaceCrackRoundBarBendingMurakami86::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: None,
            ratio: Some(self.ratio),
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let semi_tension = SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 { d: self.d, ratio: self.ratio };
        let st = semi_tension.beta(geometry)?;

        // Need to use the a-direction from the main beta as a c-direction below
        let geometry_swapped = CrackGeometry {
            a: None,
            c: geometry.a.clone(),
            ratio: Some(self.ratio),
        };

        // Note using forward direction in place of sideways
        let edge_bending = EdgeCrackStripBendingMurakami87 { b: self.d };
        let eb = edge_bending.beta(&geometry_swapped)?;

        let edge_tension = EdgeCrackStripTensionMurakami87 { b: self.d };
        let et = edge_tension.beta(&geometry_swapped)?;
        
        let result =
            st.beta_of_interest(semi_tension.direction_of_interest()) * 
            eb.beta_of_interest(edge_bending.direction_of_interest()) / 
            et.beta_of_interest(edge_tension.direction_of_interest());

        Ok(BetaResult {
            a: Some(result),
            c: None
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Semi-elliptical surface crack in a round bar in tension
/// Murakami87 p. 657
#[derive(Clone)]
struct SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {
    d: f64,
    ratio: f64,
}

impl SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {
    pub const NAME: &'static str = "serbt-murakami87";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for SemiEllipticalSurfaceCrackRoundBarTensionMurakami87 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SemiEllipticalSurfaceCrackRoundBarTensionMurakami87::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: None,
            ratio: Some(self.ratio),
        }
    }

    // There seems to be an error somewhere as this equation does not
    // reproduce the table of sample points given by murakami87. It
    // comes close. The original reference provided is murakami and
    // Tsuru 86 which is also slightly different to murakami87 in the
    // second part of the equation and in the table of sample beta
    // factors for the round bar in bending.
    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let a = geometry.get_length_a().unwrap();
        let a_on_c = self.ratio; 
        let a_on_d = a / self.d;

        let result = (1.122 - 0.230 * a_on_c - 0.901 * a_on_c.powi(2) + 0.949 * a_on_c.powi(3)
                    - 0.280 * a_on_c.powi(4))
            * (1.0 + 0.157 * a_on_d - 0.634 * a_on_d.powi(2) + 4.590 * a_on_d.powi(3)
                   - 6.628 * a_on_d.powi(4));

        Ok(BetaResult{
            a: Some(result),
            c: None
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}



/// Semi-elliptical surface crack in a round bar in tension
/// Murakami and Tsuru 86
#[derive(Clone)]
struct SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {
    d: f64,
    ratio: f64,
}

impl SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {
    pub const NAME: &'static str = "serbt-murakami86";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for SemiEllipticalSurfaceCrackRoundBarTensionMurakami86 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SemiEllipticalSurfaceCrackRoundBarTensionMurakami86::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: None,
            ratio: Some(self.ratio),
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let a = geometry.get_length_a().unwrap();
        let a_on_c = self.ratio; 
        let a_on_d = a / self.d;

        // The original reference by
        // Stress Intensity Factor equations for a semi-elliptical surface crack in a shaft under bending
        let result = (1.122 - 0.230 * a_on_c - 0.901 * a_on_c.powi(2) + 0.949 * a_on_c.powi(3) - 0.280 * a_on_c.powi(4)) *
            (1.0 + 0.314 * a_on_d - 2.536 * a_on_d.powi(2) + 36.72 * a_on_d.powi(3) - 106.048 * a_on_d.powi(4));

        Ok(BetaResult { 
            a: Some(result),
            c: None
        })         
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}



/// Edge crack in a strip in bending
/// Murakami87 p. 657
#[derive(Clone)]
struct EdgeCrackStripBendingMurakami87 {
    b: f64,
}

impl EdgeCrackStripBendingMurakami87 {
    pub const NAME: &'static str = "esb-murakami87";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for EdgeCrackStripBendingMurakami87 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &EdgeCrackStripBendingMurakami87::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c_on_b = geometry.get_length_c().unwrap() / self.b;

        let result = 1.121 - 1.199 * c_on_b + 4.775 * c_on_b.powi(2) - 1.628 * c_on_b.powi(3)
            - 7.035 * c_on_b.powi(4) + 13.27 * c_on_b.powi(5);
        
        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Edge crack in a strup in bending tension
/// Murakami87 p. 657
#[derive(Clone)]
struct EdgeCrackStripTensionMurakami87 {
    b: f64,
}

impl EdgeCrackStripTensionMurakami87 {
    pub const NAME: &'static str = "est-murakami87";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for EdgeCrackStripTensionMurakami87 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &EdgeCrackStripTensionMurakami87::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c_on_b = geometry.get_length_c().unwrap() / self.b;

        let result = 1.12 - 0.231 * c_on_b + 10.55 * c_on_b.powi(2) - 21.72 * c_on_b.powi(3)
            + 30.39 * c_on_b.powi(4);

        Ok(BetaResult { 
            a: None,
            c: Some(result)
        })
    }

    fn area(&self, _a_on_d: f64, _a_on_c: f64, c_on_b: f64, _a_on_r: f64) -> f64 {
        c_on_b
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Beta for a circular corner/edge crack with the coupon constrained to extend uniaxially
///
/// The crack transitions from a corner to an edge crack
/// (Coupon is constrained so that ends are not free to rotate).
///
/// Here the a is in the depth direction of the coupon but it should
/// equal c The coupon is constrained so it is not free to rotate
/// which is more representative of a coupon clamped in the jaws in a
/// test machine. d is the thickness of the coupon.
/// This is probably for a 160 mm radius coupon but it is not clear.
///
/// Ref. Beta Values for Low Kt Specimens by M. `McDonald`
/// DSTO Minute Air07/048 Combat Aircraft Support.
//TODO: Confirm the accuracy of the implementation since the change to c_on_b
// #[derive(Clone)]
// struct CornerCrackConstrainedTensionMcdonald07 {
//     deep: table::Table,
//     b: f64,
// }

// impl CornerCrackConstrainedTensionMcdonald07 {
//     pub const NAME: &'static str = "qcct-mcdonald07";

//     fn new(b: f64) -> Self {
//         let c = [
//             0.0, 0.0001, 0.0006, 0.0011, 0.0016, 0.0021001, 0.0026003, 0.0031004, 0.0036005,
//             0.0041007, 0.0046014, 0.005102, 0.0056027, 0.0061034, 0.0066051, 0.0071055, 0.0076062,
//             0.0081089, 0.0086145, 0.0091147, 0.00962, 0.0101257,
//         ];
        
//         let betas = vec![
//             0.709411963f64,
//             0.709411963,
//             0.710173817,
//             0.714482118,
//             0.727905309,
//             0.750298769,
//             0.778871362,
//             0.808379238,
//             0.845057485,
//             0.892278507,
//             0.952439289,
//             1.013915829,
//             1.076584454,
//             1.134948417,
//             1.188004837,
//             1.240923778,
//             1.293874445,
//             1.347036619,
//             1.374464709,
//             1.405714819,
//             1.434011108,
//             1.46362585,
//         ];
        
//         // non-dimensionalise the crack depth data by the coupon width 25 mm
//         let c_on_bs = c.iter().map(|x| x / 25.0e-3).collect::<Vec<f64>>();
        
//         Self {
//             deep: table::Table::new(vec![0.0], c_on_bs, vec![betas], false),
//             b,
//         }
//     }
// }


// impl Beta for CornerCrackConstrainedTensionMcdonald07 {
//     fn direction_of_interest(&self) -> &DirectionOfInterest {
//         &DirectionOfInterest::C
//     }

//     fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
//         CrackGeometry {
//             a: None,
//             c: Some(CrackFront {
//                 length: c,
//                 angle: PHI_C,
//             }),
//             ratio: Some(a / c),
//         }
//     }

//     fn set_b(&mut self, b: f64) {
//         self.b = b;
//     }

//     fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
//         let c = geometry.c.as_ref().unwrap().length;
//         let c_on_b = c / self.b;
//         let result = self.deep.interp(c_on_b, 0.0);
        
//         Ok(BetaResult {
//             a: None,
//             c: Some(result)
//         })
//     }

//     fn area(&self, a_on_d: f64, _a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
//         1.0 - PI * a_on_d * a_on_d
//     }

//     fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
//         Box::new(self.clone())
//     }
// }

// Additional beta factors from
// @TechReport{paul88,
//   author = 	 {J. Paul and D. Lombardo},
//   title = 	 {CRKGRW - Crack growth program users manual},
//   institution =  {Defence Science and Technology Organisation},
//   year = 	 1988,
//   month = 	 {June}}

/// Bowie solution for a single crack from a circular hole
#[derive(Clone)]
struct SingleSidedThroughCrackCircularHoleTensionBowie56 {
    r: f64,
}

impl SingleSidedThroughCrackCircularHoleTensionBowie56 {
    pub const NAME: &'static str = "ssht-bowie56";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for SingleSidedThroughCrackCircularHoleTensionBowie56 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SingleSidedThroughCrackCircularHoleTensionBowie56::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.c.as_ref().unwrap().length;
        let c_on_r = c / self.r;
        let result = 0.6762 + 0.8734 / (0.3246 + c_on_r);
        
        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}

/// Bowie solution for a double cracked circular hole
#[derive(Clone)]
struct DoubleSidedThroughCrackCircularHoleTensionBowie56 {
    r: f64,
}

impl DoubleSidedThroughCrackCircularHoleTensionBowie56 {
    pub const NAME: &'static str = "dsht-bowie56";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::C;
}

impl Beta for DoubleSidedThroughCrackCircularHoleTensionBowie56 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &DoubleSidedThroughCrackCircularHoleTensionBowie56::DEFAULT_DOI
    }

    fn geometry(&self, _a: f64, c: f64) -> CrackGeometry {
        CrackGeometry {
            a: None,
            c: Some(CrackFront {
                length: c,
                angle: PHI_C,
            }),
            ratio: None,
        }
    }

    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let c = geometry.c.as_ref().unwrap().length;
        let c_on_r = c / self.r;
        let result = 0.9439 + 0.6865 / (0.2772 + c_on_r);

        Ok(BetaResult {
            a: None,
            c: Some(result)
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Semi-circular surface crack in a round bar in tension
#[derive(Clone)]
struct SemiCircularSurfaceCrackRoundBarTensionForman86 {
    d: f64,
}

impl SemiCircularSurfaceCrackRoundBarTensionForman86 {
    pub const NAME: &'static str = "scrbt-forman86";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for SemiCircularSurfaceCrackRoundBarTensionForman86 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SemiCircularSurfaceCrackRoundBarTensionForman86::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: None,
            ratio: None,
        }
    }

    #[allow(non_snake_case)]
    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let a = geometry.get_length_a().unwrap();
        let B = FRAC_PI_2 * (a / self.d);
        let Y = 1.0 - B.sin();
        let G = 0.92 * (2.0 / PI) * (1.0 / B.cos()) * (B.tan() / B).sqrt();

        let result = G * (0.752 + 1.286 * B + 0.37 * Y.powi(3));

        Ok(BetaResult{
            a: Some(result),
            c: None
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/// Semi-circular surface crack in a round bar in bending
#[derive(Clone)]
struct SemiCircularSurfaceCrackRoundBarBendingForman86 {
    d: f64,
}

impl SemiCircularSurfaceCrackRoundBarBendingForman86 {
    pub const NAME: &'static str = "scrbb-forman86";
    pub const DEFAULT_DOI: DirectionOfInterest = DirectionOfInterest::A;
}

impl Beta for SemiCircularSurfaceCrackRoundBarBendingForman86 {
    fn direction_of_interest(&self) -> &DirectionOfInterest {
        &SemiCircularSurfaceCrackRoundBarBendingForman86::DEFAULT_DOI
    }

    fn geometry(&self, a: f64, _c: f64) -> CrackGeometry {
        CrackGeometry {
            a: Some(
                CrackFront {
                    length: a,
                    angle: PHI_A,
                }
            ),
            c: None,
            ratio: None,
        }
    }

    #[allow(non_snake_case)]
    fn beta(&self, geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
        let a = geometry.get_length_a().unwrap();
        let B = FRAC_PI_2 * (a / self.d);
        let Y = 1.0 - B.sin();
        let G = 0.92 * (2.0 / PI) * (1.0 / B.cos()) * (B.tan() / B).sqrt();

        let result = G * (0.923 + 0.199 * Y.powi(4));

        Ok(BetaResult{
            a: Some(result),
            c: None
        })
    }

    fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
        Box::new(self.clone())
    }
}


/*
/// Another beta factor for a round bar with greater depth and hopefully more accurate than Murakami.
/// It provides the beta factor at two points around the crack front.
///
/// Experimental and finite element analyses on stress intensity
/// factors of an elliptical surface crack in a circular shaft under
/// tension and bending
/// C. S. Shin and C. Q. Cai
/// International Journal of Fracture 129: 239–264, 2004
struct SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
    deep: table::Table,
    surface: table::Table,
    d: f64,
}

impl SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
    fn new(d: f64) -> Self {
        let a_on_ds = vec![
            0.067, 0.133, 0.200, 0.267, 0.333, 0.400, 0.467, 0.533, 0.600, 0.667, 0.733, 0.800
        ];
        let a_on_cs = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

        // x/h = 0.0
        let betas_deep = vec![
            vec![0.963, 0.954, 0.929, 0.878, 0.834, 0.786, 0.739, 0.692, 0.649, 0.609, 0.576],
            vec![0.897, 0.890, 0.870, 0.840, 0.801, 0.757, 0.710, 0.662, 0.618, 0.576, 0.537],
            vec![0.872, 0.866, 0.848, 0.820, 0.783, 0.739, 0.690, 0.640, 0.592, 0.547, 0.506],
            vec![0.879, 0.873, 0.856, 0.828, 0.790, 0.743, 0.692, 0.637, 0.583, 0.532, 0.486],
            vec![0.917, 0.911, 0.893, 0.863, 0.823, 0.773, 0.716, 0.654, 0.592, 0.532, 0.478],
            vec![0.991, 0.984, 0.964, 0.932, 0.888, 0.832, 0.767, 0.695, 0.621, 0.549, 0.482],
            vec![1.112, 1.104, 1.082, 1.045, 0.994, 0.930, 0.854, 0.768, 0.678, 0.588, 0.504],
            vec![1.302, 1.294, 1.268, 1.224, 1.164, 1.087, 0.995, 0.889, 0.775, 0.659, 0.550],
            vec![1.609, 1.599, 1.566, 1.512, 1.437, 1.341, 1.224, 1.088, 0.938, 0.783, 0.634],
            vec![2.126, 2.113, 2.070, 1.998, 1.899, 1.771, 1.614, 1.429, 1.222, 1.002, 0.787],
            vec![3.082, 3.063, 3.002, 2.899, 2.755, 2.570, 2.342, 2.069, 1.758, 1.421, 1.083],
            vec![5.140, 5.110, 5.011, 4.841, 4.606, 4.302, 3.923, 3.466, 2.934, 2.344, 1.737],
        ];

        // x/h = 1.0
        let betas_surface = vec![
            vec![0.486, 0.523, 0.553, 0.578, 0.596, 0.609, 0.616, 0.618, 0.613, 0.603, 0.587],
            vec![0.510, 0.548, 0.579, 0.604, 0.623, 0.635, 0.641, 0.640, 0.633, 0.619, 0.599],
            vec![0.557, 0.596, 0.629, 0.654, 0.673, 0.684, 0.689, 0.686, 0.677, 0.660, 0.637],
            vec![0.600, 0.640, 0.673, 0.699, 0.717, 0.728, 0.732, 0.728, 0.717, 0.699, 0.674],
            vec![0.654, 0.695, 0.729, 0.755, 0.774, 0.784, 0.788, 0.783, 0.771, 0.751, 0.724],
            vec![0.742, 0.786, 0.822, 0.850, 0.869, 0.880, 0.882, 0.877, 0.862, 0.840, 0.809],
            vec![0.877, 0.926, 0.966, 0.996, 1.017, 1.028, 1.029, 1.022, 1.004, 0.978, 0.941],
            vec![1.062, 1.118, 1.163, 1.197, 1.220, 1.231, 1.232, 1.222, 1.200, 1.168, 1.124],
            vec![1.326, 1.393, 1.446, 1.485, 1.511, 1.524, 1.524, 1.510, 1.483, 1.443, 1.389],
            vec![1.755, 1.838, 1.904, 1.953, 1.985, 2.000, 1.998, 1.979, 1.943, 1.890, 1.820],
            vec![2.544, 2.659, 2.749, 2.816, 2.859, 2.878, 2.873, 2.844, 2.791, 2.714, 2.614],
            vec![4.138, 4.317, 4.458, 4.560, 4.625, 4.651, 4.639, 4.589, 4.500, 4.374, 4.209],
        ];

        let deep = table::Table::new(
            a_on_cs.clone(),
            a_on_ds.clone(),
            table::transpose(&betas_deep),
            false,
        );
        let surface = table::Table::new(
            a_on_cs,
            a_on_ds,
            table::transpose(&betas_surface),
            false,
        );

        SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
            deep,
            surface,
            d,
        }
    }
}

impl Beta for SemiEllipticalSurfaceCrackRoundBarBendingShin04 {
    fn beta(
        &self,
        a: f64,
        c: f64,
        _phis: &[f64],
    ) -> Result<Vec<f64>, (Vec<f64>, String)> {
        let a_on_c = a / c;
        let a_on_d = a / self.d;
        Ok(vec![
            self.deep.interp(a_on_d, a_on_c),
            self.surface.interp(a_on_d, a_on_c),
        ])
    }

    fn area(&self, a_on_d: f64, a_on_c: f64, _c_on_b: f64, _a_on_r: f64) -> f64 {
        2.0 * a_on_d * a_on_c
    }
}

*/

pub mod testing_objects {
    use super::*;

    struct Constant {
        value: f64
    }

    impl Constant {
        pub fn new(value: f64) -> Self {
            Self {
                value,
            }
        }
    }

    impl Beta for Constant {
        fn geometry(&self, a: f64, c: f64) -> CrackGeometry {
            CrackGeometry {
                a: Some(CrackFront {
                    length: a,
                    angle: PHI_A
                }),
                c: Some(CrackFront {
                    length: c,
                    angle: PHI_C,
                }),
                ratio: None,
            }
        }

        fn beta(&self, _geometry: &CrackGeometry) -> Result<BetaResult, (BetaResult, String)> {
            Ok(BetaResult {
                a: Some(self.value),
                c: Some(self.value),
            })
        }

        fn direction_of_interest(&self) -> &DirectionOfInterest {
            todo!()
        }

        fn inner_clone(&self) -> Box<dyn Beta + Send + Sync> {
            todo!()
        }
    }

    /// Get a boxed beta which will return a constant value for both
    /// 'a' and 'c' from beta()
    /// 
    /// Currently only viable for testing in non-optimisation scenarios
    pub fn make_constant(value:f64) -> Box<dyn Beta + Send + Sync> {
        Box::new(Constant::new(value)) as Box<dyn Beta + Send + Sync>
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_beta_1d_interpolate() {
        let data = 
            vec![
                vec![0.0, 0.0],
                vec![1.0, 1.0],
                vec![2.0, 2.0],
                vec![3.0, 3.0],
                vec![4.0, 4.0],
            ];

        let table_beta_1d = TableBeta1D::new(2.0, data.clone());
        match table_beta_1d.interpolate(2.5) {
            Ok(result) => println!("{}", result),
            Err(_) => assert!(false),
        }
    }

    #[test]
    fn test_centre_crack_tension_fedderson66() {
        let centre = CentreCrackTensionFedderson66 { b: 1.0 };
        let geometry = centre.geometry(0.0, 0.0);
        let result = centre.beta(&geometry).unwrap().c.unwrap();

        assert!(result - 1.0 < std::f64::EPSILON);
    }

    #[test]
    fn test_centre_crack_tension_koiter65() {
        let centre = CentreCrackTensionKoiter65 { b: 1.0 };
        let geometry = centre.geometry(0.0, 0.0);
        let result = centre.beta(&geometry).unwrap().c.unwrap();

        assert!(result - 1.0 < std::f64::EPSILON);
    }

    // #[test]
    // fn test_mcdonald07() {
    //     let corner = CornerCrackConstrainedTensionMcdonald07::new(1.0);

    //     let mut geometry = corner.geometry(0.0, 0.0);
    //     let mut result = corner.beta(&geometry).unwrap().c.unwrap();
    //     assert!((result - 0.709411963).abs() < std::f64::EPSILON);

    //     geometry.set_length_c(0.00010 / 25.0e-3);
    //     result = corner.beta(&geometry).unwrap().c.unwrap();
    //     assert!((result - 0.709411963).abs() < std::f64::EPSILON);

    //     geometry.set_length_c(0.00987285 / 25.0e-3);
    //     result = corner.beta(&geometry).unwrap().c.unwrap();
    //     assert!((result - 1.4485189776434797).abs() < 1e-3);
    // }
 
    #[test]
    // check that the if statement arms a_on_c <> 1 produce the same answer
    fn check_quarter_arms() {
        let offset = 1e-6;

        let lower = QuarterEllipticalCornerCrackFinitePlateTensionNewman84::beta(1.0 - offset, 0.0, 1.0 / (1.0 - offset), 0.0, 0.0, 0.0);
        let upper = QuarterEllipticalCornerCrackFinitePlateTensionNewman84::beta(1.0 + offset, 0.0, 1.0 / (1.0 + offset), 0.0, 0.0, 0.0);

        assert!((lower - upper).abs() < 1e-6);

        let lower = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(1.0 - offset, 0.0, 1.0 / (1.0 - offset), 0.0, 0.0);
        let upper = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(1.0 + offset, 0.0, 1.0 / (1.0 + offset), 0.0, 0.0);

        assert!((lower - upper).abs() < 1e-6);
    }

    #[test]
    // Python
    // np.sqrt(((2/np.pi) / a_on_b) * np.tan(a_on_b * np.pi/2)) *
    // (0.752 + 2.02 * a_on_b + 0.37 * (1 - np.sin(a_on_b * np.pi/2))**3) / np.cos(a_on_b * np.pi/2))
    fn check_single_edge_notched_tension_tada73() {
        let single = SingleSidedEdgeCrackTensionTada73 { b: 1.0 };
        let geometry = single.geometry(0.0, 0.0);
        let result = single.beta(&geometry).unwrap().c.unwrap();

        println!("tada73 sent(0.0) = {}", result);

        assert!((result - 1.122).abs() < 1e-6);
    }

    // #[test]
    // fn check_finite_quarter_beta() {
    //     let quarter = QuarterCircularCornerCrackFinitePlateTensionMurakami87 { d: 1.0 };

    //     let mut geometry = quarter.geometry(0.0, 0.0);
    //     let mut result = quarter.beta(&geometry).unwrap().a.unwrap();
    //     assert!((result - 0.671604).abs() < 1e-6);

    //     geometry.set_length_a(0.5);
    //     result = quarter.beta(&geometry).unwrap().a.unwrap();
    //     assert!((result - 0.724126).abs() < 1e-6);

    //     geometry.set_length_a(1.0);
    //     result = quarter.beta(&geometry).unwrap().a.unwrap();
    //     assert!((result - 0.7622).abs() < 1e-6);
    // }

    #[test]
    fn check_semi_elliptical_infinite_anderson05() {
        let mut result = SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::beta(1.0, 0.0);
        assert!((result - 0.662541348868913).abs() < std::f64::EPSILON);

        result = SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05::beta(0.5, 0.0);
        assert!((result - 0.8959634744787476).abs() < std::f64::EPSILON);
    }

    #[test]
    fn check_semi_elliptical_newman84() {
        let mut result = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(1.0, 0.0, 1.0, 0.0, 0.0);
        assert!(result - 0.6625413488689131 < std::f64::EPSILON);

        result = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(0.5, 0.0, 2.0, 0.0, 0.0);
        assert!(result - 0.8959634744787476 < std::f64::EPSILON);
    }
/*
    #[test]
    fn print_results() {
        let phis = vec![0.0, FRAC_PI_2];
        let semi = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84 { d: todo!(), b: todo!() };
        println!(
            "SemiEllipticalSurfaceCrackFinitePlateTensionNewman84(1.0, [0, pi/2]) {:?}",
            semi.beta(1.0, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SemiEllipticalSurfaceCrackFinitePlateTensionNewman84(0.5, [0, pi/2]) {:?}",
            semi.beta(0.5, 0.0, 0.0, 0.0, &phis)
        );

        // equal a and c betas
        println!(
            "SemiEllipticalSurfaceCrackFinitePlateTensionNewman84(0.82, [0, pi/2]) {:?}",
            semi.beta(0.82, 0.0, 0.0, 0.0, &phis)
        );

        println!(
            "SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05(1.0, [0, pi/2]) {:?}",
            semi.beta(1.0, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SemiEllipticalSurfaceCrackInfinitePlateTensionAnderson05(0.5, [0, pi/2]) {:?}",
            semi.beta(0.5, 0.0, 0.0, 0.0, &phis)
        );

        let quarter = QuarterBroek86 {};
        println!("Quarterbroek86 {:?}", quarter.beta(0.0, 0.0, 0.0, 0.0, &phis));

        let quarter = QuarterCircularCornerCrackFinitePlateTensionMurakami87 { d: todo!() };
        println!(
            "QuarterCircularCornerCrackFinitePlateTensionMurakami87(a_on_t = 0) {:?}",
            quarter.beta(0.0, 0.0, 0.0, 0.0, &phis)
        );
        println!("QuarterEllipticalCornerCrackFinitePlateTensionNewman84(a_on_c=1.0, a_on_d=0.0, c_on_b=0, [0, pi/2]) {:?}", quarter.beta(1.0, 0.0, 0.0, 0.0, &phis));
        println!("QuarterEllipticalCornerCrackFinitePlateTensionNewman84(a_on_c=0.5, a_on_d=0.0, c_on_b=0, [0, pi/2]) {:?}",   quarter.beta(0.5, 0.0, 0.0, 0.0, &phis));

        let elliptical = EllipticalEmbeddedCrackFinitePlateTensionNewman84 { d: todo!(), b: todo!() };
        println!("EllipticalEmbeddedCrackFinitePlateTensionNewman84(a_on_c=1.0, a_on_d=0.0, c_on_b=0, [0, pi/2]) {:?}",  elliptical.beta(1.0, 0.0, 0.0, 0.0, &phis));

        let single = SingleSidedEdgeCrackTensionTada73 { b: todo!() };
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 0) {:?}",
            single.beta(0.0, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 1e-6) {:?}",
            single.beta(1e-6, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 0.1) {:?}",
            single.beta(0.1, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 0.5) {:?}",
            single.beta(0.5, 0.0, 0.0, 0.0, &phis)
        );
        println!(
            "SingleSidedEdgeCrackTensionTada73(a_on_d = 1.2) {:?}",
            single.beta(1.2, 0.0, 0.0, 0.0, &phis)
        );

        let double = DoubleSidedEdgeCrackTensionTada73 { b: todo!() };
        println!(
            "DoubleSidedEdgeCrackTensionTada73(a_on_d = 0) {:?}",
            double.beta(0.0, 0.0, 0.0, 0.0, &phis)
        );

        let compact_tension_tada73 = CompactCoupon {
            width: 1.0,
            depth: 1.0,
        };
        println!(
            "compact_tension_tada73(a_on_d = 0.5, b=1, w = 1) {:?}",
            compact_tension_tada73.beta(0.5, 1.0, 1.0, 0.0, &phis)
        );
    }

    #[test]
    fn check_semi_elliptical_surface_crack_round_bar_bending_shin04() {
        let phis = vec![0.0, FRAC_PI_2];
        let semi = SemiEllipticalSurfaceCrackRoundBarBendingShin04::new(todo!());

        let result = semi.beta(0.067, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.963).abs() < 1e-3);

        let result = semi.beta(0.2, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.872).abs() < 1e-3);

        let result = semi.beta(0.5, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 1.19672).abs() < 1e-3);

        let result = semi.beta(0.8, 0.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 5.140).abs() < 1e-3);

        let result = semi.beta(0.067, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.576).abs() < 1e-3);

        let result = semi.beta(0.2, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.506).abs() < 1e-3);

        let result = semi.beta(0.5, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.524).abs() < 1e-3);

        let result = semi.beta(0.8, 1.0, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 1.737).abs() < 1e-3);

        let result = semi.beta(0.067, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.81).abs() < 1e-3);

        let result = semi.beta(0.2, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 0.761).abs() < 1e-3);

        let result = semi.beta(0.5, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 1.0348).abs() < 1e-3);

        let result = semi.beta(0.8, 0.45, 0.0, 0.0, &phis)[0];
        println!("Result: {}", result);
        assert!((result - 4.454).abs() < 1e-3);
    }

    #[test]
    fn manual_check_print_round_bar_table_murakami87() {
        // prints the beta table from murakami P. 658.
        // The first row of the table for a/d = 0 is correct but there is a slight error thereafter.
        // This implies there is an error with the a/d effect (a/d > 0)

        // This is the generated version of the beta table using the code
        // below. The results in this table are slightly different from the
        // that provided on P.658 of Murakami. Since this beta is based on the
        // following three betas:
        //
        // semi_elliptical_surface_crack_round_bar_tension_murakami87(*a_on_d, *crack.a_on_c);
        // edge_crack_strip_bending_murakami87(*a_on_d);
        // edge_crack_strip_tension_murakami87(*a_on_d);
        //
        // then there is either a mistake in my implementation of one or more of
        // these three or in the table that was produced by Murakami.

        // However, the formula given by murakami must be in error since the
        // bending of a crack must start with term of 1.12, but this is for the
        // equation containing the crack aspect ratio beta = b/a which is the
        // aspect ratio of the crack.

        //                                           a/c
        // a/d  | 0.000  0.100  0.200  0.300  0.400  0.500  0.600  0.700  0.800  0.900  1.000
        //      -----------------------------------------------------------------------------
        // 0    | 1.123  1.092  1.048  0.996  0.940  0.884  0.829  0.778  0.733  0.694  0.661
        // 0.01 | 1.114  1.084  1.040  0.989  0.933  0.877  0.823  0.773  0.728  0.689  0.656
        // 0.02 | 1.105  1.074  1.031  0.980  0.925  0.869  0.816  0.766  0.721  0.683  0.650
        // 0.03 | 1.094  1.064  1.021  0.970  0.916  0.861  0.808  0.758  0.714  0.676  0.643
        // 0.04 | 1.082  1.052  1.010  0.960  0.906  0.852  0.799  0.750  0.707  0.669  0.637
        // 0.05 | 1.070  1.041  0.999  0.949  0.896  0.842  0.790  0.742  0.699  0.661  0.630
        // 0.06 | 1.058  1.028  0.987  0.938  0.885  0.832  0.781  0.733  0.691  0.654  0.622
        // 0.07 | 1.045  1.016  0.975  0.927  0.875  0.822  0.771  0.724  0.682  0.646  0.615
        // 0.08 | 1.032  1.003  0.963  0.915  0.864  0.812  0.762  0.715  0.674  0.637  0.607
        // 0.09 | 1.018  0.990  0.950  0.903  0.853  0.801  0.752  0.706  0.665  0.629  0.599
        // 0.1  | 1.005  0.977  0.938  0.892  0.842  0.791  0.742  0.697  0.656  0.621  0.591
        // 0.11 | 0.992  0.965  0.926  0.880  0.831  0.781  0.733  0.688  0.648  0.613  0.584
        // 0.12 | 0.979  0.952  0.914  0.869  0.820  0.771  0.723  0.679  0.640  0.605  0.576
        // 0.13 | 0.967  0.940  0.902  0.858  0.810  0.761  0.714  0.670  0.631  0.597  0.569
        // 0.14 | 0.955  0.928  0.891  0.847  0.799  0.751  0.705  0.662  0.623  0.590  0.562
        // 0.15 | 0.943  0.917  0.880  0.836  0.789  0.742  0.696  0.654  0.616  0.583  0.555
        // 0.16 | 0.931  0.905  0.869  0.826  0.780  0.733  0.688  0.646  0.608  0.575  0.548
        // 0.17 | 0.920  0.895  0.859  0.816  0.770  0.724  0.679  0.638  0.601  0.569  0.541
        // 0.18 | 0.909  0.884  0.849  0.807  0.761  0.716  0.671  0.630  0.594  0.562  0.535
        // 0.19 | 0.899  0.874  0.839  0.797  0.753  0.707  0.664  0.623  0.587  0.556  0.529
        // 0.2  | 0.889  0.864  0.830  0.789  0.744  0.700  0.656  0.616  0.580  0.549  0.523
        // 0.21 | 0.879  0.855  0.821  0.780  0.736  0.692  0.649  0.610  0.574  0.543  0.517
        // 0.22 | 0.870  0.846  0.812  0.772  0.728  0.685  0.642  0.603  0.568  0.538  0.512
        // 0.23 | 0.861  0.837  0.804  0.764  0.721  0.677  0.636  0.597  0.562  0.532  0.506
        // 0.24 | 0.852  0.829  0.795  0.756  0.714  0.671  0.629  0.591  0.556  0.527  0.501
        // 0.25 | 0.844  0.820  0.787  0.748  0.706  0.664  0.623  0.585  0.551  0.521  0.496
        
        let table_beta = SemiEllipticalSurfaceCrackRoundBarBendingMurakami87 { d: todo!() };
        println!("{}", table_beta.as_table());
    }

    #[test]
    fn manual_check_print_round_bar_table_murakami86() {
        // prints the beta table from Murakami and Tsuru 86
        let table_beta = SemiEllipticalSurfaceCrackRoundBarBendingMurakami86 { d: todo!() };
        println!("{}", table_beta.as_table());
    }
*/
    #[test]
    // check that for the same conditions the quarter elliptical
    // always produces a higher beta than the semi elliptical
    fn check_semi_quarter_elliptical_crack_newman84() {
        let a_on_ds = [0.0f64, 0.2, 0.5, 0.8];
        let a_on_cs = [0.2f64, 0.5, 0.8];
        let phi = 0.0;
        let c_on_b = 0.0;

        for &a_on_d in &a_on_ds {
            for &a_on_c in &a_on_cs {
                println!("a_on_d {}, a_on_c {}", a_on_d, a_on_c);

                let s = SemiEllipticalSurfaceCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, 1.0 / a_on_c, c_on_b, phi);
                let q = QuarterEllipticalCornerCrackFinitePlateTensionNewman84::beta(a_on_c, a_on_d, 1.0 / a_on_c, c_on_b, 0.0, phi);

                println!("semi {:?}, quarter {:?}", s, q);
                assert!(q > s);
            }
        }
    }
}
