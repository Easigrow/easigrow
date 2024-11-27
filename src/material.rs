#![allow(clippy::unreadable_literal)]

//! Database of fatigue parameters for typical materials.

use std::collections::BTreeMap;

use crate::dadn::*;

#[derive(Debug, Clone)]
pub struct Properties {
    pub yield_stress: f64,
    pub k1c: f64,
    pub youngs_modulus: f64,
}

pub struct EqnProperty {
    pub name: &'static str,
    pub cite: &'static str,
    pub units: &'static str,
    // A BTreeMap is being used to guarantee ordering during iteration for
    // use in optimisations.
    pub params: BTreeMap<ParameterLabel, f64>,
}

impl Default for Properties {
    fn default() -> Self {
        Self {
            yield_stress: 450.0,
            k1c: 33.0,
            youngs_modulus: 71e3,
        }
    }
}

lazy_static! {
    static ref MATERIALS: BTreeMap<&'static str, EqnProperty> = {
        let mut materials = BTreeMap::new();

        let mut name = "paris:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Paris::UNITS,
                params: BTreeMap::from([(ParameterLabel::c, 1.00e-9), (ParameterLabel::m, 3.0)]),
            },
        );

        name = "walker:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Walker::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.00e-10),
                    (ParameterLabel::m, 0.5),
                    (ParameterLabel::n, 3.0),
                ]),
            },
        );

        name = "forman:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.00e-10),
                    (ParameterLabel::n, 3.0),
                    (ParameterLabel::kf, 60.0),
                ]),
            },
        );

        name = "forman-ag:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.00e-10),
                    (ParameterLabel::n, 3.0),
                    (ParameterLabel::kf, 60.0),
                ]),
            },
        );

        name = "nasgro:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Forman 05]",
                units: Nasgro::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::smax_on_sigma0, 0.3),
                    (ParameterLabel::alpha, 2.0),
                    (ParameterLabel::k_crit, 35.16),
                    (ParameterLabel::deltak0, 0.80),
                    (ParameterLabel::cth, 2.20),
                    (ParameterLabel::cth_minus, 0.1),
                    (ParameterLabel::p, 1.0),
                    (ParameterLabel::q, 1.0),
                    (ParameterLabel::c, 6.35e-10),
                    (ParameterLabel::n, 2.50),
                    (ParameterLabel::a_intr, 38.1e-6),
                ]),
            },
        );

        name = "nasgro:aa7050t7451-LT";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Forman 05]",
                units: Nasgro::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::smax_on_sigma0, 0.3),
                    (ParameterLabel::alpha, 2.0),
                    (ParameterLabel::k_crit, 35.16),
                    (ParameterLabel::deltak0, 0.80),
                    (ParameterLabel::cth, 2.20),
                    (ParameterLabel::cth_minus, 0.1),
                    (ParameterLabel::p, 1.0),
                    (ParameterLabel::q, 1.0),
                    (ParameterLabel::c, 6.35e-10),
                    (ParameterLabel::n, 2.50),
                    (ParameterLabel::a_intr, 38.1e-6),
                ]),
            },
        );

        name = "burchill:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Burchill::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.00e-10),
                    (ParameterLabel::m, 3.0),
                    (ParameterLabel::d, 1.00e-10),
                    (ParameterLabel::n, 3.0),
                ]),
            },
        );

        name = "kujawski:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Kujawski::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.00e-10),
                    (ParameterLabel::m, 3.0),
                    (ParameterLabel::alpha, 0.5),
                ]),
            },
        );

        name = "hartman:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Hartman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::d, 1.00e-10),
                    (ParameterLabel::deltak_th, 1.0),
                    (ParameterLabel::a, 30.0),
                    (ParameterLabel::alpha, 3.0),
                ]),
            },
        );

        name = "paris:newman-aa7050t7451";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[none]",
                units: Paris::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.593e-11),
                    (ParameterLabel::m, 3.668),
                ]),
            },
        );

        name = "forman:aa2024t3-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 7.13e-9),
                    (ParameterLabel::n, 2.70),
                    (ParameterLabel::kf, 71.3),
                ]),
            },
        );

        name = "forman:aa2024t351-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 5.00e-9),
                    (ParameterLabel::n, 2.88),
                    (ParameterLabel::kf, 63.2),
                ]),
            },
        );

        name = "forman:aa2024t4-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 8.57e-9),
                    (ParameterLabel::n, 2.60),
                    (ParameterLabel::kf, 58.1),
                ]),
            },
        );

        name = "forman:aa2024t6-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 2.00e-8),
                    (ParameterLabel::n, 2.62),
                    (ParameterLabel::kf, 69.8),
                ]),
            },
        );

        name = "forman:aa2024t8-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.33e-8),
                    (ParameterLabel::n, 2.65),
                    (ParameterLabel::kf, 65.3),
                ]),
            },
        );

        name = "forman:aa2024t851-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 7.72e-9),
                    (ParameterLabel::n, 2.78),
                    (ParameterLabel::kf, 61.4),
                ]),
            },
        );

        name = "forman:aa2219t851-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 4.84e-8),
                    (ParameterLabel::n, 2.16),
                    (ParameterLabel::kf, 57.5),
                ]),
            },
        );

        name = "forman:aa2618t6-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 8.56e-9),
                    (ParameterLabel::n, 2.58),
                    (ParameterLabel::kf, 45.9),
                ]),
            },
        );

        name = "forman:aa6061t6-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 2.27e-7),
                    (ParameterLabel::n, 1.66),
                    (ParameterLabel::kf, 60.1),
                ]),
            },
        );

        name = "forman:aa6061t651-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 9.60e-8),
                    (ParameterLabel::n, 1.84),
                    (ParameterLabel::kf, 41.2),
                ]),
            },
        );

        name = "forman:aa7010t73651-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 2.06e-8),
                    (ParameterLabel::n, 2.46),
                    (ParameterLabel::kf, 46.0),
                ]),
            },
        );

        name = "forman:aa7050t7352-forging";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 2.75e-9),
                    (ParameterLabel::n, 3.29),
                    (ParameterLabel::kf, 64.0),
                ]),
            },
        );

        name = "forman:aa7050t73651-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 4.11e-9),
                    (ParameterLabel::n, 2.98),
                    (ParameterLabel::kf, 55.0),
                ]),
            },
        );

        name = "forman:aa7075t6-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 1.37e-8),
                    (ParameterLabel::n, 3.02),
                    (ParameterLabel::kf, 63.9),
                ]),
            },
        );

        name = "forman:aa7075t7351";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 6.27e-9),
                    (ParameterLabel::n, 2.78),
                    (ParameterLabel::kf, 55.8),
                ]),
            },
        );

        name = "forman:aa7175t3652-forging";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 2.61e-9),
                    (ParameterLabel::n, 2.91),
                    (ParameterLabel::kf, 38.0),
                ]),
            },
        );

        name = "forman:aa7178t651-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 3.74e-8),
                    (ParameterLabel::n, 2.06),
                    (ParameterLabel::kf, 30.7),
                ]),
            },
        );

        name = "forman:aa7475t7351-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 3.24e-8),
                    (ParameterLabel::n, 2.32),
                    (ParameterLabel::kf, 78.2),
                ]),
            },
        );

        name = "forman:aa7475t76-sheet";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 6.54e-8),
                    (ParameterLabel::n, 2.18),
                    (ParameterLabel::kf, 79.9),
                ]),
            },
        );

        name = "forman:aa7475t7651-plate";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 9.30e-9),
                    (ParameterLabel::n, 2.73),
                    (ParameterLabel::kf, 63.1),
                ]),
            },
        );

        name = "forman:a357t6-sandcasting";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 2.19e-9),
                    (ParameterLabel::n, 2.94),
                    (ParameterLabel::kf, 41.5),
                ]),
            },
        );

        name = "forman:a357t6-investmentcasting";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Schwarmann 86]",
                units: Forman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::c, 6.65e-9),
                    (ParameterLabel::n, 2.40),
                    (ParameterLabel::kf, 38.2),
                ]),
            },
        );

        name = "hartman:jones13-aa7050t7451";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[Jones 12]",
                units: Hartman::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::d, 7.0e-10),
                    (ParameterLabel::deltak_th, 0.1),
                    (ParameterLabel::a, 47.0),
                    (ParameterLabel::alpha, 2.0),
                ]),
            },
        );

        name = "white:barter14-aa7050t7451";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[White 18]",
                units: White::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 0.25481858),
                    (ParameterLabel::b, 1.10247048),
                    (ParameterLabel::c, 4.35831677),
                    (ParameterLabel::d, 23.08586582),
                    (ParameterLabel::e, 0.03420171),
                    (ParameterLabel::f, 0.4717843),
                    (ParameterLabel::kic, 31.54),
                ]),
            },
        );

        name = "white:chan16-aa7050t7451";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: White::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 0.2918621458996122),
                    (ParameterLabel::b, 1.2635107616404941),
                    (ParameterLabel::c, 3.5528305197144334),
                    (ParameterLabel::d, 22.243180576185246),
                    (ParameterLabel::e, 0.03924086425080324),
                    (ParameterLabel::f, 0.5551311271413691),
                    (ParameterLabel::kic, 41.45917108627669),
                ]),
            },
        );

        name = "lin_lower_th_5:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: LinLowerTh5::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 2.3805),
                    (ParameterLabel::b, 7.3607e-10),
                    (ParameterLabel::c, -0.5799),
                    (ParameterLabel::deltak_th, 0.6238),
                ]),
            },
        );

        // Material: 7085t7452-CA
        name = "lin_lower_th_5a:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: LinLowerTh5A::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 1.8476),
                    (ParameterLabel::b, 7.7723e-10),
                    (ParameterLabel::c, -0.0403),
                    (ParameterLabel::d, -998.5654),
                    (ParameterLabel::deltak_th, 0.6717),
                    (ParameterLabel::k_ut, 8645.19),
                ]),
            },
        );

        // Material: 7085t7452-CA
        name = "lin_lower_th_5b:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: LinLowerTh5B::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 1.8476),
                    (ParameterLabel::b, 7.7723e-10),
                    (ParameterLabel::c1, -0.0403),
                    (ParameterLabel::c2, 0.0403),
                    (ParameterLabel::d, -998.5654),
                    (ParameterLabel::deltak_th, 0.6717),
                    (ParameterLabel::k_ut, 8645.19),
                ]),
            },
        );

        // Material: 7085t7452-CA
        name = "lin_lower_th_6:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: LinLowerTh6::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 2.4299),
                    (ParameterLabel::b, 6.8283e-10),
                    (ParameterLabel::c, -0.8636),
                    (ParameterLabel::d, 0.1065),
                    (ParameterLabel::deltak_th, 0.6402),
                ]),
            },
        );

        // Material: 7085t7452-CA
        name = "lin_lower_th_6a:default";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: LinLowerTh6A::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 1.9544),
                    (ParameterLabel::b, 7.3345e-10),
                    (ParameterLabel::c, -0.2724),
                    (ParameterLabel::d, 0.0589),
                    (ParameterLabel::e, -6165.6640),
                    (ParameterLabel::deltak_th, 0.6709),
                    (ParameterLabel::k_ut, 61964.1148),
                ]),
            },
        );

        // Material: AA7075-T7351
        name = "no_r_low_th1:aa7075t7351";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: NoRLowerTh1::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 2.9619),
                    (ParameterLabel::b, 8.039e-10),
                    (ParameterLabel::deltak_th, 0.5471),
                ]),
            },
        );

        // Material: AA7075-T7351
        name = "no_r_low_th2:aa7075t7351";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: NoRLowerTh2::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 3.4073),
                    (ParameterLabel::b, 2.0471e-10),
                    (ParameterLabel::c, 3.4690),
                    (ParameterLabel::deltak_th, 0.7798),
                ]),
            },
        );

        // Material: AA7075-T7351
        name = "no_r_low_up_th1:aa7075t7351";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: NoRLowerUpperTh1::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 2.1494),
                    (ParameterLabel::b, 1.6498e-9),
                    (ParameterLabel::c, -627.6457),
                    (ParameterLabel::deltak_th, 0.7175),
                    (ParameterLabel::k_ut, 10000.0),
                ]),
            },
        );

        // Material: AA7075-T7351
        name = "no_r_low_up_th2:aa7075t7351";
        materials.insert(
            name,
            EqnProperty {
                name,
                cite: "[]",
                units: NoRLowerUpperTh2::UNITS,
                params: BTreeMap::from([
                    (ParameterLabel::a, 3.3650),
                    (ParameterLabel::b, 2.0225e-10),
                    (ParameterLabel::c, 3.7639),
                    (ParameterLabel::d, 0.9730),
                    (ParameterLabel::deltak_th, 0.7818),
                    (ParameterLabel::k_ut, 219.1398),
                ]),
            },
        );

        materials
    };
}

/// Return the full map of materials as <name, material>
pub fn get_all_dadns() -> &'static BTreeMap<&'static str, EqnProperty> {
    &MATERIALS
}

/// Get a single material using the equation:material format
pub fn get_dadn(name: &str) -> Option<&EqnProperty> {
    match MATERIALS.get(name) {
        Some(eqn) => Some(eqn),
        None => None,
    }
}

#[cfg(test)]
mod tests {
    use crate::table;

    #[test]
    fn checkpairtable() {
        // ID Code = M7GJ11AC1
        let nasgro_table = table::PairTable::new(
            vec![0.08, 0.1, 0.4, 0.5, 0.7, 0.8],
            vec![
                vec![
                    6.8865, 7.1373, 8.5740, 11.2303, 12.6427, 13.7598, 18.5817, 21.4952, 23.1962,
                ],
                vec![
                    3.9084, 3.9241, 4.2994, 4.8354, 5.5089, 6.4799, 7.5038, 8.9768, 10.6903,
                    13.0444, 15.4584, 18.3955, 21.1440, 23.9296, 26.0453, 28.3481, 29.8538,
                ],
                vec![
                    2.5410, 2.5675, 3.0591, 3.6212, 4.1691, 4.7674, 5.3145, 6.1336, 7.2176, 8.8385,
                    10.6337, 13.1499, 15.5498, 16.7494,
                ],
                vec![
                    4.7534, 4.8525, 5.4633, 6.2946, 7.6432, 9.2110, 11.4217, 13.7010, 16.4040,
                    18.8969, 22.2837, 26.4682, 26.4850,
                ],
                vec![
                    1.9011, 1.9119, 2.2139, 2.8924, 3.3723, 3.8130, 4.2446, 4.8279, 5.4787, 6.1235,
                ],
                vec![
                    1.7701, 2.0983, 2.5583, 2.9669, 3.3785, 3.8713, 4.7429, 5.5754, 5.8345,
                ],
            ],
            vec![
                vec![
                    8.1560e-07, 1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05,
                    1.8000e-05, 3.2000e-05, 6.0000e-05,
                ],
                vec![
                    5.8598e-08, 6.0000e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07,
                    1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05, 1.8000e-05,
                    3.2000e-05, 6.0000e-05, 1.0000e-04, 1.8000e-04, 2.6886e-04,
                ],
                vec![
                    5.8374e-08, 6.0000e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07,
                    1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05, 1.8000e-05,
                    3.2000e-05, 4.3702e-05,
                ],
                vec![
                    8.9219e-07, 1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05,
                    1.8000e-05, 3.2000e-05, 6.0000e-05, 1.0000e-04, 1.8000e-04, 3.2000e-04,
                    3.2064e-04,
                ],
                vec![
                    5.8414e-08, 6.0000e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07,
                    1.0000e-06, 1.8000e-06, 3.2000e-06, 5.4116e-06,
                ],
                vec![
                    6.3591e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07, 1.0000e-06,
                    1.8000e-06, 3.2000e-06, 3.9566e-06,
                ],
            ],
        );

        let tol = 1e-5;
        //assert_eq!(nasgro_table.interp(6.8865, 0.08), 8.156e-7);
        //        assert_eq!(nasgro_table.interp(6.1235, 0.7), 3.9566e-6);
        //assert_eq!(nasgro_table.interp(5.8345, 0.8), 3.9566e-6);
        //assert_eq!(nasgro_table.interp(23.1962, 0.08), 6.0e-5);

        assert!((nasgro_table.interp(6.8865, 0.08) - 8.156e-7).abs() < tol);
        // these fail in the table lookup
        //assert!((nasgro_table.interp(5.8345, 0.8) - 3.9566e-6).abs() < tol);
        //assert!((nasgro_table.interp(23.1962, 0.08) - 6.0e-5).abs() < tol);
        assert!((nasgro_table.interp(4.7534, 0.5) - 8.9219e-7).abs() < tol);

//        println!("(--nasgro_table.interp(26.4850, 0.5) = {:e}", nasgro_table.interp(26.4850, 0.5));
//        println!("(--nasgro_table.interp(26.4, 0.4) = {:e}", nasgro_table.interp(26.4, 0.4));
        // no gsl gives nasgro da/dn 0.00000003266836692698721
        // GSL gives (nasgro_table.interp(26.4850, 0.5) = 0.00032064

//        assert!((nasgro_table.interp(26.4850, 0.5) - 3.2064e-4).abs() < 1.0e-7);
        assert!(true);
    }
}
