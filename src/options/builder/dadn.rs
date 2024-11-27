use crate::options::EasiOptions;
use fatigue::{
    dadn::{self, relabel_parameters, DaDn},
    material, table, COMMENT,
};
use log::error;

pub fn get_dadn(options: &EasiOptions, params: &[f64]) -> Box<dyn DaDn + Send + Sync> {
    // Create a map of the new values to use in instantiating the dadn model
    // This relies on the params provided to this method to be sorted in the
    // same order as the options map
    let remapped_params = match relabel_parameters(params, &options.dadn) {
        Ok(result) => result,
        Err(why) => {
            error!("Error: {}", why);
            std::process::exit(1)
        }
    };

    let deltak_th = if remapped_params.contains_key(&dadn::ParameterLabel::deltak_th) {
        *remapped_params.get(&dadn::ParameterLabel::deltak_th).unwrap()
    } else {
        options.deltak_th.unwrap_or(0.0)
    };

    let kmax_th = if options.kmax_th_zero {
        0.0
    } else {
        deltak_th
    };

    let dadn_options = dadn::Options {
        rmax: options.rmax,
        rmin: options.rmin,
        deltak_th,
        kneg: options.kneg,
        kmax_th,
        kmax_th_follows_deltak_th: !options.kmax_th_zero,
    };

    match dadn::make_model(&options.dadn, &remapped_params, dadn_options) {
        Ok(equation) => equation,
        Err(why) => {
            error!("Error: {}", why);
            std::process::exit(1)
        }
    }
}

pub fn get_dadn_params(options: &EasiOptions) -> Result<(Vec<f64>, String), String> {
    let mut params = options.params.values().cloned().collect::<Vec<f64>>();
    let mut message = format!("{}Using provided dadn parameters", COMMENT);
    if params.is_empty() {
        // extract out the appropriate material parameters from a file
        params = if options.dadn.starts_with("file:") {
            let filename = options.dadn.trim_start_matches("file:");
            message = format!(
                "{}No parameters given, using the dk values in the dadn file {}",
                COMMENT, filename
            );
            let table = table::Table::read_file(filename, true);
            // collapse down the dks and use these as the parameters for optimising
            table.variables()
        // or from the internal database.
        } else {
            match material::get_dadn(&options.dadn) {
                Some(result) => {
                    message = format!(
                        "{}No parameters given, obtaining from material library for {}",
                        COMMENT, options.dadn
                    );
                    result.params.values().cloned().collect()
                }
                None => {
                    message = format!("Unknown dadn '{}', or no dadn material or parameters have been provided", options.dadn);
                    return Err(message);
                }
            }
        }
    }

    Ok((params, message))
}
