// NonNan numbers Floating point numbers cannot be sorted because they
// can be Nan.  We introduce a new type of number guaranteed not to be
// an NaN to allow floats to be sorted or the maximum found.

use std::cmp::Ordering;

#[derive(PartialEq, Clone, Debug)]
pub struct NonNan(f64);

impl NonNan {
    pub fn new(val: f64) -> Option<NonNan> {
        if val.is_nan() {
            None
        } else {
            Some(NonNan(val))
        }
    }
}

impl Eq for NonNan {}
impl Ord for NonNan {
    fn cmp(&self, other: &NonNan) -> Ordering {
        self.0.partial_cmp(&other.0).unwrap()
    }
}

impl PartialOrd for NonNan {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
