use fatigue::{beta, material};

static HIGHLIGHTS: &str = 
"    - Load sequence modification including truncation and clipping
    - Fatigue cycle counting and cycle modification including truncation
    - Linear-Elastic Fracture Mechanics-based fatigue crack growth prediction
    - Built-in beta factor solutions and da/dn vs. dK models
    - Optimisation of fatigue crack growth model parameters for improved
      predictions under realistic spectrum loading conditions
    - Prediction of the fracture surface appearance caused by fatigue under
      spectrum loading
";

static UNITS: &str =
"The units employed in Easigrow should be consistent with those of the da/dN
versus dK data. The default model parameters provided by Easigrow relate dK in
MPa sqrt(m) to da/dN in m/cycle. Therefore, when employing these default
parameters the scaled load spectrum input should be in MPa, while crack lengths
and geometry should be in meters (m). Easigrow can also be used with imperial
units, where ksi and inches (in) can be used with a da/dN versus dK relationship
that predicts growth in in/cycle from dK in ksi sqrt(in).

In general, Easigrow calculates K from a stress input. However, for the
compact-tension specimen geometry, the beta factor should be used in concert
with load, not stress. In this case, the scaled spectrum should be in
meganewtons (MN) when employing metric units and kilopounds (kip) when employing
imperial units. 
";

/// Prints out lists of data. Sort of an extended help.
pub fn print_list() {
    // List of the valid names for producing output
    let output = [
        ("block", "block number"),
        ("loadline_no", "load line number"),
        ("cycle_no", "cycle number"),
        ("a/c", "crack aspect ratio at the end of the current cycle"),
        ("a/d", "cracked fraction of forward distance at the end of the current cycle"),
        ("c/b", "cracked fraction of sideways distance at the end of the current cycle"),
        ("kmax_a", "kmax for the current cycle calculated in the a-direction"),
        ("kmax_c", "kmax for the current cycle calculated in the c-direction"),
        ("kmin_a", "kmin for the current cycle calculated in the a-direction"),
        ("kmin_c", "kmin for the current cycle calculated in the c-direction"),
        ("dk_a", "dK for the current cycle calculated in the a-direction"),
        ("dk_c", "dK for the current cycle calculated in the c-direction"),
        ("R", "current cycle's R"),
        ("beta_a", "beta at a calculated before FCG for the current cycle"),
        ("beta_c", "beta at c calculated before FCG for the current cycle"),
        ("a", "crack length in the a-direction at the end of the current cycle"),
        ("c", "crack length in the c-direction at the end of the current cycle"),
        ("a_prior", "crack length in the a-direction at the start of the current cycle"),
        ("c_prior", "crack length in the c-direction at the start of the current cycle"),
        ("da", "fatigue crack growth increment in the a-direction for the current cycle"),
        ("dc", "fatigue crack growth increment in the c-direction for the current cycle"),
        ("peak", "scaled peak (maximum) stress for the current cycle"),
        ("valley", "scaled valley (minimum) stress for the current cycle"),
    ];

    let formats = [
        (
            "Crack growth measurement files",
"These contain measurements typically collected from fatigue coupon tests. They
are used to estimate measured fatigue crack growth rates, which are the
benchmark that a model is optimised to match.

The crack growth measurement file is in the following format:


<block> <a>
...

or

<line> <block> <a>
...

or 

<a>
...
                   
The same format must be used for an entire file. Here <line> represents the
corresponding line no. (starting at 0) of the load sequence, and <block> is the
repeat of the load sequence (also starting at zero) when the measured crack
length was reached. Easigrow uses the difference between the measurements and
times on successive lines to estimate crack growth rates. It is recommended that
successive measurements are separated by whole blocks. For the third format
shown above, Easigrow assumes that 1 block separates the measurements on
successive lines. It is not uncommon for crack length measurements to be made
along different paths when collecting measurements via fractography. For an
Easigrow crack growth measurement file, blank lines are used separate
measurements along different paths. Easigrow only calculates crack growth rates
using measurements from the same path.
",
        ),
        (
            "Optimisation files",
"Each line in an optimisation file specifies a set of measurements that are used
to estimate fatigue crack growth rates, which Easigrow will optimise a model to
match. The remainder of the line contains commands that instruct Easigrow on how
to predict these crack growth rates. These instructions include specifying the
load spectrum, geometry, and beta factor solution appropriate for the fatigue
coupon test responsible for generating the measurements. The ‘easigrow’ command
is omitted from these commands.

The format of the optimisation file is:

--crack_infile <Measurement file 1> <Easigrow commands to predict measurement set 1>
--crack_infile <Measurement file 2> <Easigrow commands to predict measurement set 2>
...
",
        ),
        (
            "Tabular beta factor solution files",
            "Two options allow the user to provide beta factor solutions via tabular data:
--beta file-1d:<FILE>
--beta file-2d:<FILE>

1-D beta factor solution files
------------------------------

The beta factor solution is provided as a function of the ratio between the
crack length and the distance between the crack origin and the free-edge the
crack is growing towards (c/sideways). The format required is as follows:

# Comment describing the contents of the file
c_on_sideways1 beta1
c_on_sideways2 beta2
c_on_sideways3 beta3
...

2-D beta factor solution files
------------------------------

The 2-dimensional tabular option provides a beta factor for the a-direction only,
and is only valid before the crack grows through the thickness (i.e., a < forward).
The required format consists of multiple columns, separated by whitespace,
specifying the beta factor as a function of the ratio a/forward for different
a/c values. The format required is as follows:

a_on_c1 a_on_c2 a_on_c3 ...
a_on_forward1 beta_11 beta_12 beta_13 ...
a_on_forward2 beta_21 beta_22 beta_23 ...
a_on_forward3 beta_23 beta_32 beta_33 ...
...

It is important to note that because the beta factor is only calculated in the
a-direction, fatigue crack growth can only be predicted in this direction. Since
crack growth cannot be predicted in the c-direction, it is not possible to
accurately update the a/c ratio. Therefore, the 2-dimensional tabular beta
factor solution employs a fixed a/c value throughout the FCG prediction, based
on the provided --astart and --cstart values.
",
        ),
        (
            "Tabular da/dN versus dK data file",
"The first row defines the R levels at which dK values are supplied, while the
first column lists the da/dN values for which dK is defined. Each subsequent
column then contains the dK values corresponding to the listed da/dN values for
each R ratio. The format required is as follows:

# Comment describing the contents of the file
R1 R2 R3 ...
da/dN_1 dK1_R1 dK1_R2 dK1_R3 ...
da/dN_2 dK2_R1 dK2_R2 dK2_R3 ...
da/dN_3 dK3_R1 dK3_R2 dK3_R3 ...
...

For example:

-1 -0.5 0 0.5
1.00E-08 4.004 4.215 4.531 3.692
1.00E-06 25.49 26.40 27.82 19.23
1.00E-04 80.48 81.94 84.44 50.82 
",
        ),
    ];

    let biblio = [
["[Broek 86]", "Broek, D. (1986) Elementary Engineering Fracture Mechanics. 4th ed.
           The Hague, Martinus Nijhoff"
        ],

["[Newman 86]", "Newman, J. C. and Raju, I. S. (1986) Stress-Intensity Factor
            Equations for Cracks in Three-Dimensional Finite Bodies Subjected to
            Tension and Bending Loads. In: Computational Methods in the
            Mechanics of Fracture. Elsevier Science Publishers 312-334"],

["[Anderson 05]", "Anderson, T. L. (2005) Fracture Mechanics: Fundamentals and
              Applications. 3rd ed, CRC Press"],
        
["[Tada 73]", "Tada, H., Paris, P. C. and Irwin, G. R. (1973) The Stress Analysis of
          Cracks Handbook, Del. Research Corporation"],

["[ASTM 97]", "Standard Test Method for Plane-Strain Fracture Toughness of Metallic
          Materials. (1997) ASTM E-399-90, ASTM International"],

["[Fedderson 66]", "Fedderson, C. E. (1966) Discussion to: \"Plane Strain Crack
               Toughness Testing of High Strength Metallic Materials\" by W. F.
               Brown and J. E. Srawley. ASTM STP 410, Philadelphia, ASTM
               International"],

["[Bowie 56]", "Bowie, O. L. (1956) Analysis of a Infinite Plate Containing Radial
           Cracks Originating at the Boundary of an Internal Circular Hole.
           Journal of Mathematics and Physics 35 60-71"],
        
["[Murakami 87]", "Stress Intensity Factors Handbook (1987) Murakami, Y. ed., Pergamon"],

["[Murakami 86]", "Murakami, Y. and Tsuru, T. (1987) Stress Intesity Factor Equations
              for a Semi-Elliptical Surface Crack in a Shaft under Bending. In: 
              Stress Intensity Factors Handbook. Soc. Mater. Sci., Japan No. 
              87-0164 B (p2062-2065)"],

["[Forman 86]", "Forman, R. G. and Shivakumar, V. (1986) Growth Behavior of Surface
            Cracks in the Circumferential Plane of Solid and Hollow Cylinders.
            In: Underwood, J. H., et al. (eds.) Fracture Mechanics: Seventeenth
            Volume. ASTM STP 1131, Philadelphia 519-546"],

["[Schwarmann 86]", "Schwarmann, L. (1986) Material data of high strength aluminium
                alloys for durability evaluation of structures: fatigue strength,
                crack propagation, fracture toughness. Dusseldorf, Germany,
                Aluminium-Verlag"],

["[Jones 12]", "Jones, R., Molent, L. and Walker, K. (2012) Fatigue crack growth in a
           diverse range of materials. International Journal of Fatigue 40 43-50"],
        
["[Forman 05]", "Forman, R. G., et al. (2005) Fatigue Crack Database For Damage
            Tolerance Analysis. DOT/FAA/AR-05/15, Office of Aviation Research,
            Washington, D. C., Federal Aviation Authority"],

["[White 18]", "White, P., Mongru, D. and Molent, L. (2018) A Crack Growth Based
           Individual Aircraft Monitoring Method Utilizing a Damage Metric.
           Structural Health Monitoring 17 (5) 1178-1191"],

    ];

    // Set up a new counter to automatically label the section headers
    let mut header = Counter::new();

    header.section("Program Highlights");
    print!("{}", HIGHLIGHTS);

    header.section("Units");
    print!("{}", UNITS);

    header.section("Output Parameters");

    for &(abbrev, descrip) in &output {
        println!("{:20} {}", abbrev, descrip);
    }

    header.section("Beta Models");
    // List of the beta equations that are implemented
    let betas = beta::get_all_betas();
    println!("{:20} {:4} {:30} Description", "Name", "DOI*", "Parameters**");
    for beta in &betas {
        println!(
            "{:20} {:4} {:30} {} {}",
            beta.name,
            beta.direction_of_interest.to_lowercase(),
            beta.args.to_string(),
            beta.summary,
            beta.cite
        );
    }
    println!("\n*DOI: Direction of interest");
    println!("\n**Key to user inputs:");
    println!("f = forward; s = sideways; aa = a_angle; ca = c_angle; as = astart; cs = cstart; r = radius");

    header.section("Cycle Counting  Methods");

    println!(
"A load sequence needs to be transformed into a series of distinct fatigue cycles
so that cycle-by-cycle fatigue crack growth can be predicted. These cycles can
either be directly input from a file or derived from a turning point load
sequence via cycle counting. Easigrow offers two cycle counting methods:
rainflow (the default) and tension.

rainflow  Smaller interruption cycles within larger cycles are successively
          extracted until all peak and valley loads have been paired into cycles.

tension   Valley loads are paired with the peak load that follows.

After cycles counting, Easigrow orders the extracted cycles according to where
the peaks in each occurred within the original load sequence. ");

    header.section("da/dN versus dK Models");
    println!(
"Built-in da/dN versus dK models are specified according to EQUATION:material.

EQUATION = the name of the da/dN equation.
material = specifies equation parameters for the chosen material.

If parameter values are specified via the --parameters option, these will be
used instead of the standard library values.\n");

    println!("{:35} {:20} Coefficients", "Name", "Ref.");
    let materials = material::get_all_dadns();
    for (name, mat) in materials.iter() {
        print!("{:35} {:20} ", name, mat.cite);
        for (label, value) in &mat.params {
            if (*value).abs() < 1.0 {
                print!("{}={:e} ", label, value);
            } else {
                print!("{}={} ", label, value);
            }
            
        }
        println!();
    }
    println!(
        "{:35} {:20} {:30}",
        "file:FILE", " ", "Read FILE of tabular dadn data."
    );

    header.section("File formats");

    for &(file, form) in &formats {
        header.subsection(file);
        println!("{}", form);
    }

    header.section("References");

    for bib in &biblio {
        println!("{} {}\n", bib[0], bib[1]);
    }

    println!();
}

struct Counter {
    section: usize,
    subsection: usize,
}

impl Counter {
    fn new() -> Counter {
        Counter {
            section: 0,
            subsection: 0,
        }
    }

    // print as a header
    fn section(&mut self, head: &str) {
        self.section += 1;
        let header = format!("{}. {}", self.section, head);
        println!("\n{}", header);
        // Underline
        for _ in 0..header.len() {
            print!("=");
        }
        println!("\n");
    }

    fn subsection(&mut self, head: &str) {
        self.subsection += 1;
        let header = format!("{}.{}. {}", self.section, self.subsection, head);
        println!("{}", header);
        // Underline
        for _ in 0..header.len() {
            print!("-");
        }
        println!("\n");
    }
}
