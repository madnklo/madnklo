
extern crate madnklo;
extern crate rand;

use madnklo::phase_space_generator::FlatPhaseSpaceGenerator;
use madnklo::vector::LorentzVector;
use rand::Rng;

fn main() {
    let mut res = vec![LorentzVector::default(); 8];
    let mut rng = rand::thread_rng();

    let fpsg = FlatPhaseSpaceGenerator::new((0..8).map(|i| 100. + 10. * i as f64).collect::<Vec<_>>());
    for x in 0..1 {
        //let rands: Vec<f64> = (0..res.len() * 3).map(|_| rng.gen()).collect();
        let rands = vec![
            0.9346490925737034,
            0.7640013788037885,
            0.6036304724712833,
            0.952669596156114,
            0.910539279844903,
            0.8384921489436372,
            0.14241075965631556,
            0.9144834752658628,
            0.10476932994117161,
            0.6341194437942385,
            0.3765031864448587,
            0.5179499127936523,
            0.09265998612840654,
            0.4878250698298594,
            0.5319210066391838,
            0.6550897942162808,
            0.3655625679728105,
            0.23444679221853548,
            0.09297003102257273,
            0.16332138769108173,
        ];

        let v = fpsg.generate(
            5000.,
            &rands,
            &mut res,
        );
        println!("weight={:e}", v);

        for (i, p) in res.iter().enumerate() {
            println!("p{}={:e}", i + 1, p);
        }
    }
}
