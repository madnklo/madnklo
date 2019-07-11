
use std::f64::consts::{FRAC_PI_2, PI};
use vector::LorentzVector;

pub trait PhaseSpaceGenerator {
    fn generate(&mut self, e_cm: f64, x: &[f64], ps: &mut [LorentzVector<f64>]) -> f64;
}

pub struct FlatPhaseSpaceGenerator {
    volume_factors: Vec<f64>,
    masses: Vec<f64>,
    r: Vec<f64>, // rescaled input
}

impl FlatPhaseSpaceGenerator {
    const MAX_EXTERNAL: usize = 20;
    const PRECISION: f64 = 1e-15;
    const EPSILON_BORDER: f64 = 1e-10;
    const MAXIMUM_DERIVATIVE: [f64; FlatPhaseSpaceGenerator::MAX_EXTERNAL] = [
        0., 0., 0.5, 0.7, 0.75, 0.8, 0.805, 0.81, 0.82, 0.83, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9,
        0.9, 0.9, 0.9, 0.9,
    ];

    pub fn new(masses: Vec<f64>) -> FlatPhaseSpaceGenerator {
        let mut volume_factors = vec![0., 0.];

        for n in 2..FlatPhaseSpaceGenerator::MAX_EXTERNAL {
            let mut f = 1.;
            for i in 2..=n - 2 {
                f *= i as f64;
            }

            volume_factors.push(FRAC_PI_2.powi(n as i32 - 1) / f.powi(2) / (n - 1) as f64);
        }

        let r = vec![0f64; masses.len() * 3 - 4];

        FlatPhaseSpaceGenerator {
            volume_factors,
            masses,
            r
        }
    }

    #[inline]
    fn rho(mp: f64, m: f64, mps: f64) -> f64 {
        let mp_sq = mp * mp;
        ((mp_sq - (m + mps) * (m + mps)) * (mp_sq - (m - mps) * (m - mps))).sqrt() / mp_sq * 0.125
    }

    /// Solve `p*x^(p-1)-(p-1)*x^p-r=0`
    fn get_u(p: usize, r: f64) -> f64 {
        if r == 0. || r == 1.0 {
            return r;
        }

        debug_assert!(p > 1);
        if p == 2 {
            let r = 1. - (1. - r).sqrt();
            return r;
        }

        // numerically find root by Newton's method
        let mut x = FlatPhaseSpaceGenerator::MAXIMUM_DERIVATIVE[p];
        loop {
            let xp = x.powi(p as i32 - 2);
            let eval = -r + xp * x * (p as f64 - p as f64 * x + x);
            if eval.abs() < Self::PRECISION {
                break;
            }

            let dx = (p as f64 - 1.) * (p as f64 - p as f64 * x) * xp;
            x -= eval / dx;
        }

        debug_assert!(x >= 0. && x <= 1.0);
        x
    }
}

impl PhaseSpaceGenerator for FlatPhaseSpaceGenerator {
    fn generate(&mut self, e_cm: f64, x: &[f64], ps: &mut [LorentzVector<f64>]) -> f64 {
        let mut q = LorentzVector::from_args(e_cm, 0., 0., 0.);
        let mut mass_sum = self.masses.iter().sum::<f64>();
        let mut m = q.square().sqrt() - mass_sum;
        let n = ps.len();
        debug_assert!(n < FlatPhaseSpaceGenerator::MAX_EXTERNAL);
        let mut weight = self.volume_factors[n] * m.powi(2 * n as i32 - 4);

        for (rr, xr) in self.r.iter_mut().zip(x) {
            *rr = xr.max(FlatPhaseSpaceGenerator::EPSILON_BORDER).min(1. - FlatPhaseSpaceGenerator::EPSILON_BORDER);
        }
        
        for i in 0..n - 1 {
            let mi = m + mass_sum; // compute the intermediate mass

            // note: we take the sqrt of u, since Mi^2=u2*..*u_i * M^2
            // in the paper this relation is incorrect
            // additionally, u is fixed to 0 for i=n-2
            let u = if i == n - 2 {
                0.
            } else {
                // TODO: change index convention to i*3
                FlatPhaseSpaceGenerator::get_u(n - i - 1, self.r[i]).sqrt()
            };

            let rho =
                FlatPhaseSpaceGenerator::rho(mi, m * u + mass_sum - self.masses[i], self.masses[i]);
            let qi = 4. * mi * rho;

            if i == n - 2 {
                weight *= 8. * FlatPhaseSpaceGenerator::rho(mi, self.masses[i + 1], self.masses[i]);
            } else {
                weight *= rho / FlatPhaseSpaceGenerator::rho(m, m * u, 0.)
                    * (m * u + mass_sum - self.masses[i])
                    / (m * u);
            }

            m *= u;

            let cos_theta = 2. * self.r[n - 2 + 2 * i] - 1.;
            let sin_theta = (1. - cos_theta * cos_theta).sqrt();
            let phi = 2. * PI * self.r[n - 1 + 2 * i];
            let (sin_phi, cos_phi) = phi.sin_cos();

            ps[i] = LorentzVector::from_args(
                qi.hypot(self.masses[i]),
                qi * cos_phi * sin_theta,
                qi * sin_phi * sin_theta,
                qi * cos_theta,
            );

            ps[i] = ps[i].boost(&(q / q.t));
            q = q - ps[i];

            mass_sum -= self.masses[i];
        }
        ps[n - 1] = q;

        weight
    }
}
