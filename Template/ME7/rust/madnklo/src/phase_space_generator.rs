use std::f64::consts::{FRAC_PI_2, PI};
use vector::LorentzVector;

pub trait PhaseSpaceGenerator {
    fn get_PS_point(
        &mut self,
        x: &[f64],
        ps: &mut [LorentzVector<f64>],
    ) -> (f64, (f64, Option<f64>), (f64, Option<f64>));

    fn get_dimensions(&self) -> usize;
}

pub struct FlatPhaseSpaceGenerator {
    volume_factors: Vec<f64>,
    pub masses: (Vec<f64>, Vec<f64>),
    r: Vec<f64>, // rescaled input
    collider_energy: f64,
    beam_type: (isize, isize),
    correlated_beam_convolution: bool,
    is_beam_factorization_active: (bool, bool),
}

impl FlatPhaseSpaceGenerator {
    const MAX_EXTERNAL: usize = 20;
    const PRECISION: f64 = 1e-15;
    const EPSILON_BORDER: f64 = 1e-10;
    const MAXIMUM_DERIVATIVE: [f64; FlatPhaseSpaceGenerator::MAX_EXTERNAL] = [
        0., 0., 0.5, 0.7, 0.75, 0.8, 0.805, 0.81, 0.82, 0.83, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9,
        0.9, 0.9, 0.9, 0.9,
    ];

    // The lowest value that the center of mass energy can take.
    // We take here 1 GeV, as anyway below this non-perturbative effects dominate
    // and factorization does not make sense anymore
    const ABSOLUTE_E_CM_MIN: f64 = 1.;

    pub fn new(
        n_initial: usize,
        masses: (Vec<f64>, Vec<f64>),
        collider_energy: f64,
        beam_type: (isize, isize),
        correlated_beam_convolution: bool,
        is_beam_factorization_active: (bool, bool),
    ) -> FlatPhaseSpaceGenerator {
        if n_initial != 2 {
            unimplemented!("Only n_initial = 2 is supported at the moment");
        }

        let mut volume_factors = vec![0., 0.];

        for n in 2..FlatPhaseSpaceGenerator::MAX_EXTERNAL {
            let mut f = 1.;
            for i in 2..=n - 2 {
                f *= i as f64;
            }

            volume_factors.push(FRAC_PI_2.powi(n as i32 - 1) / f.powi(2) / (n - 1) as f64);
        }

        let r = vec![0f64; masses.1.len() * 3 - 4];

        FlatPhaseSpaceGenerator {
            volume_factors,
            masses,
            r,
            collider_energy,
            beam_type,
            correlated_beam_convolution,
            is_beam_factorization_active,
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

    pub fn generate(&self, e_cm: f64, x: &[f64], ps: &mut [LorentzVector<f64>]) -> f64 {
        let mut q = LorentzVector::from_args(e_cm, 0., 0., 0.);
        let mut mass_sum = self.masses.1.iter().sum::<f64>();
        let mut m = q.square().sqrt() - mass_sum;
        let n = ps.len();
        debug_assert!(n < FlatPhaseSpaceGenerator::MAX_EXTERNAL);
        let mut weight = self.volume_factors[n] * m.powi(2 * n as i32 - 4);

        for i in 0..n - 1 {
            let mi = m + mass_sum; // compute the intermediate mass

            // note: we take the sqrt of u, since Mi^2=u2*..*u_i * M^2
            // in the paper this relation is incorrect
            // additionally, u is fixed to 0 for i=n-2
            let u = if i == n - 2 {
                0.
            } else {
                // TODO: change index convention to i*3
                FlatPhaseSpaceGenerator::get_u(n - i - 1, x[i]).sqrt()
            };

            let rho = FlatPhaseSpaceGenerator::rho(
                mi,
                m * u + mass_sum - self.masses.1[i],
                self.masses.1[i],
            );
            let qi = 4. * mi * rho;

            if i == n - 2 {
                weight *=
                    8. * FlatPhaseSpaceGenerator::rho(mi, self.masses.1[i + 1], self.masses.1[i]);
            } else {
                weight *= rho / FlatPhaseSpaceGenerator::rho(m, m * u, 0.)
                    * (m * u + mass_sum - self.masses.1[i])
                    / (m * u);
            }

            m *= u;

            let cos_theta = 2. * x[n - 2 + 2 * i] - 1.;
            let sin_theta = (1. - cos_theta * cos_theta).sqrt();
            let phi = 2. * PI * x[n - 1 + 2 * i];
            let (sin_phi, cos_phi) = phi.sin_cos();

            ps[i] = LorentzVector::from_args(
                qi.hypot(self.masses.1[i]),
                qi * cos_phi * sin_theta,
                qi * sin_phi * sin_theta,
                qi * cos_theta,
            );

            ps[i] = ps[i].boost(&(q / q.t));
            q = q - ps[i];

            mass_sum -= self.masses.1[i];
        }
        ps[n - 1] = q;

        weight
    }
}

impl PhaseSpaceGenerator for FlatPhaseSpaceGenerator {
    #[inline]
    fn get_dimensions(&self) -> usize {
        self.masses.1.len() * 3 - 4
            + if self.correlated_beam_convolution {
                1
            } else {
                0
            }
            + if self.is_beam_factorization_active.0 {
                1
            } else {
                0
            }
            + if self.is_beam_factorization_active.1 {
                1
            } else {
                0
            }
            + if self.beam_type == (0, 0) {
                0
            } else {
                if self.masses.0.len() == 2 && self.masses.1.len() == 1 {
                    1
                } else {
                    2
                }
            }
    }

    fn get_PS_point(
        &mut self,
        x: &[f64],
        ps: &mut [LorentzVector<f64>],
    ) -> (f64, (f64, Option<f64>), (f64, Option<f64>)) {
        // check if we got provided the right number of random variables
        debug_assert_eq!(x.len(), self.get_dimensions());
        debug_assert_eq!(ps.len(), self.masses.0.len() + self.masses.1.len());
        debug_assert_eq!(self.masses.0.len(), 2); // we only support 2 -> n

        // move away from the border points
        for (rr, xr) in self.r.iter_mut().zip(x) {
            *rr = xr
                .max(FlatPhaseSpaceGenerator::EPSILON_BORDER)
                .min(1. - FlatPhaseSpaceGenerator::EPSILON_BORDER);
        }

        // Generate the ISR collinear factorization convolution variables xi<i> if
        // necessary. In order for the + distributions of the PDF counterterms and integrated
        // collinear ISR counterterms to hit the PDF only (and not the matrix elements or
        // observables functions), a change of variable is necessary: xb_1' = xb_1 * xi1
        let ext_start_index = x.len() - (self.masses.1.len() * 3 - 4);
        let (xi1, xi2) = if self.correlated_beam_convolution {
            // Both xi1 and xi2 must be set equal then
            (
                Some(self.r[ext_start_index - 1]),
                Some(self.r[ext_start_index - 1]),
            )
        } else {
            match self.is_beam_factorization_active {
                (true, true) => (
                    Some(self.r[ext_start_index - 2]),
                    Some(self.r[ext_start_index - 1]),
                ),
                (false, true) => (None, Some(self.r[ext_start_index - 1])),
                (true, false) => (Some(self.r[ext_start_index - 1]), None),
                (false, false) => (None, None),
            }
        };

        let mut wgt = 1.;
        let (E_cm, xb_1, xb_2) = match self.beam_type {
            (0, 0) => (self.collider_energy, 1., 1.),
            (1, 1) | (-1, 1) | (1, -1) | (-1, -1) => {
                let mut PDF_ycm = self.r[0];

                // We generate the PDF from two variables \tau = x1*x2 and ycm = 1/2 * log(x1/x2), so that:
                //  x_1 = sqrt(tau) * exp(+ycm)
                //  x_2 = sqrt(tau) * exp(-ycm)
                // The jacobian of this transformation is 1.
                let tot_final_state_masses = self.masses.1.iter().sum::<f64>();
                if tot_final_state_masses > self.collider_energy {
                    panic!("Collider energy is not large enough, there is no phase-space left.");
                }

                // Keep a hard cut at 1 GeV, which is the default for absolute_Ecm_min
                let tau_min = (tot_final_state_masses
                    .max(FlatPhaseSpaceGenerator::ABSOLUTE_E_CM_MIN)
                    / self.collider_energy)
                    .powi(2);
                let tau_max = 1.0;

                let PDF_tau = if self.masses.0.len() == 2 && self.masses.1.len() == 1 {
                    // Account for the \delta(xb_1*xb_2*s - m_h**2) and corresponding y_cm matching to unit volume
                    wgt /= self.collider_energy.powi(2);
                    // Here tau is fixed by the \delta(xb_1*xb_2*s - m_h**2) which sets tau to
                    tau_min
                } else {
                    // Including the corresponding Jacobian
                    wgt *= (tau_max - tau_min);
                    // Rescale tau appropriately
                    tau_min + (tau_max - tau_min) * self.r[1]
                };

                // And we can now rescale ycm appropriately
                let ycm_min = 0.5 * PDF_tau.ln();
                let ycm_max = -ycm_min;
                PDF_ycm = ycm_min + (ycm_max - ycm_min) * PDF_ycm;
                // and account for the corresponding Jacobian
                wgt *= (ycm_max - ycm_min);

                let xb_1 = PDF_tau.sqrt() * PDF_ycm.exp();
                let xb_2 = PDF_tau.sqrt() * -PDF_ycm.exp();
                // /!\ The mass of initial state momenta is neglected here.
                let E_cm = (xb_1 * xb_2).sqrt() * self.collider_energy;
                (E_cm, xb_1, xb_2)
            }
            _ => unimplemented!("Unknown beam configuration"),
        };

        wgt *= self.generate(E_cm, &self.r[ext_start_index..], &mut ps[2..]);

        // set the initial state
        // TODO: only do once?
        if self.masses.0[0] == 0. || self.masses.0[1] == 0. {
            ps[0] = LorentzVector::from_args(E_cm / 2.0, 0., 0., E_cm / 2.0);
            ps[1] = LorentzVector::from_args(E_cm / 2.0, 0., 0., -E_cm / 2.0);
        } else {
            let E_cm_sq = E_cm.powi(2);
            let M1sq = self.masses.0[0].powi(2);
            let M2sq = self.masses.0[1].powi(2);
            let E1 = (E_cm * E_cm + M1sq - M2sq) / E_cm;
            let E2 = (E_cm * E_cm - M1sq + M2sq) / E_cm;
            let Z = (E_cm_sq * E_cm_sq - 2. * E_cm_sq * M1sq - 2. * E_cm_sq * M2sq + M1sq * M1sq
                - 2. * M1sq * M2sq
                + M2sq * M2sq)
                .sqrt()
                / E_cm;
            ps[0] = LorentzVector::from_args(E1 / 2.0, 0., 0., Z / 2.0);
            ps[1] = LorentzVector::from_args(E2 / 2.0, 0., 0., -Z / 2.0);
        }

        (wgt, (xb_1, xi1), (xb_2, xi2))
    }
}
