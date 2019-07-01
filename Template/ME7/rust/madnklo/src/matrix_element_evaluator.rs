use libc::{c_char, c_double, c_int, c_longlong, c_void};
use std::ffi::CString;
use vector::LorentzVector;

pub trait MatrixElement {
    unsafe fn c_me_accessor_hook(p: *const c_double, hel: *const c_int, user_alphas: *const c_double, ans: *mut c_double) where Self: Sized;
    unsafe fn c_initialise(p: *const c_char)  where Self: Sized;
}

pub struct MatrixElementEvaluator<T> where T: MatrixElement + Sized {
    external_momenta_plain: [[f64; 4]; 6],//MatrixElementEvaluator::NUM_EXTERNALS],
    matrix_element: T,
}

impl<T> MatrixElementEvaluator<T> where T: MatrixElement + Sized {
    pub const NUM_EXTERNALS: usize = 6;

    pub fn new(card_filename: &str, matrix_element: T) -> MatrixElementEvaluator<T> {
        let mut a = [' ' as i8; 512];
        for (xa, c) in a.iter_mut().zip(card_filename.chars()) {
            *xa = c as i8;
        }

        unsafe {
            T::c_initialise(&a[0] as *const c_char);
        }

        MatrixElementEvaluator {
            external_momenta_plain: [[0.; 4]; 6],//MatrixElementEvaluator::NUM_EXTERNALS],
            matrix_element,
        }
    }

    pub fn evaluate(&mut self, external_momenta: &[LorentzVector<f64>], helicity: i32, alpha_s: f64) -> f64 {
        for (e_p, e) in self.external_momenta_plain.iter_mut().zip(external_momenta) {
            e_p[0] = e.t;
            e_p[1] = e.x;
            e_p[2] = e.y;
            e_p[3] = e.z;
        }

        let mut ans = 0.0f64;
        unsafe {
            T::c_me_accessor_hook(
                &self.external_momenta_plain[0] as *const c_double,
                &helicity as *const c_int,
                &alpha_s as *const c_double,
                &mut ans as *mut c_double,
            );
        }
        ans
    }

}
