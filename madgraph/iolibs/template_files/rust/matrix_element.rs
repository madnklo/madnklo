
mod %(matrix_element_id)sMatrixElementMod {
    use libc::{c_char, c_double, c_int};

    #[link(name = "%(matrix_element_lib)s", kind = "static")]
    extern "C" {
        #[link_name = "c_%(prefix)sme_accessor_hook"]
        pub fn me_accessor_hook(p: *const c_double, hel: *const c_int, user_alphas: *const c_double, ans: *mut c_double);
        #[link_name = "c_%(prefix)sinitialise"]
        pub fn initialise(p: *const c_char);
    }
}

pub struct %(matrix_element_id)sMatrixElement {}
impl MatrixElement for %(matrix_element_id)sMatrixElement {
    #[inline]
    unsafe fn me_accessor_hook(p: *const c_double, hel: *const c_int, user_alphas: *const c_double, ans: *mut c_double) {
        %(matrix_element_id)sMatrixElementMod::me_accessor_hook(p, hel, user_alphas, ans);
    }
    
    #[inline]
    unsafe fn initialise(p: *const c_char) {
        %(matrix_element_id)sMatrixElementMod::initialise(p);
    }
}
