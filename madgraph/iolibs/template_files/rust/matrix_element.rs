
mod %(matrix_element_id)sMatrixElementMod {
    use libc::{c_char, c_double, c_int, c_longlong, c_void};

    #[link(name = "%(matrix_element_lib)s", kind = "static")]
    extern "C" {
        #[link_name = "c_%(prefix)sme_accessor_hook"]
        pub fn c_me_accessor_hook(p: *const c_double, hel: *const c_int, user_alphas: *const c_double, ans: *mut c_double);
        #[link_name = "c_%(prefix)sinitialise"]
        pub fn c_initialise(p: *const c_char);
    }
}

pub struct %(matrix_element_id)sMatrixElement {}
impl MatrixElement for %(matrix_element_id)sMatrixElement {
    #[inline]
    unsafe fn c_me_accessor_hook(p: *const c_double, hel: *const c_int, user_alphas: *const c_double, ans: *mut c_double) {
        unsafe {
            %(matrix_element_id)sMatrixElementMod::c_me_accessor_hook(p, hel, user_alphas, ans);
        }
    }
    
    #[inline]
    unsafe fn c_initialise(p: *const c_char) {
        unsafe {
            %(matrix_element_id)sMatrixElementMod::c_initialise(p);
        }
    }
}
