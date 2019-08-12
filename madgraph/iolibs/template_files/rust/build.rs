use std::env;

fn main() {
    let proc_root = env::var("PROC_ROOT").unwrap_or_else(
        |_| { panic!("The environment variable PROC_ROOT should be defined."); }
    );

    let madnklo_root = env::var("MADNKLO_ROOT").unwrap_or_else(
        |_| { panic!("The environment variable MADNKLO_ROOT should be defined."); }
    );

    println!("cargo:rustc-link-search={}/lib", proc_root);
    println!("cargo:rustc-link-search={}/lib/matrix_elements", proc_root);
    println!("cargo:rustc-link-search={}/lib/helas", proc_root);    
    println!("cargo:rustc-link-lib=static=model");
    println!("cargo:rustc-link-lib=gfortran");

    println!("cargo:rustc-link-search={}/vendor/fjcore", madnklo_root);
    println!("cargo:rustc-link-search={}/HEPTools/lhapdf6/lib", madnklo_root);
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=gcc");
}
