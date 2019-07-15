use std::env;

fn main() {
    let proc_root = env::var("PROC_ROOT").unwrap_or_else(
        |_| { panic!("The environment variable PROC_ROOT should be defined."); }
    );

    println!("cargo:rustc-link-search={}/lib", proc_root);
    println!("cargo:rustc-link-search={}/lib/matrix_elements", proc_root);
    println!("cargo:rustc-link-lib=static=model");
    println!("cargo:rustc-link-lib=gfortran");

    println!("cargo:rustc-link-search=/home/ben/Sync/Research/madnklo/HEPTools/lhapdf6/lib"); // TODO: make dynamic
}
