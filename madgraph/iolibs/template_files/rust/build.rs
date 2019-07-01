use std::env;

fn main() {
    let working_dir = env::current_dir().unwrap();
    let proc_root = env::var("PROC_ROOT").unwrap_or_else(
        |_| { println!("No PROC_ROOT defined. Using working directory"); 
            working_dir.to_str().unwrap().to_owned() 
    });

    println!("cargo:rustc-link-search={}/lib", proc_root);
    println!("cargo:rustc-link-search={}/lib/matrix_elements", proc_root);
    println!("cargo:rustc-link-lib=static=model");
    println!("cargo:rustc-link-lib=gfortran");
}
