fn main() {
    cc::Build::new()
        .file("archt/kbessel.c")
        .flag("-Wno-incompatible-pointer-types")
        .compile("archt-c");
}
