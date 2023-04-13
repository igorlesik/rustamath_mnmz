# Rustamath. Library of minimization functions.

![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)
[![CI](https://github.com/igorlesik/rustamath_mnmz/actions/workflows/test.yml/badge.svg)](https://github.com/igorlesik/rustamath_mnmz/actions/workflows/test.yml)
[![crates.io version][crates-io-shields]][crates-io]
[![docs.rs][docs-rs-shields]][docs-rs]

[crates-io]: https://crates.io/crates/rustamath_mnmz
[crates-io-shields]: https://img.shields.io/crates/v/rustamath_mnmz.svg
[docs-rs]: https://docs.rs/rustamath_mnmz
[docs-rs-shields]: https://img.shields.io/badge/docs.rs-rustdoc-green.svg


Task of minimization: for given function _f_ that depends on one or more independent
variables, find the value of those variables where _f_ takes on a minimum value.

Supported methods:

- One Dimension
  * Bracket a Minimum
  * Golden Section Search
  * Brent’s Method
  * Brent’s Method using First Derivative
- Multidimensions
  * Downhill Simplex Method

## Example of Downhill Simplex search

```rust
fn test_paraboloid() {
    //  Paraboloid center at (1,2), scale factors (10, 20), minimum value 30
    let p = vec![1.0, 2.0, 10.0, 20.0, 30.0];

    let paraboloid = |x: &[f64]|  {
        // Paraboloid centered on (p[0],p[1]), with scale factors (p[2],p[3]) and minimum p[4]
        p[2] * (x[0] - p[0]) * (x[0] - p[0]) + p[3] * (x[1] - p[1]) * (x[1] - p[1]) + p[4]
    };

    let (min, fmin, nr_iterations) = amoeba(paraboloid, &[100.0, -100.0], 1.1, 1.0e-9, 100);

    println!("min: {}, {} fmin: {fmin} iterations: {nr_iterations}", min[0], min[1]);

    assert_float_absolute_eq!(min[0], 1.0, 1.0e-4);
    assert_float_absolute_eq!(min[1], 2.0, 1.0e-4);
    assert_float_absolute_eq!(fmin,  30.0, 1.0e-4);
}
```

Output:

```console
min: 0.999933263302534, 1.9999850642280714 fmin: 30.000000072002226 iterations: 78
```

## Example of Brent’s Method Using First Derivative

```rust
fn test_cosine() {
    use super::{golden_section_search, brent_search};

    // Minimum at Pi when x ∈ [0, 2*Pi].
    let cosine = |x: f64| (x.cos(), -(x.sin()));

    let ranges = vec![(0.01, 1.0)];

    for range in ranges {
        let (xmin, f, nr_iterations) =
            brent_df_search(cosine, range.0, range.1, 0.0, 0);

        let (xmin_golden, _, nr_iterations_golden) =
            golden_section_search(|x| cosine(x).0, range.0, range.1, 0.0, 0);

        let (xmin_brent, _, nr_iterations_brent) =
            brent_search(|x| cosine(x).0, range.0, range.1, 0.0, 0);

        println!("xmin: {:.8} f(xmin): {:6.2} iterations: {} vs brent {} vs golden {}",
            xmin, f, nr_iterations, nr_iterations_brent, nr_iterations_golden
        );

        assert_float_relative_eq!(xmin, std::f64::consts::PI, 1.0e-8);
        assert_float_relative_eq!(xmin_brent, std::f64::consts::PI, 1.0e-8);
        assert_float_relative_eq!(xmin_golden, std::f64::consts::PI, 1.0e-8);
    }
}
```

Output:

```console
xmin: 3.14159265 f(xmin):  -1.00 iterations: 4 vs brent 8 vs golden 36
```