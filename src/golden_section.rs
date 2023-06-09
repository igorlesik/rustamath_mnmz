//! Golden section search for a minimum.
//!
//! (c) 2023 Igor Lesik
//! MIT license
//!
//! References:
//!
//! 1. William H. Press - Numerical recipes, the art of scientific computing.
//!   Cambridge University Press (2007).
//!
use super::bracket::{find_bracket, shft3, shft2};

/// Smallest tolerance.
///
/// See book "Numerical recipes, the art of scientific computing."
/// sqrt(f64 precision 10^16), by Tailor series for `f(x+eps)`
const MIN_TOLERANCE: f64 = 3.0e-8_f64;

/// Golden section search for a minimum.
///
/// - William H. Press - Numerical recipes, the art of scientific computing.
///   Cambridge University Press (2007).
///
/// Given a function `f`, and given a bracketing triplet of abscissas a, b, c
/// (such that b is between a and c, and f(x) is less than both f(a) and f(c)),
/// this routine performs a golden section search for the minimum,
/// isolating it to a fractional precision of about `tolerance`.
///
/// # Example
///
/// ```
/// use rustamath_mnmz::golden_section_search;
/// use assert_float_eq::*;
/// // Roots 1.0 and 2.0, minimum at 1.5.
/// let poly2 = |x: f64| (x-1.0)*(x-2.0);
/// let ranges = vec![(10.0, 20.0), (20.0, 10.0), (-10.0, 0.0)];
/// for range in ranges {
///     let (xmin, f, nr_iterations) = golden_section_search(poly2, range.0, range.1, 0.0, 0);
///     println!("MIN: {:.8} f(xmin): {:6.2} iterations:{}",
///         xmin, f, nr_iterations
///     );
///     assert_float_relative_eq!(xmin, 1.5, 1.0e-8);
/// }
///
pub fn golden_section_search<F: Fn (f64) -> f64>(
    fun: F,
    a: f64,
    b: f64,
    tol: f64,
    max_iterations: usize
) -> (f64, f64, usize)
{
    let tol = tol.max(MIN_TOLERANCE);
    let max_iterations = if max_iterations < 1 { 500 } else { max_iterations.min(1000) };
    const R: f64 = 0.61803399_f64;
    const C: f64 = 1.0 - R; // The golden ratios.

    let bracket = find_bracket(&fun, a, b);
    let a = bracket.a;
    let b = bracket.b;
    let c = bracket.c;

    // At any given time we will keep track of four points, x0,x1,x2,x3.
    let mut x1: f64;
    let mut x2: f64;
    let mut x0 = a;
    let mut x3 = c;

    // Make x0 to x1 the smaller segment, and fill in the new point to be tried.
    if (c-b).abs() > (b-a).abs() {
        x1 = b;
        x2 = b + C*(c-b);
    } else {
        x2 = b;
        x1 = b - C*(b-a);
    }

    // The initial function evaluations. Note that we never need to evaluate
    // the function at the original endpoints.
    let mut f1 = fun(x1);
    let mut f2 = fun(x2);
    let mut nr_iterations: usize = 0;

    while (x3-x0).abs() > tol*(x1.abs() + x2.abs()) {
        if f2 < f1 {
            let d = R*x2 + C*x3;
            shft3(&mut x0, &mut x1, &mut x2, d);
            shft2(&mut f1, &mut f2, fun(x2));
        }
        else {
            let d = R*x1 + C*x0;
            shft3(&mut x3, &mut x2, &mut x1, d);
            shft2(&mut f2, &mut f1, fun(x1));
        }
        nr_iterations += 1;

        if nr_iterations >= max_iterations { break; }

        // @igor: saw/non-smooth functions demostrate that `tol*(x1.abs() + x2.abs())`
        // gets smaller faster than `(x3-x0).abs()` preventing the conversion;
        // here if we see that x3 is close to x0 and f1 to f2 we force the exit.
        if nr_iterations > 10 && (x3-x0).abs() < tol && (f1 - f2).abs() < tol { break; }
    }

    // Output the best of the two current values.
    if f1 < f2 {
        (x1, f1, nr_iterations)
    }
    else {
        (x2, f2, nr_iterations)
    }
}

#[cfg(test)]
#[test]
fn test_poly2() {
    // Roots 1.0 and 2.0, minimum at 1.5.
    let poly2 = |x: f64| (x-1.0)*(x-2.0);

    let ranges = vec![(10.0, 20.0), (20.0, 10.0), (-10.0, 0.0),
        (-2000.0, -1000.0), (-10_000.0, 30_000.0), (0.0001, 0.0002), (-0.00001, 1.4999)];

    for range in ranges {
        let (xmin, f, nr_iterations) = golden_section_search(poly2, range.0, range.1, 0.0, 0);

        println!("MIN: {:.8} f(xmin): {:6.2} iterations:{}",
            xmin, f, nr_iterations
        );

        assert_float_relative_eq!(xmin, 1.5, 1.0e-8);
    }
}

#[cfg(test)]
#[test]
fn test_cosine() {
    // Minimum at Pi on [0, 2*Pi].
    let cosine = |x: f64| x.cos();

    let ranges = vec![(0.01, 1.0)];

    for range in ranges {
        let (xmin, f, nr_iterations) = golden_section_search(cosine, range.0, range.1, 0.0, 0);

        println!("MIN: {:.8} f(xmin): {:6.2} iterations:{}",
            xmin, f, nr_iterations
        );

        assert_float_relative_eq!(xmin, std::f64::consts::PI, 1.0e-8);
    }
}

#[cfg(test)]
#[test]
fn test_saw() {
    let saw = |x: f64| if  x >= 0.0 { x*x*x } else { -x / 1000.0 } ;

    let ranges = vec![(10.0, 20.0), (20.0, 10.0), (-10.0, 0.0),
        (-2000.0, -1000.0), (-10_000.0, 30_000.0), (0.0001, 0.0002), (-0.00001, 1.4999)];

        for range in ranges {
            let (xmin, f, nr_iterations) = golden_section_search(saw, range.0, range.1, 1.0e-5, 0);

            println!("MIN: {:.8} f(xmin): {:6.2} iterations:{}",
                xmin, f, nr_iterations
            );

            assert_float_absolute_eq!(xmin, 0.0, 1.0e-5);
        }
}