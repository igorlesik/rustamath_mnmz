//! Brent's method with First Derivative for searching for a minimum.
//!
//! (c) 2023 Igor Lesik
//! MIT license
//!
//! References:
//!
//! 1. William H. Press - Numerical recipes, the art of scientific computing.
//!   Cambridge University Press (2007).
//!
use super::bracket::{find_bracket, mov3};

/// Smallest tolerance.
///
/// See book "Numerical recipes, the art of scientific computing."
/// sqrt(f64 precision 10^16), by Tailor series for `f(x+eps)`
const MIN_TOLERANCE: f64 = 3.0e-8_f64;

/// Brent's method to search for a minimum.
///
/// - William H. Press - Numerical recipes, the art of scientific computing.
///   Cambridge University Press (2007).
///
///
/// # Example
///
/// ```
/// use rustamath_mnmz::{brent_df_search, brent_search, golden_section_search};
/// use assert_float_eq::*;
/// let cosine = |x: f64| (x.cos(), -(x.sin())); // Minimum at Pi when x ∈ [0, 2*Pi].
/// let range = (0.01, 1.0);
///
/// let (xmin, f, nr_iterations) = brent_df_search(cosine, range.0, range.1, 0.0, 0);
/// let (_, _, nr_iterations_brent) = brent_search(|x| cosine(x).0, range.0, range.1, 0.0, 0);
/// let (_, _, nr_iterations_golden) = golden_section_search(|x| cosine(x).0, range.0, range.1, 0.0, 0);
///
/// println!("xmin: {:.8} f(xmin): {:6.2} iterations: {} vs brent {} golden {}",
///     xmin, f, nr_iterations, nr_iterations_brent, nr_iterations_golden);
///
/// assert_float_relative_eq!(xmin, std::f64::consts::PI, 1.0e-8);
/// ```
pub fn brent_df_search<F: Fn (f64) -> (f64, f64)>(
    fun: F,
    a: f64,
    b: f64,
    tol: f64,
    max_iterations: usize
) -> (f64, f64, usize)
{
    let tol = tol.max(MIN_TOLERANCE);
    let max_iterations = if max_iterations < 1 { 500 } else { max_iterations.min(1000) };

       // ZEPS is a small number that protects against trying to achieve
    // fractional accuracy for a minimum that happens to be exactly zero.
    // https://doc.rust-lang.org/std/primitive.f64.html#associatedconstant.EPSILON
    const ZEPS: f64 = f64::EPSILON * 1.0e-3;

    let bracket = find_bracket(|x| fun(x).0, a, b);
    let ax = bracket.a;
    let _b = bracket.b;
    let c = bracket.c;

    // a and b must be in ascending order, but input abscissas need not be.
    let mut a = if ax < c { ax } else { c };
    let mut b = if ax > c { ax } else { c };

    // This will be the distance moved on the step before last.
    let mut e: f64 = 0.0;
    let mut d: f64 = 0.0;

    let mut x = b; let mut w = b; let mut v = b;

    let (mut fx, mut dx) = fun(x);
    let mut fw = fx;
    let mut fv = fx;
    let mut dw = dx;
    let mut dv = dx;

    let mut nr_iterations: usize = 0;

    for _i in 0..max_iterations {
        // test if we done
        let xm = 0.5 * (a+b);
        let tol1 = tol * x.abs() + ZEPS;
        let tol2 = 2.0 * (tol1 + ZEPS);

        if (x - xm).abs() <= (tol2 - 0.5*(b - a)) {
            break;
        }

        if e.abs() > tol1 {
            let mut d1 = 2.0 * (b-a); // Initialize these d's to an out-of-bracket value.
            let mut d2 = d1;
            if dw != dx { d1 = (w-x)*dx/(dx-dw); } // Secant method with one point.
            if dv != dx { d2 = (v-x)*dx/(dx-dv); }

            // Which of these two estimates of d shall we take?
            // We will insist that they be within the bracket,
            // and on the side pointed to by the derivative at x:
            let u1 = x + d1;
            let u2 = x + d2;
            let ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
            let ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;

            let olde = e; // Movement on the step before last.
            e = d;
            // Take only an acceptable d, and if both are acceptable, then take
            // the smallest one.
            if ok1 || ok2 {
                if ok1 && ok2 {
                    d = if d1.abs() < d2.abs() { d1 } else { d2 };
                } else if ok1 {
                    d = d1;
                } else {
                    d = d2;
                }

                if d.abs() <= (0.5*olde).abs() {
                    let u = x + d;
                    if u-a < tol2 || b-u < tol2 {
                        d = tol1.copysign(xm-x);
                    }
                }
                else { // Bisect, not golden section.
                    // Decide which segment by the sign of the derivative.
                    e = if dx >= 0.0 { a-x } else { b-x };
                    d = 0.5 * e;
                }
            }
            else {
                e = if dx >= 0.0 { a-x } else { b-x };
                d = 0.5 * e;
            }
        }
        else {
            e = if dx >= 0.0 { a-x } else { b-x };
            d = 0.5 * e;
        }

        let u: f64;
        let fu: f64;

        if d.abs() >= tol1 {
            u = x + d;
            (fu, _) = fun(u);
        }
        else {
            u = x + tol1.copysign(d);
            (fu, _) = fun(u);
            // If the minimum step in the downhill direction takes us uphill,
            // then we are done.
            if fu > fx  {
                break;
            }
        }

        let du: f64;

        (_, du) = fun(u);
        if fu <= fx {
            if u >= x { a = x; } else { b = x; }
            mov3(&mut v, &mut fv, &mut dv, w, fw, dw);
            mov3(&mut w, &mut fw, &mut dw, x, fx, dx);
            mov3(&mut x, &mut fx, &mut dx, u, fu, du);
        }
        else {
            if u < x { a = u; } else { b = u; }
            if fu <= fw || w == x {
                mov3(&mut v, &mut fv, &mut dv, w, fw, dw);
                mov3(&mut w, &mut fw, &mut dw, u, fu, du);
            }
            else if fu < fv || v == x || v == w {
                mov3(&mut v, &mut fv, &mut dv, u, fu, du);
            }
        }
    
        nr_iterations += 1;
    }

    (x, fx, nr_iterations)
}

#[cfg(test)]
#[test]
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

        println!("xmin: {:.8} f(xmin): {:6.2} iterations: {} vs brent {} golden {}",
            xmin, f, nr_iterations, nr_iterations_brent, nr_iterations_golden
        );

        assert_float_relative_eq!(xmin, std::f64::consts::PI, 1.0e-8);
        assert_float_relative_eq!(xmin_brent, std::f64::consts::PI, 1.0e-8);
        assert_float_relative_eq!(xmin_golden, std::f64::consts::PI, 1.0e-8);
    }
}

#[cfg(test)]
#[test]
fn test_poly2() {
    use super::{golden_section_search, brent_search};

    // Roots 1.0 and 2.0, minimum at 1.5.
    let poly2 = |x: f64| ((x-1.0)*(x-2.0), 2.0*x-3.0);

    let ranges = vec![(10.0, 20.0), (20.0, 10.0), (-10.0, 0.0),
        (-2000.0, -1000.0), (-10_000.0, 30_000.0), (0.0001, 0.0002), (-0.00001, 1.4999)];

    for range in ranges {
        let (xmin, f, nr_iterations) = brent_df_search(poly2, range.0, range.1, 0.0, 0);

        let (xmin_golden, _, nr_iterations_golden) =
            golden_section_search(|x| poly2(x).0, range.0, range.1, 0.0, 0);

        let (xmin_brent, _, nr_iterations_brent) =
            brent_search(|x| poly2(x).0, range.0, range.1, 0.0, 0);

        println!("xmin: {:.8} f(xmin): {:6.2} iterations: {} vs brent {} vs golden {}",
            xmin, f, nr_iterations, nr_iterations_brent, nr_iterations_golden
        );

        assert_float_relative_eq!(xmin, 1.5, 1.0e-8);
        assert_float_relative_eq!(xmin_brent, 1.5, 1.0e-8);
        assert_float_relative_eq!(xmin_golden, 1.5, 1.0e-8);
    }
}
