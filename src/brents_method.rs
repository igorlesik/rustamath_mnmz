//! Brent's method for searching for a minimum.
//!
//! (c) 2023 Igor Lesik
//! MIT license
//!
//! References:
//!
//! 1. William H. Press - Numerical recipes, the art of scientific computing.
//!   Cambridge University Press (2007).
//!
use super::bracket::{FunToMnmz, find_bracket, shft3};

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
pub fn brent_search(
    fun: FunToMnmz,
    a: f64,
    b: f64,
    tol: f64,
    max_iterations: usize
) -> (f64, f64, usize)
{
    let tol = tol.max(MIN_TOLERANCE);
    let max_iterations = if max_iterations < 1 { 500 } else { max_iterations.min(1000) };
    const RGOLD: f64 = 0.61803399_f64;
    const CGOLD: f64 = 1.0 - RGOLD; // The golden ratios.

    // ZEPS is a small number that protects against trying to achieve
    // fractional accuracy for a minimum that happens to be exactly zero.
    // https://doc.rust-lang.org/std/primitive.f64.html#associatedconstant.EPSILON
    const ZEPS: f64 = f64::EPSILON * 1.0e-3;

    let bracket = find_bracket(fun, a, b);
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

    let mut fx = fun(x);
    let mut fw = fx;
    let mut fv = fx;

    let mut nr_iterations: usize = 0;

    for _i in 0..max_iterations {
        // test if we done
        let xm = 0.5 * (a+b);
        let tol1 = tol * x.abs() + ZEPS;
        let tol2 = 2.0 * (tol1 + ZEPS);

        if (x - xm).abs() <= (tol2 - 0.5*(b - a)) { break; }

        // @igor force exit
        if nr_iterations > 100 && (b - a).abs() < tol { break; }

        // Construct a trial parabolic fit.
        if e.abs() > tol1 {
            let r = (x-w)*(fx-fv);
            let q = (x-v)*(fx-fw);
            let p = (x-v)*q-(x-w)*r;
            let q = 2.0*(q-r);
            let p = if q > 0.0 { -p } else { p };
            let q = q.abs();
            let etemp = e;
            e = d;

            // determine the acceptability of the parabolic fit
            if p.abs() >= (0.5*q*etemp).abs() || p <= q*(a-x) || p >= q*(b-x) {
                // take the golden section step into the larger of the two segments.
                e = if x >= xm  { a-x } else { b-x };
                d = CGOLD * e;
            }
            else {
                d = p / q; // Take the parabolic step.
                let u = x + d;
                if (u-a) < tol2 || (b-u) < tol2 {
                    d = tol1.copysign(xm-x);
                }
            }
        }
        else {
            e = if x >= xm  { a-x } else { b-x };
            d = CGOLD * e;
        }

        let u = if d.abs() >= tol1 { x+d } else { x + tol1.copysign(d) };
        let fu = fun(u);

        if fu <= fx {
            if u >= x { a = x; } else { b = x; }
            shft3(&mut v, &mut w, &mut x, u);
            shft3(&mut fv, &mut fw, &mut fx, fu);
        }
        else {
            if u < x { a = u; } else { b = u; }
            if fu <= fw || w == x {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }

        nr_iterations += 1;
    }

    (x, fx, nr_iterations)
}

#[cfg(test)]
#[test]
fn test_poly2() {
    use super::golden_section::golden_section_search;

    // Roots 1.0 and 2.0, minimum at 1.5.
    let poly2 = |x: f64| (x-1.0)*(x-2.0);

    let ranges = vec![(10.0, 20.0), (20.0, 10.0), (-10.0, 0.0),
        (-2000.0, -1000.0), (-10_000.0, 30_000.0), (0.0001, 0.0002), (-0.00001, 1.4999)];

    for range in ranges {
        let (xmin, f, nr_iterations) = brent_search(poly2, range.0, range.1, 0.0, 0);

        let (_, _, nr_iterations_golden) = golden_section_search(poly2, range.0, range.1, 0.0, 0);

        println!("MIN: {:.8} f(xmin): {:6.2} iterations: {} vs golden {}",
            xmin, f, nr_iterations, nr_iterations_golden
        );

        assert_float_relative_eq!(xmin, 1.5, 1.0e-8);
    }
}

#[cfg(test)]
#[test]
fn test_cosine() {
    use super::golden_section::golden_section_search;

    // Minimum at Pi on [0, 2*Pi].
    let cosine = |x: f64| x.cos();

    let ranges = vec![(0.01, 1.0)];

    for range in ranges {
        let (xmin, f, nr_iterations) = brent_search(cosine, range.0, range.1, 0.0, 0);

        let (_, _, nr_iterations_golden) = golden_section_search(cosine, range.0, range.1, 0.0, 0);

        println!("MIN: {:.8} f(xmin): {:6.2} iterations: {} vs golden {}",
            xmin, f, nr_iterations, nr_iterations_golden
        );

        assert_float_relative_eq!(xmin, std::f64::consts::PI, 1.0e-8);
    }
}

#[cfg(test)]
#[test]
fn test_saw() {
    use super::golden_section::golden_section_search;

    let saw = |x: f64| if  x >= 0.0 { x*x*x } else { -x / 1000.0 } ;

    let ranges = vec![(10.0, 20.0), (20.0, 10.0), (-10.0, 0.0),
        (-2000.0, -1000.0), (-10_000.0, 30_000.0), (0.0001, 0.0002), (-0.00001, 1.4999)];

        for range in ranges {
            let (xmin, f, nr_iterations) = brent_search(saw, range.0, range.1, 1.0e-5, 0);

            let (xmin_gold, f_gold, nr_iterations_gold) = golden_section_search(saw, range.0, range.1, 1.0e-5, 0);

            println!("MIN: {:.8} f(xmin): {:6.2} iterations: {} vs golden {:.8} f(xmin): {:6.2} {}",
                xmin, f, nr_iterations, xmin_gold, f_gold, nr_iterations_gold
            );

            assert_float_absolute_eq!(xmin, 0.0, 1.0e-5);
        }
}