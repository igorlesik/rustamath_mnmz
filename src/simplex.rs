//! Downhill Simplex Method in Multidimensions.
//!
//! (c) 2023 Igor Lesik
//! MIT license
//!
//! References:
//!
//! 1. William H. Press - Numerical recipes, the art of scientific computing.
//!   Cambridge University Press (2007).
//!


/// Downhill Simplex Method in Multidimensions.
///
/// References:
///
/// 1. William H. Press - Numerical recipes, the art of scientific computing.
///   Cambridge University Press (2007).
///
/// Multidimensional minimization of the function `fun(x)`, where `x[0..ndim-1]`
/// is a vector in `ndim` dimensions, by the downhill simplex method of Nelder and Mead.
/// The initial simplex is specified as in equation `Pi = P0 + de` by a `point[0..ndim-1]` and a
/// constant displacement `step_delta` along each coordinate direction.
/// Returned is the location of the minimum.
///
pub fn amoeba<F: Fn (&[f64]) -> f64>(
    _fun: F,
    point: &[f64],
    step_delta: f64,
    max_iterations: usize
) -> f64
{
    let ndim = point.len();
    let mut dels = Vec::<f64>::new();
    dels.resize(ndim, step_delta);

    const TINY: f64 = 1.0e-10;

    let mut y = Vec::<f64>::new();
    y.resize(ndim + 1, 0.0);

    for _i in 0..max_iterations {
        let ilo = 0;
        // First we must determine which point is the highest (worst), next-highest, and
        // lowest (best), by looping over the points in the simplex.
        let ihi  = if y[0] > y[1] { 0 } else { 1 };
        let inhi = if y[0] > y[1] { 1 } else { 0 };
    }

    0.0
}



//test
//since we know that f(x,y) = x2 + y4 is always positive
//except at the original where it is zero,
//we conclude that a local minimum occurs at (0, 0).

//https://real-statistics.com/other-mathematical-topics/function-maximum-minimum/local-maxima-minima-multivariate/

//https://www.gnu.org/software/gsl/doc/html/multimin.html