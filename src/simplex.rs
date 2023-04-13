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

struct Matrix {
    pub nrows: usize,
    pub ncols: usize,
    pub v: Vec<f64>,
}

impl Matrix {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        let mut m = Matrix {
            nrows,
            ncols,
            v : Vec::<f64>::with_capacity(nrows * ncols),
        };
        m.v.resize(nrows * ncols, 0.0);
        m
    }

    #[inline] pub fn vpos(&self, row: usize, col: usize) -> usize {
        row*self.ncols + col
    }

    #[inline] pub fn get(&self, row: usize, col: usize) -> f64 {
        self.v[self.vpos(row, col)]
    }

    #[inline] pub fn set(&mut self, row: usize, col: usize, new_val: f64) -> &mut Self {
        let vpos = self.vpos(row, col);
        self.v[vpos] = new_val;
        self
    }

    #[inline] pub fn swap(&mut self, row1: usize, col1: usize, row2: usize, col2: usize) {
        let vpos1 = self.vpos(row1, col1);
        let vpos2 = self.vpos(row2, col2);
        self.v.swap(vpos1, vpos2);
    }

    #[inline] pub fn get_psum(&self, psum: &mut [f64]) {
        #[allow(clippy::needless_range_loop)]
        for j in 0..self.ncols {
            let mut sum = 0.0;
            for i in 0..self.nrows {
                sum += self.get(i, j);
            }
            psum[j] = sum;
        }
    }
}

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
///
/// Returned is the location of the minimum.
///
pub fn amoeba<F: Fn (&[f64]) -> f64>(
    fun: F,
    point: &[f64],
    step_delta: f64,
    ftol: f64,
    max_iterations: usize
) -> (Vec<f64>, f64, usize)
{
    const MIN_TOLERANCE: f64 = 1.0e-10_f64; // can be as small as f64 precision
    let ftol = ftol.max(MIN_TOLERANCE);

    let ndim = point.len();
    let mut dels = Vec::<f64>::new();
    dels.resize(ndim, step_delta);

    const TINY: f64 = 1.0e-10;

    let mut p = Matrix::new(ndim+1, ndim);
    #[allow(clippy::needless_range_loop)]
    for i in 0..ndim+1 {
        for j in 0..ndim {
            p.set(i, j, point[j]);
            if i != 0 { let x = p.get(i, i-1) + dels[i-1]; p.set(i, i-1, x); }
        }
    }

    let mpts = ndim + 1;
    let mut y = Vec::<f64>::new();
    y.resize(ndim + 1, 0.0);

    let mut psum = Vec::<f64>::new();
    psum.resize(ndim , 0.0);
    let mut pmin = Vec::<f64>::new();
    pmin.resize(ndim , 0.0);
    let mut x = Vec::<f64>::new();
    x.resize(ndim , 0.0);
    let mut ptry = Vec::<f64>::new();
    ptry.resize(ndim , 0.0);

    #[allow(clippy::needless_range_loop)]
    for i in 0..mpts {
        for j in 0..ndim {
            x[j] = p.get(i, j);
        }
        y[i] = fun(&x);
    }

    let mut fmin: f64 = y[0];

    p.get_psum(&mut psum);

    let mut nr_iterations: usize = 1;

    for _i in 0..max_iterations {
        let mut ilo = 0;
        // First we must determine which point is the highest (worst), next-highest, and
        // lowest (best), by looping over the points in the simplex.
        let mut ihi  = if y[0] > y[1] { 0 } else { 1 };
        let mut inhi = if y[0] > y[1] { 1 } else { 0 };

        for i in 0..mpts {
            if y[i] <= y[ilo] {
                ilo = i;
            }
            if y[i] > y[ihi] {
                inhi = ihi;
                ihi = i;
            }
            else if y[i] > y[inhi] && i != ihi {
                inhi = i;
            }
        }

        let rtol = 2.0 * (y[ihi] - y[ilo]).abs()
            / (y[ihi].abs() + y[ilo].abs() + TINY);

        fmin = y[0];

        // Compute the fractional range from highest to lowest and return if satisfactory.
        if rtol < ftol {
            //If returning, put best point and value in slot 0.
            y.swap(0, ilo);
            #[allow(clippy::needless_range_loop)]
            for i in 0..ndim {
                p.swap(0, i, ilo, i);
                pmin[i] = p.get(0, i);
            }
            //fmin = y[0];
            break;
        }



        //nfunc += 2;

        // Begin a new iteration. First extrapolate by a factor -1 through the face of the
        // simplex across from the high point, i.e., reflect the simplex from the high point.
        let mut ytry = amoeba_try(&mut p, &mut y, &mut psum, ihi, -1.0, &fun, &mut ptry);

        if ytry <= y[ilo] {
            // Gives a result better than the best point, so try an additional extrapolation
            // by a factor 2.
            /*ytry =*/ amoeba_try(&mut p, &mut y, &mut psum, ihi, 2.0, &fun, &mut ptry);
        }
        else if ytry >= y[inhi] {
            // The reflected point is worse than the second-highest, so look for an intermediate
            // lower point, i.e., do a one-dimensional contraction.
            let ysave = y[ihi];
            ytry = amoeba_try(&mut p, &mut y, &mut psum, ihi, 0.5, &fun, &mut ptry);
            if ytry >= ysave {
                // Can’t seem to get rid of that high point.
                // Better contract around the lowest (best) point.
                #[allow(clippy::needless_range_loop)]
                for i in 0..mpts {
                    if i != ilo {
                        for j in 0..ndim {
                            psum[j] = 0.5 * (p.get(i, j) + p.get(ilo, j));
                            p.set(i, j, psum[j]);
                        }
                        y[i] = fun(&psum);
                    }
                }
                //nfunc += ndim; // Keep track of function evaluations.
                p.get_psum(&mut psum); // Recompute psum.
            }
        }
        //else {
        //    --nfunc; // Correct the evaluation count.
        //}

        nr_iterations += 1;
    }

    (pmin, fmin, nr_iterations)
}

// Helper function: Extrapolates by a factor fac through the face of the simplex across from
// the high point, tries it, and replaces the high point if the new point is better.
fn amoeba_try<F: Fn (&[f64]) -> f64>(
    p: &mut Matrix,
    y: &mut [f64],
    psum: &mut [f64],
    ihi: usize,
    fac: f64,
    fun: F,
    ptry: &mut [f64] // size ndim
) -> f64
{
    let ndim = p.ncols;

    let fac1 = (1.0 - fac) / (ndim as f64);
    let fac2 = fac1 - fac;

    for j in 0..ndim {
        ptry[j] = psum[j] * fac1 - p.get(ihi, j) * fac2;
    }

    let ytry = fun(ptry); // Evaluate the function at the trial point.

    // If it’s better than the highest, then replace the highest.
    if ytry < y[ihi] {
        y[ihi] = ytry;
        for j in 0..ndim {
            psum[j] += ptry[j] - p.get(ihi, j);
            p.set(ihi, j, ptry[j]);
        }
    }

    ytry
}

#[cfg(test)]
#[test]
fn test_x2_y4() {
    // Since we know that `f(x,y) = x^2 + y^4` is always positive except at the original where it is zero,
    // we conclude that a local minimum occurs at (0, 0).
    fn x2_y4(x: &[f64]) -> f64 {
        x[0]*x[0] + x[1]*x[1]*x[1]*x[1]
    }

    let (min, fmin, nr_iterations) = amoeba(x2_y4, &[100.0, -100.0], 1.0, 1.0e-8, 100);

    println!("min: {}, {} fmin: {fmin} iterations: {nr_iterations}", min[0], min[1]);

    assert_float_absolute_eq!(min[0], 0.0, 1.0e-4);
    assert_float_absolute_eq!(min[1], 0.0, 1.0e-4);
}

#[cfg(test)]
#[test]
fn test_x2_y2_xy() {
    fn x2_y4_xy(x: &[f64]) -> f64 {
        x[0]*x[0] + x[1]*x[1] - 2.0*x[0]
    }

    let (min, fmin, nr_iterations) = amoeba(x2_y4_xy, &[10.0, 10.0], 0.1, 1.0e-9, 100);

    println!("min: {}, {} fmin: {fmin} iterations: {nr_iterations}", min[0], min[1]);

    assert_float_absolute_eq!(min[0], 1.0, 1.0e-4);
    assert_float_absolute_eq!(min[1], 0.0, 1.0e-4);
}


//https://www.gnu.org/software/gsl/doc/html/multimin.html