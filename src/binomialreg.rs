// use std::iter::zip;


fn logfactorial(x: usize) -> f64{
    (1..(x+1)).map(|q|(q as f64).ln()).sum()
}

fn log_binomial_coeff(N: usize, k: usize) -> f64{
    logfactorial(N) - logfactorial(k) - logfactorial(N-k)
}

fn binomial_loglike(x: usize, N: usize, p: f64) -> f64{
    // for a single datapoint
    let logcoeff = log_binomial_coeff(N, x);
    let loglike = (x as f64) * p.ln() + ((N-x) as f64) * (1.0-p).ln() + logcoeff;
    loglike
}

use itertools::izip;

use crate::utils;
pub fn phantom_binomial_regression(z: &[usize], mr: &[usize], r:&[usize]) -> (f64, Vec<f64>, Vec<f64>){
    // z is the number of non-chimeric molecules at aplification r
    // mr is the total number of mulecules at amplification r
    // amplification r
    assert_eq!(z.len() , r.len());

    let prange: Vec<f64> = (0..1000).map(|x| (x as f64)/1000.0).collect();

    let mut loglike_range: Vec<f64> = Vec::new();
    for p in prange.iter(){

        // each (z,mr,r) tuple corresponds to a z = Binomial(N=mr , p^r)
        let mut loglike = 0.0;
        for (zi, mri, ri) in izip!(z, mr, r){
            let ri_f64 = *ri as f64;
            loglike += binomial_loglike(*zi, *mri, p.powf(ri_f64));
        }
        loglike_range.push(loglike);
    }

    let (ix_max, _loglike_max) = utils::argsort::argmax_float(&loglike_range);
    let pmax = prange[ix_max];
    (pmax, prange, loglike_range)
}

#[cfg(test)]
mod test{
    use statrs::assert_almost_eq;

    use super::{logfactorial, log_binomial_coeff, binomial_loglike};
    #[test]
    fn test_logfac(){

        let tolerance = 0.000000001;
        assert_almost_eq!(logfactorial(0), 0.0, tolerance);
        assert_almost_eq!(logfactorial(1), 0.0, tolerance);
        assert_almost_eq!(logfactorial(2), 2_f64.ln(), tolerance);
        assert_almost_eq!(logfactorial(3), 6_f64.ln(), tolerance);
    }

    #[test]
    fn test_log_binomial_coeff(){
        let tolerance = 0.000000001;
        
        assert_almost_eq!(log_binomial_coeff(5, 1), 5_f64.ln(), tolerance);
        assert_almost_eq!(log_binomial_coeff(5, 5), 1_f64.ln(), tolerance);
    }

    #[test]
    fn test_binomial_loglike(){
        assert_eq!(
            binomial_loglike(1, 1, 0.5 ),
            0.5_f64.ln()
        )
    }
}