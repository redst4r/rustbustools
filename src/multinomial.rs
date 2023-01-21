use statrs::distribution::{Multinomial, Binomial};  // warning the statrs::Binomial has very slow sampling (sum of Bernullis)
use probability;  // use probability::distribution::Binomial instead, which does inverse cdf sampling
use probability::distribution::Sample;
use probability::prelude::*;
use rand;
use rand::distributions::Distribution;


pub fn multinomial_sample_statrs(n: u64, pvec: Vec<f64>) -> Vec<f64>{
    let mut r = rand::thread_rng();
    let mut x :Vec<f64> = Vec::new();

    // normalize the pvec
    let _sum: f64 = pvec.iter().sum();
    let pvec_norm: Vec<f64> = pvec.iter().map(|x| x/_sum).collect();
    let dim = pvec_norm.len();

    let mut remaining_p = 1.0;
    let mut remaining_n = n;
    for (counter, p) in pvec_norm.iter().enumerate(){
        // binomial with BIN(p/remaining_p, remaining_n)
        if remaining_n == 0{
            x.push(0.0);
        }
        else if counter == dim - 1{
            //last element, p will be 1 (or due to errors, a litte >1)
            x.push(remaining_n as f64);
        }
        else{
            let _ptmp = p / remaining_p;
            if !(0.0 < _ptmp && _ptmp < 1.0){
                println!("{:?}", &pvec_norm[counter..]);
                panic!("0<{}<1, counter {}", _ptmp, counter);
            }            
            let b = Binomial::new(_ptmp, remaining_n).expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);

            x.push(b);
            remaining_n -= b as u64;
            remaining_p -= p;
        }
    }
    x
}

pub fn multinomial_sample(n: u64, pvec: Vec<f64>) -> Vec<f64>{
    /*
    my own multinomial sampling, using the fact that all marginals are binomial

    statrs version does the same algorithm, but relies internally on a statrs::distribution::Binomial
    which is extremely slow.
    */
    let mut source = source::default(42);

    let mut x :Vec<f64> = Vec::new();

    // normalize the pvec
    let _sum: f64 = pvec.iter().sum();
    let pvec_norm: Vec<f64> = pvec.iter().map(|x| x/_sum).collect();
    let dim = pvec_norm.len();

    let mut remaining_p = 1.0;
    let mut remaining_n = n;
    for (counter, p) in pvec_norm.iter().enumerate(){
        // binomial with BIN(p/remaining_p, remaining_n)
        if remaining_n == 0{
            x.push(0.0);
        }
        else if counter == dim - 1{
            //lastt else if element, p will be 1 (or due to errors, a litte >1)
            x.push(remaining_n as f64);
        }
        else{
            let _ptmp = p / remaining_p;

            if !(0.0 < _ptmp && _ptmp < 1.0){
                println!("{:?}", &pvec_norm[counter..]);
                panic!("0<{}<1, counter {}", _ptmp, counter);
            }
            let btmp = probability::distribution::Binomial::new(remaining_n as usize, _ptmp ).sample(&mut source); //.expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);
            let b = btmp as f64;


            x.push(b);
            remaining_n -= b as u64;
            remaining_p -= p;
        }
    }
    x
}

// #[test]
pub fn test_multinomial(dim: i32){
    let n = 1000;
    let p = (1..dim).map(|x| x as f64).collect();
    // println!("{:?}", p);

    let x = multinomial_sample(n, p);
    let n2: f64 = x.iter().sum();
    // println!("{:?} {}", x, n2);
}
// #[test]
pub fn test_multinomial_stats(dim: i32){
    let n = 1000;
    let p: Vec<f64> = (1..dim).map(|x| x as f64).collect();
    // println!("{:?}", *p);

    let mut r = rand::thread_rng();
    let n = Multinomial::new(&p, n).unwrap();
    let x = n.sample(&mut r);
    let n2: f64 = x.iter().sum();
    // println!("{:?} {}", x, n2);

}
