use statrs::distribution::{Multinomial, Binomial};  // warning the statrs::Binomial has very slow sampling (sum of Bernullis)
use probability;  // use probability::distribution::Binomial instead, which does inverse cdf sampling
use probability::distribution::Sample;
use probability::prelude::*;
use rand;
use rand::distributions::Distribution;
use probability::source::Xorshift128Plus;


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
            // let b = Binomial::new(_ptmp, remaining_n).expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);
            let bin = Binomial::new(_ptmp, remaining_n).unwrap_or_else(|_| panic!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n));
            
            // expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);
            let b = bin.sample(&mut r);

            x.push(b);
            remaining_n -= b as u64;
            remaining_p -= p;
        }
    }
    x
}

pub fn multinomial_sample(n: u64, pvec: &Vec<f64>, source: &mut Xorshift128Plus) -> Vec<f64>{
    /*
    my own multinomial sampling, using the fact that all marginals are binomial

    statrs version does the same algorithm, but relies internally on a statrs::distribution::Binomial
    which is extremely slow.
    */
    let dim = pvec.len();
    let mut x :Vec<f64> = Vec::with_capacity(dim);

    // normalize the pvec
    let _sum: f64 = pvec.iter().sum();
    let pvec_norm: Vec<f64> = pvec.iter().map(|x| x/_sum).collect();

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
            let btmp = probability::distribution::Binomial::new(remaining_n as usize, _ptmp ).sample(source); //.expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);
            let b = btmp as f64;


            x.push(b);
            remaining_n -= b as u64;
            remaining_p -= p;
        }
    }
    x
}

// use std::slice::binary_search;

pub fn multinomial_sample_binary_search(n: u64, pvec: &Vec<f64>, source: &mut Xorshift128Plus) -> Vec<f64>{

    // initialize to zero
    let mut x: Vec<f64> = std::iter::repeat(0.0).take(pvec.len()).collect();

    // normalize the pvec
    let _sum: f64 = pvec.iter().sum();
    let pvec_norm: Vec<f64> = pvec.iter().map(|x| x/_sum).collect();

    let mut _acc = 0.0;
    let cdf: Vec<f64> = pvec_norm
        .iter()
        .map(|x| {
            _acc += x;
            _acc
        })
        .collect();

    // println!("{:?}", cdf);

    let rv_uniform = Uniform::new(0.0, 1.0);
    for _i in 0..n{
        let r = rv_uniform.sample(source);
        
        // bit weird logic here: we probably wont find the exact value, so we'll never get Result::Ok
        // but instead Result:Err will tell us the location wher to be inserted!
        // 
        let res = cdf.binary_search_by(|v| {
            v.partial_cmp(&r).unwrap()});

        let interval = match res{
            Ok(index_exact_match) => {println!("Should not happen {}", index_exact_match); index_exact_match}
            // [0.1, 0.30000000000000004, 0.6000000000000001, 1.0]
            // r=0.12 ->  index= 1 i.e. the intervale 0.1-03
            Err(first_index_larger_than_r) => first_index_larger_than_r
        };
        // println!("{} {}", r, interval);
        x[interval] += 1.0;
         
    }
    x
    // println!("{:?}", x);

}

#[test]
pub fn test_multinomial_binary(){
    use std::time::Instant;

    // Realistic sizes
    // dim = 91_355_762 records,
    // N   = 124_166_898 counts

    // dim = 50_000_000;
    // n = 100_000_000;
    // takes about 380sec currently

    let mut random_source = source::default(4);
    let dim = 50_000_000;
    let n = 100_000_000;
    let p = (1..dim).map(|x| x as f64).collect();
    // println!("{:?}", p);
    let now = Instant::now();
    let x = multinomial_sample_binary_search(n, &p, &mut random_source);
    let elapsed_time = now.elapsed();
     println!("Time for binary search N={} dim={}: {} secs", n, dim, elapsed_time.as_secs());

     let s:f64 = x.iter().fold(0.0, |acc, el| acc+*el);
     assert_eq!(n as f64, s)
}

// #[test]
pub fn test_multinomial(dim: i32){

    let mut random_source = source::default(42);

    let n = 1000;
    let p = (1..dim).map(|x| x as f64).collect();
    // println!("{:?}", p);

    let x = multinomial_sample(n, &p, &mut random_source);
    let _n2: f64 = x.iter().sum();
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
    let _n2: f64 = x.iter().sum();
    // println!("{:?} {}", x, n2);

}

// #[test]
pub fn multinomial_speed_eval(){
    use std::fs::File;
    use std::time::Instant;
    use core::ops::{Add, Div, Sub};
    use std::fmt::Debug;
    
    pub fn linspace<T>(x0: T, xend: T, n: u16) -> Vec<T>
    where
        T: Sub<Output = T> + Add<Output = T> + Div<Output = T> + Clone + Debug,
        u16: TryInto<T> + TryInto<usize>,
        <u16 as TryInto<T>>::Error: Debug,
    {
        let segments: T = (n - 1)
            .try_into()
            .expect("requested number of elements did not fit into T");
        let n_size: usize = n.try_into()
            .expect("requested number of elements exceeds usize");
            
        let dx = (xend - x0.clone()) / segments;
    
        let mut x = vec![x0; n_size];
    
        for i in 1..n_size {
            x[i] = x[i - 1].clone() + dx.clone();
        }
    
        x
    }
    
    // let dims = vec![1_000, 10_000, 100_000, 1_000_000];
    // let N: Vec<usize> = vec![1_000, 10_000, 100_000, 1_000_000];

    let dims: Vec<usize> = linspace(4.0,7.0, 20).iter().map(|x| 10_f64.powf(*x) as usize).collect();
    let N: Vec<usize> = linspace(4.0,7.0, 20).iter().map(|x| 10_f64.powf(*x) as usize).collect();

    println!("{:?}", dims);

    let mut random_source = source::default(42);

    struct Res {
        n: usize,
        d: usize,
        time: f32
    }

    let mut results:Vec<Res> = Vec::new();
    for n in N{
        for d in dims.clone(){

            println!("n={} d={}", n, d);
            let pvec: Vec<_> = (1..d).map(|x| x as f64).collect();

            let now = Instant::now();
            let _x = multinomial_sample_binary_search(n as u64, &pvec, &mut random_source);
            let elapsed_time = now.elapsed().as_secs_f32();

            let r = Res{n, d, time:elapsed_time};
            results.push(r)
            }
    }

    let mut fh = File::create("/tmp/multinom_speed.csv").unwrap();
    use std::io::Write;
    writeln!(fh, "N,d,time").unwrap();

    for r in results{
        writeln!(fh, "{},{},{}", r.n, r.d, r.time).unwrap();
    }
}
