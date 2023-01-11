use std::io::Write;
use std::io::IoSlice;
use std::thread;
use std::path::Path;
use std::error::Error;
use rustfft::{FftPlanner, num_complex::Complex};
use fast_float::parse;
use rayon::prelude::*;

use crate::SAMPLE;
use crate::AVG_LEN;
use crate::THREADS;
use crate::T_RES;
use crate::F_RES;

use crate::NotePeak;
use crate::NotesToIndex;
use crate::Peaks;

pub fn process_of_elimination(data: &Vec<Vec<Vec<f32>>>) -> Vec<Vec<Vec<f32>>>{
    let mut out = data.clone();

    for t in 0..out[0][0].len(){
        //PRVIH PET NOTA E žice mogu samo da proizvede ta žica
        //tkd ovde bezbedno može da se radi prosta eliminacija
        let mut max = 0.0;
        let mut index = 0;

        for i in 0..5{
            if out[5][i][t] > max{
                max = out[5][i][t];
                index = i;
            }
        }

        for i in 0..5{
            if i != index{
                out[5][i][t] = 0.0;
            }
        }

        //isto važi i za poslednje note e žice
        let mut max = 0.0;
        let mut index = 0;
        for i in 15..20{
            if out[0][i][t] > max{
                max = out[0][i][t];
                index = i;
            }
        }

        for i in 15..20{
            if i != index{
                out[0][i][t] = 0.0;
            }
        }
    }

    out
}

pub fn schmitt(data: &Vec<Vec<Vec<f32>>>,
               high: f32,
               low: f32) -> Vec<Vec<Vec<f32>>>{

    let mut out = data.clone();

    for string in 0..data.len(){
        for note in 0..data[0].len(){
            let max = data[string][note].iter().max_by(|a, b| a.total_cmp(b)).unwrap();

            let high_thresh = high;
            let low_thresh = low;

            out[string][note][0] = 0.0;
            for t in 1..data[string][note].len(){
                if out[string][note][t-1] == 0.0{
                    if data[string][note][t] > high_thresh{
                        out[string][note][t] = 1.0;
                    }else{
                        out[string][note][t] = 0.0;
                    }
                }else if out[string][note][t-1] == 1.0{ 
                    if data[string][note][t] < low_thresh{
                        out[string][note][t] = 0.0;
                    }else{
                        out[string][note][t] = 1.0;
                    }
                }else{
                    out[string][note][t] = 0.0;
                }
            }
        }
    }

        out
}

pub fn to_dirac(data: &Vec<Vec<Vec<f32>>>) -> Vec<Vec<Vec<f32>>>{
        let mut out = vec![vec![Vec::new();data[0].len()];6];
    
    for string in 0..data.len(){
    for note in 0..data[0].len(){
        for t in 0..data[string][note].len()-1{
            if data[string][note][t+1] > data[string][note][t]{
                out[string][note].push(1.0);
            }else{
                out[string][note].push(0.0);
            }
        }
    }
    }

    out
}

pub fn block_average_decemation(data: &Vec<f32>, avg_len: usize) -> Vec<f32>{
    let mut out = Vec::new();
        for t in 0..data.len()/avg_len - 1{
            let mut avg = data[t * avg_len..(t + 1) * avg_len].iter().sum::<f32>() / avg_len as f32;
            out.push(avg);
        }

    out
}

pub fn block_max_decemation(data: &Vec<f32>, avg_len: usize) -> Vec<f32>{
    //HORRIFICLLY SLOW CODE
    //WORKS THO
    //UPDATE TO USE SORT AND A STACK OF MAX VALUES
    let mut out = Vec::new();
        
        for t in 0..data.len()/avg_len{
            let max = *data[t * avg_len..(t + 1) * avg_len].iter().max_by(|a, b| a.total_cmp(&b)).unwrap();
            out.push(max);
        }

    out
}

pub fn find_peaks(data: &Vec<f32>, threshold_ratio: f32) -> Vec<NotePeak>{
    let mut local_peaks: Vec<NotePeak> = Vec::new();
    let threshold = data.iter().max_by(|a, b| a.total_cmp(b)).unwrap() * threshold_ratio;

    for t in 1..data.len()-1{
        if (data[t] > data[t-1]) && (data[t] > data[t+1]){
            if(data[t] > threshold){
                local_peaks.push(NotePeak{time: t as f32 * T_RES, ampl: data[t], index: t}); 
            }
        }
    }

    local_peaks 
}


