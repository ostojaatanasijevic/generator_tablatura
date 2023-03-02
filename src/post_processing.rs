use fast_float::parse;
use rayon::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
use std::error::Error;
use std::f32::consts::PI;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::string;
use std::thread;

use crate::cli::Args;
use crate::fft;
use crate::fourier::calculate_window_function;
use crate::AVG_LEN;
use crate::THREADS;
use crate::T_RES;

pub fn threaded_fft_fir_filtering(
    note_intensity: Vec<Vec<Vec<f32>>>,
    h: &Vec<f32>,
    args: &Args,
) -> Vec<Vec<Vec<f32>>> {
    let mut note_intensity = note_intensity;
    let window = crate::fourier::calculate_window_function(h.len(), &args.window_function); // blackman je bolji od hann
    let mut handles = vec![];
    for i in 0..6 {
        let h = h.clone();
        let mut string_data = note_intensity.remove(0);
        let window = window.clone();
        let dec_len = args.decemation_len;

        handles.push(thread::spawn(move || {
            let mut out_string_data = Vec::new();

            for n in 0..crate::DIFF_TABLE[i] {
                let mut note_data = string_data.remove(0);
                let temp = fft::convolve(&note_data, &h, &window, None); // applying low pass filter
                out_string_data.push(block_max_decemation(&temp, dec_len));
            }

            (i, out_string_data)
        }));
    }

    let mut note_intensity = vec![Vec::new(); 6];
    for handle in handles {
        let tmp = handle.join().unwrap();
        note_intensity[tmp.0] = tmp.1;
    }

    for wire in 1..6 {
        for i in 0..20 - crate::DIFF_TABLE[wire] {
            let temp = note_intensity[wire - 1][i].clone();
            note_intensity[wire].push(temp);
        }
    }

    note_intensity
}

pub fn eliminate_by_string(string_data: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let mut out = vec![vec![0.0; string_data[0].len()]; string_data.len()];

    for t in 0..string_data[0].len() {
        //find max in this moment
        let mut max = string_data[0][t];
        let mut max_note = 0;
        for note in 1..string_data.len() {
            if string_data[note][t] > max {
                max = string_data[note][t];
                max_note = note;
            }
        }

        out[max_note][t] = max;
    }

    out
}

//pub fn butterworth(n: usize, attenuation_db: f32) -> (Vec<f32>, Vec<f32>) {}

pub fn lp_filter(fp: f32, n: usize) -> Vec<f32> {
    let mut out = vec![0.0; 2 * n + 1];

    for t in 0..n {
        let m = n - t;
        out[t] = ((fp * m as f32).sin()) / (m as f32 * PI);
    }
    out[n] = fp / PI;
    for t in 1..n {
        out[t + n] = ((fp * t as f32).sin()) / (t as f32 * PI);
    }

    out
}

pub fn hp_filter(fp: f32, n: usize) -> Vec<f32> {
    let mut out = vec![0.0; 2 * n + 1];

    for t in 0..n - 1 {
        let m = n - t;
        out[t] = ((3.14 * m as f32).sin() - (3.14 * fp * m as f32).sin()) / (m as f32 * 3.14)
    }
    out[n] = 1.0 - fp;
    for t in 1..n {
        out[t + n] = ((3.14 * t as f32).sin() - (3.14 * fp * t as f32).sin()) / (t as f32 * 3.14)
    }

    out
}

pub fn fir_filter(h: &Vec<f32>, input: &Vec<f32>) -> Vec<f32> {
    let len = h.len();
    let mut x = vec![0.0; len];
    x.append(&mut input.clone());

    let mut out = vec![0.0; input.len() + len];

    for i in 0..input.len() {
        for ai in 0..len - 1 {
            out[i + len] += h[ai] * x[i + len - ai];
        }
    }

    out
}

pub fn iir_filter(b: &Vec<f32>, a: &Vec<f32>, input: &Vec<f32>) -> Vec<f32> {
    let mut out = vec![0.0; input.len() + b.len()];

    for i in b.len()..input.len() {
        for bi in 0..b.len() {
            out[i] -= b[bi] * out[i - bi];
        }

        for ai in 0..a.len() {
            out[i] += a[ai] * input[i - ai];
        }
    }

    out
}

pub fn process_of_elimination(data: &Vec<Vec<Vec<f32>>>) -> Vec<Vec<Vec<f32>>> {
    let mut out = data.clone();

    for t in 0..out[0][0].len() {
        //PRVIH PET NOTA E žice mogu samo da proizvede ta žica
        //tkd ovde bezbedno može da se radi prosta eliminacija
        let mut max = 0.0;
        let mut index = 0;

        for i in 0..5 {
            if out[5][i][t] > max {
                max = out[5][i][t];
                index = i;
            }
        }

        for i in 0..5 {
            if i != index {
                out[5][i][t] = 0.0;
            }
        }

        //isto važi i za poslednje note e žice
        let mut max = 0.0;
        let mut index = 0;
        for i in 15..20 {
            if out[0][i][t] > max {
                max = out[0][i][t];
                index = i;
            }
        }

        for i in 15..20 {
            if i != index {
                out[0][i][t] = 0.0;
            }
        }
    }

    out
}

pub fn binary_schmitt(data: &Vec<Vec<Vec<f32>>>, high: f32, low: f32) -> Vec<Vec<Vec<f32>>> {
    let mut out = data.clone();

    for string in 0..data.len() {
        for note in 0..data[0].len() {
            let max = data[string][note]
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap();

            let high_thresh = high;
            let low_thresh = low;

            out[string][note][0] = 0.0;
            for t in 1..data[string][note].len() {
                if out[string][note][t - 1] == 0.0 {
                    if data[string][note][t] > high_thresh {
                        out[string][note][t] = 1.0;
                    } else {
                        out[string][note][t] = 0.0;
                    }
                } else if out[string][note][t - 1] == 1.0 {
                    if data[string][note][t] < low_thresh {
                        out[string][note][t] = 0.0;
                    } else {
                        out[string][note][t] = 1.0;
                    }
                } else {
                    out[string][note][t] = 0.0;
                }
            }
        }
    }

    out
}
pub fn schmitt(data: &Vec<Vec<Vec<f32>>>, high: f32, low: f32) -> Vec<Vec<Vec<f32>>> {
    let mut out = data.clone();

    for string in 0..data.len() {
        for note in 0..data[0].len() {
            let max = data[string][note]
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap();

            let high_thresh = high;
            let low_thresh = low;

            out[string][note][0] = 0.0;
            for t in 1..data[string][note].len() {
                if out[string][note][t - 1] == 0.0 {
                    if data[string][note][t] <= high_thresh {
                        out[string][note][t] = 0.0;
                    }
                } else if out[string][note][t - 1] > 0.0 {
                    if data[string][note][t] < low_thresh {
                        out[string][note][t] = 0.0;
                    }
                } else {
                    out[string][note][t] = 0.0;
                }
            }
        }
    }

    out
}

pub fn to_dirac(data: &Vec<Vec<Vec<f32>>>) -> Vec<Vec<Vec<f32>>> {
    let mut out = vec![vec![Vec::new(); data[0].len()]; 6];

    for string in 0..data.len() {
        for note in 0..data[0].len() {
            for t in 0..data[string][note].len() - 1 {
                if data[string][note][t + 1] > data[string][note][t] {
                    out[string][note].push(1.0);
                } else {
                    out[string][note].push(0.0);
                }
            }
        }
    }

    out
}

pub fn decemation(data: &Vec<f32>, n: usize) -> Vec<f32> {
    let mut out: Vec<f32> = Vec::new();
    for t in 0..(data.len() / n) {
        out.push(data[t * n]);
    }
    out
}

pub fn block_average_decemation(data: &Vec<f32>, avg_len: usize) -> Vec<f32> {
    let mut out = Vec::new();
    for t in 0..data.len() / avg_len - 1 {
        let mut avg = data[t * avg_len..(t + 1) * avg_len].iter().sum::<f32>() / avg_len as f32;
        out.push(avg);
    }

    out
}

pub fn block_max_decemation(data: &Vec<f32>, avg_len: usize) -> Vec<f32> {
    //HORRIFICLLY SLOW CODE
    //WORKS THO
    //UPDATE TO USE SORT AND A STACK OF MAX VALUES
    let mut out = Vec::new();

    for t in 0..data.len() / avg_len {
        let max = *data[t * avg_len..(t + 1) * avg_len]
            .iter()
            .max_by(|a, b| a.total_cmp(&b))
            .unwrap();
        out.push(max);
    }

    out
}

#[derive(Debug)]
pub struct NotePeak {
    pub time: f32,
    pub ampl: f32,
    pub index: usize,
}

pub fn find_peaks(data: &Vec<f32>, threshold: f32) -> Vec<NotePeak> {
    let mut local_peaks: Vec<NotePeak> = Vec::new();
    //let threshold = data.iter().max_by(|a, b| a.total_cmp(b)).unwrap() * threshold_ratio;

    for t in 1..data.len() - 1 {
        if (data[t] > data[t - 1]) && (data[t] > data[t + 1]) {
            if (data[t] > threshold) {
                local_peaks.push(NotePeak {
                    time: t as f32 * T_RES,
                    ampl: data[t],
                    index: t,
                });
            }
        }
    }

    local_peaks
}

fn generate_volume_map(song: &Vec<i16>, decemation_len: usize) -> Vec<f32> {
    // first abs
    let mut before: Vec<f32> = Vec::new();

    for sample in song.iter() {
        before.push(sample.abs() as f32);
    }

    let out = block_max_decemation(&before, decemation_len);
    out
}

fn vector_multiply(a: &Vec<f32>, b: &Vec<f32>) -> Vec<f32> {
    let mut out: Vec<f32> = Vec::new();
    let len = std::cmp::min(a.len(), b.len());
    println!("first len {}, second {}", a.len(), b.len());

    for t in 0..len {
        out.push(a[t] * b[t]);
    }

    out
}
