use fast_float::parse;
use plotters::prelude::*;
use rayon::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
use std::error::Error;
use std::f32::consts::PI;
use std::fs::read_to_string;
use std::fs::File;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::thread;

use crate::offset_table;
use crate::post_processing::block_average_decemation;
use crate::post_processing::block_max_decemation;
use crate::post_processing::fir_filter;
use crate::Note;
use crate::AVG_LEN;
use crate::HERZ;
use crate::NFFT;
use crate::SAMPLE;
use crate::STRINGS;
use crate::THREADS;

const BROJ_ZICA: usize = 6;
const BROJ_PRAGOVA: usize = 20;

pub fn convolve(
    input_large: &Vec<i16>,
    input_small: &Vec<Complex<f32>>,
    window: &Vec<f32>,
    convolution_type: &str,
) -> Vec<f32> {
    let mut start_index: usize = 0;
    let mut stop_index: usize = SAMPLE;
    let mut nfft: usize = SAMPLE;

    if convolution_type == "circular" {
        start_index = 0;
        stop_index = SAMPLE;
        nfft = SAMPLE;
    } else if convolution_type == "save" {
        start_index = SAMPLE;
        stop_index = SAMPLE * 2;
        nfft = SAMPLE * 2;
    } else if convolution_type == "add" {
        start_index = 0;
        stop_index = SAMPLE * 2;
        nfft = SAMPLE * 2;
    }

    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let ifft = planner.plan_fft_inverse(nfft);

    let chunk_lenght = input_large.len();
    let num_of_chunks: usize = chunk_lenght / SAMPLE;

    let mut out = vec![0.0; input_large.len() + nfft];

    for c in 0..num_of_chunks {
        let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft];

        for i in 0..SAMPLE {
            pesma_fft[i].re = (input_large[c * SAMPLE + i] as f32) * window[i] / 65536.0;
        }

        fft.process(&mut pesma_fft);

        let mut s_buffer: Vec<Complex<f32>> = pesma_fft
            .iter()
            .zip(input_small.iter())
            .map(|(x, y)| x * y.conj())
            .collect();

        ifft.process(&mut s_buffer);

        let current: Vec<f32> = s_buffer[start_index..stop_index]
            .iter()
            .map(|a| a.norm())
            .collect();

        for i in 0..current.len() {
            out[c * SAMPLE + i] += current[i];
        }
    }

    out
}
