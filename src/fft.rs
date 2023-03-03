use crate::fourier;
use crate::harmonics::Note;
use fast_float::parse;
use plotters::prelude::*;
use rayon::prelude::*;

use realfft::RealFftPlanner;
use rustfft::{num_complex::Complex, FftPlanner};

use std::error::Error;
use std::f32::consts::PI;
use std::fs::read_to_string;
use std::fs::File;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::thread;

use crate::cli::Args;
use crate::fourier::unsinkable;
use crate::post_processing::block_average_decemation;
use crate::post_processing::block_max_decemation;
use crate::post_processing::fir_filter;
use crate::AVG_LEN;
use crate::HERZ;
use crate::OFFSET_TABLE;
use crate::STRINGS;
use crate::THREADS;

const BROJ_ZICA: usize = 6;
const BROJ_PRAGOVA: usize = 20;

pub fn threaded_interlaced_convolution_realfft(
    song: &Vec<i16>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
    window: &Vec<f32>,
    args: &Args,
) -> Vec<Vec<Vec<f32>>> {
    let mut note_intensity = Vec::new();
    let mut handles = vec![];
    for i in 0..6 {
        let song = song.clone();
        let sample_notes = sample_notes[i].clone();
        let window = window.clone();
        let conv_type = args.conv_type.clone();

        handles.push(thread::spawn(move || {
            let mut notes_on_string = vec![Vec::new(); 20];
            for n in 0..(crate::DIFF_TABLE[i]) {
                notes_on_string[n] = interlaced_convolution_realfft(
                    &song,
                    &sample_notes[n],
                    &window,
                    Some(window.len()),
                    512,
                );
            }

            notes_on_string
        }));
    }

    let mut joined_data: Vec<f32> = Vec::new();
    for handle in handles {
        let tmp = handle.join().unwrap();
        note_intensity.push(tmp);
    }

    for wire in 1..6 {
        for i in 0..20 - crate::DIFF_TABLE[wire] {
            note_intensity[wire][i + crate::DIFF_TABLE[wire]] = note_intensity[wire - 1][i].clone();
        }
    }

    note_intensity
}

pub fn threaded_interlaced_convolution(
    song: &Vec<i16>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
    window: &Vec<f32>,
    args: &Args,
) -> Vec<Vec<Vec<f32>>> {
    let mut note_intensity = Vec::new();
    let mut handles = vec![];
    for i in 0..6 {
        let song = song.clone();
        let sample_notes = sample_notes[i].clone();
        let window = window.clone();
        let conv_type = args.conv_type.clone();

        handles.push(thread::spawn(move || {
            let mut notes_on_string = vec![Vec::new(); 20];
            for n in 0..(crate::DIFF_TABLE[i]) {
                notes_on_string[n] =
                    interlaced_convolution(&song, &sample_notes[n], &window, None, 512);
            }

            notes_on_string
        }));
    }

    let mut joined_data: Vec<f32> = Vec::new();
    for handle in handles {
        let tmp = handle.join().unwrap();
        note_intensity.push(tmp);
    }

    for wire in 1..6 {
        for i in 0..20 - crate::DIFF_TABLE[wire] {
            note_intensity[wire][i + crate::DIFF_TABLE[wire]] = note_intensity[wire - 1][i].clone();
        }
    }

    note_intensity
}

pub fn redundant_threaded_interlaced_convolution(
    song: &Vec<i16>,
    sample_notes: &Vec<Vec<Vec<i16>>>,
    window: &Vec<f32>,
    args: &Args,
) -> Vec<Vec<Vec<f32>>> {
    let mut note_intensity = Vec::new();
    let mut handles = vec![];
    for i in 0..6 {
        let song = song.clone();
        let sample_notes = sample_notes[i].clone();
        let window = window.clone();
        let conv_type = args.conv_type.clone();

        handles.push(thread::spawn(move || {
            let mut notes_on_string = vec![Vec::new(); 20];
            for n in 0..20 {
                notes_on_string[n] =
                    interlaced_convolution(&song, &sample_notes[n], &window, None, 512);
            }

            notes_on_string
        }));
    }

    let mut joined_data: Vec<f32> = Vec::new();
    for handle in handles {
        let tmp = handle.join().unwrap();
        note_intensity.push(tmp);
    }

    note_intensity
}

pub fn convolve<T: unsinkable, U: unsinkable>(
    input_large: &Vec<T>,
    input_small: &Vec<U>,
    window: &Vec<f32>,
    sample_len: Option<usize>,
) -> Vec<f32> {
    let s: usize;
    match sample_len {
        None => s = input_small.len(),
        Some(sample_len) => s = sample_len,
    }

    //check power of 2
    let sample_len = s;
    let next_pow = 2.0_f32.powf((s as f32).log2().ceil()) as usize;

    let mut nfft: usize = next_pow * 2;

    let mut planner = RealFftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let ifft = planner.plan_fft_inverse(nfft);

    //PROCESS SMALL CHUNK
    let mut h = vec![Complex { re: 0.0, im: 0.0 }; nfft / 2 + 1];
    let mut h_r = vec![0.0; nfft];

    let window = fourier::calculate_window_function(input_small.len(), "blackman");
    for i in 0..sample_len {
        h_r[i] = (input_small[i].to_float()) * window[i] / 65536.0;
    }

    fft.process(&mut h_r, &mut h);

    let chunk_lenght = input_large.len();
    let num_of_chunks: usize = chunk_lenght / sample_len;

    let mut out = vec![0.0; input_large.len() + nfft];

    for c in 0..num_of_chunks {
        let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft / 2 + 1];
        let mut pesma_fft_r = vec![0.0; nfft];

        for i in 0..sample_len {
            pesma_fft_r[i] = (input_large[c * sample_len + i].to_float()) * window[i] / 65536.0;
        }

        fft.process(&mut pesma_fft_r, &mut pesma_fft);

        let mut s_buffer: Vec<Complex<f32>> = pesma_fft
            .iter()
            .zip(h.iter())
            .map(|(x, y)| x * y.conj())
            .collect();

        let mut current = vec![0.0; nfft];
        ifft.process(&mut s_buffer, &mut current);

        for i in 0..current.len() {
            out[c * sample_len + i] += current[i];
        }
    }

    out
}

pub fn interlaced_convolution_realfft<T: unsinkable, U: unsinkable>(
    input_large: &Vec<T>,
    input_small: &Vec<U>,
    window: &Vec<f32>,
    sample_len: Option<usize>,
    inter_size: usize,
) -> Vec<f32> {
    //check power of 2
    let next_pow = sample_len.unwrap();

    let sample_len = input_small.len();
    let mut start_index: usize = 0;
    let mut stop_index: usize = next_pow;
    let mut nfft: usize = next_pow;

    let mut realplanner = RealFftPlanner::<f32>::new();
    let fft = realplanner.plan_fft_forward(nfft);
    let ifft = realplanner.plan_fft_inverse(nfft);

    //PROCESS SMALL CHUNK
    let mut h = vec![Complex { re: 0.0, im: 0.0 }; nfft / 2 + 1];
    let mut h_r = vec![0.0; nfft];

    //potrebno samo ako se radi koherentno odabiranje, što neću raditi hahahahaahah
    //let window = fourier::calculate_window_function(sample_len, "blackman");
    for i in 0..sample_len {
        h_r[i] = (input_small[i].to_float()) * window[i] / 65536.0;
    }

    fft.process(&mut h_r, &mut h);

    let chunk_lenght = input_large.len();
    let num_of_chunks: usize = chunk_lenght / inter_size - nfft / inter_size;

    let mut out: Vec<f32> = Vec::new();

    for c in 0..num_of_chunks {
        let mut pesma_fft_r = vec![0.0; nfft];
        let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft / 2 + 1];

        for i in 0..sample_len {
            pesma_fft_r[i] = (input_large[c * inter_size + i].to_float()) * window[i] / 65536.0;
        }

        fft.process(&mut pesma_fft_r, &mut pesma_fft);

        let mut s_buffer: Vec<Complex<f32>> = pesma_fft
            .iter()
            .zip(h.iter())
            .map(|(x, y)| x * y.conj() / 10000000.0)
            .collect();

        let mut current = vec![0.0; nfft];
        ifft.process(&mut s_buffer, &mut current);
        out.extend(block_max_decemation(&current, current.len()));
    }
    out
}

pub fn interlaced_convolution<T: unsinkable, U: unsinkable>(
    input_large: &Vec<T>,
    input_small: &Vec<U>,
    window: &Vec<f32>,
    sample_len: Option<usize>,
    inter_size: usize,
) -> Vec<f32> {
    let s: usize;
    match sample_len {
        None => s = input_small.len(),
        Some(sample_len) => s = sample_len,
    }
    let sample_len = s;

    let mut start_index: usize = 0;
    let mut stop_index: usize = sample_len;
    let mut nfft: usize = sample_len;

    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    let ifft = planner.plan_fft_inverse(nfft);

    //PROCESS SMALL CHUNK
    let mut h = vec![Complex { re: 0.0, im: 0.0 }; nfft];

    for i in 0..sample_len {
        h[i].re = (input_small[i].to_float()) * window[i] / 65536.0;
    }

    fft.process(&mut h);

    let chunk_lenght = input_large.len();
    let num_of_chunks: usize = chunk_lenght / inter_size - sample_len / inter_size;

    let mut out: Vec<f32> = Vec::new();

    for c in 0..num_of_chunks {
        let mut pesma_fft = vec![Complex { re: 0.0, im: 0.0 }; nfft];

        for i in 0..sample_len {
            pesma_fft[i].re = (input_large[c * inter_size + i].to_float()) * window[i] / 65536.0;
        }

        fft.process(&mut pesma_fft);

        let mut s_buffer: Vec<Complex<f32>> = pesma_fft
            .iter()
            .zip(h.iter())
            .map(|(x, y)| x * y.conj())
            .collect();

        ifft.process(&mut s_buffer);

        let current: Vec<f32> = s_buffer[start_index..stop_index]
            .iter()
            .map(|a| a.norm()) //real abs is the same as norm
            .collect();
        out.extend(current);
    }

    out
}
