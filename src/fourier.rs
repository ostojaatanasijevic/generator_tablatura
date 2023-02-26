#![allow(warnings, unused)]

use fast_float::parse;
use plotters::prelude::*;
use rayon::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
use std::error::Error;
use std::f32::consts::PI;
use std::fs::File;
use std::io::IoSlice;
use std::io::Write;
use std::path::Path;
use std::thread;

use crate::offset_table;
use crate::post_processing::block_average_decemation;
use crate::post_processing::block_max_decemation;
use crate::Note;
use crate::AVG_LEN;
use crate::HERZ;
use crate::STRINGS;
use crate::THREADS;

const broj_zica: usize = 6;
const broj_pragova: usize = 20;

pub trait unsinkable {
    fn to_float(&self) -> f32;
}

impl unsinkable for i16 {
    fn to_float(&self) -> f32 {
        *self as f32
    }
}
impl unsinkable for f32 {
    fn to_float(&self) -> f32 {
        *self as f32
    }
}

fn num_to_float<T: unsinkable>(input: T) -> f32 {
    input.to_float()
}

pub fn average(data: &Vec<f32>) -> f32 {
    let mut out = 0.0;
    for bit in data {
        out += bit;
    }
    out
}

pub fn calculate_sample_note(
    note: &Note,
    window: &Vec<f32>,
    nfft: usize,
    padd: usize,
) -> Vec<Complex<f32>> {
    let data = open_sample_note(&note);

    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft + padd);

    let mut fft_data = vec![Complex { re: 0.0, im: 0.0 }; nfft + padd];
    for i in 0..nfft {
        fft_data[i].re = (data[i] as f32) * window[i] / 65536.0;
    }

    fft.process(&mut fft_data);

    //println!("{} sample note calculated!", &note.name);
    fft_data
}

//HANN FUNKCIJA MOŽE SE RAČUNATI KAO Y[k] = X[K] - 1/2(X[K-1] + X[K+1])
//ovim se izbegava operacija množenja jer se samo šiftuje i sabira
//ne pamti se ni prozorska funkcija

pub fn calculate_window_function(n: usize, wt: &str) -> Vec<f32> {
    let mut out: Vec<f32> = vec![];
    let a = 0.543478261;
    if wt == "blackman" {
        for i in 0..n {
            out.push((a - (1.0 - a) * ((2.0 * PI * (i as f32) / n as f32) as f32).cos()));
        }
    }
    if wt == "hann" {
        for i in 0..n {
            out.push((0.5 - 0.5 * ((2.0 * PI * (i as f32) / n as f32) as f32).cos()));
        }
    }
    if wt == "rect" {
        for i in 0..n {
            out.push(1.0);
        }
    }
    if wt == "kaiser" {
        for i in 0..n {
            out.push(1.0);
        }
    }
    if wt == "flattop" {
        for i in 0..n {
            out.push(1.0);
        }
    }
    out
}

pub fn open_song(filename: &str, seconds: f32) -> Vec<i16> {
    let song_samples = (seconds * 44100.0) as usize;

    let mut file =
        File::open(Path::new(filename)).expect(&format!("Can't open file named {filename}"));
    let (_, raw_data) =
        wav::read(&mut file).expect(&format!("Can't read file named {filename}, im retarded"));
    let data = raw_data
        .as_sixteen()
        .expect(&format!("Wav file : {filename} isn't 16 bit!"))[0..song_samples]
        .to_vec();

    println!("song loaded!");
    data
}

pub fn open_sample_note(note: &Note) -> Vec<i16> {
    let filename = format!("midi/{}/{}.wav", &note.name[0..1], &note.name);

    let mut file =
        File::open(Path::new(&filename[..])).expect(&format!("Can't open file named {filename}"));
    let (_, raw_data) =
        wav::read(&mut file).expect(&format!("Can't read file, im retarded: ~{filename}~"));
    let data = raw_data
        .as_sixteen()
        .expect(&format!("Wav file : {filename} isn't 16 bit!"));

    data.to_vec()
}

pub fn open_sample_notes(sample_len: usize) -> Vec<Vec<Vec<i16>>> {
    let start: usize = 0;
    let stop: usize = 20;
    let mut samples_fft = vec![vec![Vec::new(); stop - start]; 6];

    for string in 0..STRINGS.len() {
        for note in start..stop {
            //use format!
            //let filename = format!("midi/{}/{}{}.wav", STRINGS[string], STRINGS[string], &note);

            let filename = format!(
                "pure_sine_samples/audiocheck.net_sin_{}Hz_-3dBFS_3s.wav",
                HERZ[offset_table[5 - string] + note]
            );

            let mut file = File::open(Path::new(&filename))
                .expect(&format!("Can't open file named {filename}"));
            let (_, raw_data) =
                wav::read(&mut file).expect(&format!("Can't read file, im retarded: ~{filename}~"));
            let data = raw_data
                .as_sixteen()
                .expect(&format!("Waw file: {} is not 16 bit", &filename))[0..sample_len]
                .to_vec();

            samples_fft[string][note] = data;
        }
    }

    println!("samples calculated!");
    samples_fft
}
