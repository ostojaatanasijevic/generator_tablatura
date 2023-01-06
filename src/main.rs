#![allow(warnings, unused)]

mod misc;
mod fourier;
mod post_processing;
mod plot;

use std::io::Write;
use std::io::IoSlice;
use std::time::Instant;
use std::fs::File;
use std::thread;
use std::path::Path;
use std::error::Error;
use std::f32::consts::PI;
use rustfft::{FftPlanner, num_complex::Complex};
use plotters::prelude::*;
use fast_float::parse;
use misc::hertz_to_notes;

use rayon::prelude::*;

pub struct NotesToIndex{
    note: i32,
    index: i32,
}

pub struct Peaks{
    freq: f32,
    ampl: f32,
}

#[derive(Debug)]
pub struct NotePeak{
    time: f32,
    ampl: f32,
    index: usize,
}

pub const THREADS: usize = 8; 
pub const INTERPOL: usize = 2; 
pub const SAMPLE: usize = 8192;
pub const F_RES: f32 = 2.0 * 44100.0 / SAMPLE as f32;
pub const T_RES: f32 = 1.0 / (44100.0 * 2.0);
pub const NFFT: usize = SAMPLE;
pub const AVG_LEN: usize = NFFT / 2; // mora da može da deli NFFT, da ne bi cureli podaci
pub const STRINGS: [&str; 6] = ["e","B","G","D","A","E"];

fn main(){
    //ADD THREAD DETECTION FOR INDIVIDUAL CPUs

    //let lookup = hertz_to_notes();
    let window = fourier::calculate_window_function(SAMPLE, "hann");
    let sample_ffts = fourier::calculate_sample_ffts(&window);

    let song_data = fourier::open_song("songs/jesen_stize_dunjo_moja.wav", 0.3);
    //MOŽDA TI NI NE TREBA KONVOLICIJA; MOŽDA JE DTFT DOVOLJAN ZA POSTIZANJE PRECIZNE PROCENE
    //FREKVENCIJE
    /* 
    let note_intensity = fourier::dtft(&song_data, &window);
    plot::single_plot_data_norm(note_intensity); 
    */ 
    
    // UZMI U OBZIR KAKO PROZORSKA FUNKCIJA UTIČE NA INTENZITET SIGNALA
    let note_intensity = fourier::threaded_dtft_and_conv(&song_data, &sample_ffts, &window, "save");
    /* 
    let schmitt_out = post_processing::schmitt(&note_intensity, 8.0, 4.0);
    let diracs = post_processing::to_dirac(&schmitt_out); // rising edge detection
    
    for t in 0..diracs[0].len(){
        for string in 0..diracs.len(){
        for note in 0..diracs[0].len(){
            if diracs[string][note][t] == 1.0{
                println!("{} {}", t as f32 * T_RES * AVG_LEN as f32, misc::index_to_note(note));
            }
        }
        }
    }
    */
    // AKO SU NA ISTOJ ŽICI 2 NOTE onda PROCES eliminacije
    plot::plot_data_norm(note_intensity); 
    
    //parse_peaks(&note_intensity);
}
