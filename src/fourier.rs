#![allow(warnings, unused)]

use std::io::Write;
use std::io::IoSlice;
use std::fs::File;
use std::path::Path;
use std::error::Error;
use rustfft::{FftPlanner, num_complex::Complex};
use plotters::prelude::*;
use fast_float::parse;
use std::thread;
use std::f32::consts::PI;
use rayon::prelude::*;

use crate::SAMPLE;
use crate::AVG_LEN;
use crate::NFFT;
use crate::THREADS;
use crate::STRINGS;
use crate::Note;
use crate::post_processing::block_average_decemation;
use crate::post_processing::block_max_decemation;

const broj_zica: usize = 6;
const broj_pragova: usize = 20;

pub fn threaded_note_convolution(song: &Vec<i16>,
                        note: &Note,
                        window_type: &str,
                        convolution_type: &str) -> Vec<f32>{
    
    let mut next_power_of_two: usize = SAMPLE;
    let mut chunks_of_the_song: Vec<Vec<i16>> = vec![];
    let mut chunk_lenght = (song.len() / THREADS) / next_power_of_two;
    chunk_lenght = chunk_lenght * next_power_of_two;

    for c in 0..THREADS{
        chunks_of_the_song.push(song[c*chunk_lenght..(c+1)*chunk_lenght].to_vec());
    }

    let mut handles = vec![]; 
    for i in 0..THREADS{
        let song_chunk = chunks_of_the_song[i].clone();
        let n = note.clone();
        let win = String::from(window_type);
        let ct = String::from(convolution_type);
        handles.push(thread::spawn(move || {
            note_convolution(&song_chunk, &n, &win, &ct)
        }));
    }

    let mut joined_data: Vec<f32> = Vec::new();
    for handle in handles{ 
        let tmp = handle.join().unwrap(); 
        
        joined_data.extend(&tmp);
    }

    joined_data
}

pub fn note_convolution(song: &Vec<i16>,
                        note: &Note,
                        window_type: &str,
                        convolution_type: &str) -> Vec<f32>{
     
    let mut next_power_of_two: usize = SAMPLE;
    // calculate nfft based on the note
    let frequency = note.freq;
    let mut sample_len = 44100.0 / frequency;
    // get as close to 2^n as possible
    let mult_factor = (next_power_of_two as f32 / sample_len) as i32;
    sample_len = sample_len * mult_factor as f32;
    let sample_len = sample_len as usize;
    let padding_len = next_power_of_two - sample_len;
   
    let window = calculate_window_function(sample_len , window_type);
    let note_sample = calculate_sample_note(&note, &window, sample_len, padding_len);

    let start_index = 0;
    let stop_index = sample_len;

    let mut final_buffer = Vec::<f32>::new();
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(next_power_of_two);
    let ifft = planner.plan_fft_inverse(next_power_of_two);
     
    let chunk_lenght = song.len();
    let num_of_chunks: usize = chunk_lenght / sample_len;
    let mut output = vec![0.0; song.len() + sample_len * 2]; 

    for c in 0..num_of_chunks{
        let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; next_power_of_two];
      
        for i in 0..sample_len{
            pesma_fft[i].re = (song[c*sample_len + i] as f32) * window[i] / 65536.0;
        }
        
        fft.process(&mut pesma_fft);
    
        let mut s_buffer = vec![Complex{ re: 0.0, im: 0.0}; next_power_of_two];
        s_buffer = pesma_fft.iter().zip(note_sample.iter())
                    .map(|(x,y)| x*y.conj()).collect();
     
        ifft.process(&mut s_buffer);

        let current: Vec<f32> = s_buffer[start_index..stop_index].iter().map(|a| a.norm()).collect();
        for t in 0..current.len(){
            output[c*sample_len + t] += current[t];
        } 
        
        // decemate here; save RAM
        let decemeted= block_max_decemation(&output[c*sample_len..(c+1)*sample_len].to_vec(), sample_len / 4); 
        final_buffer.extend(decemeted);
    }

  //println!("dtft and conv done!");
  final_buffer 
}

pub trait unsinkable{
    fn to_float(&self) -> f32;
}

impl unsinkable for i16{
    fn to_float(&self) -> f32 { 
        *self as f32
    }
}

fn num_to_float<T: unsinkable>(input: T) -> f32{
    input.to_float()
}

pub fn dtft<T: unsinkable>(input_chunk: &Vec<T>,
            window: &Vec<f32>,
            nfft: usize
            ) -> Vec<f32>{ 
    
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft);
    
    let c: usize = 0;
    let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; nfft];
      
    for i in 0..SAMPLE{
        pesma_fft[i].re = (input_chunk[c*SAMPLE + i].to_float()) * window[i];
    }
  
    fft.process(&mut pesma_fft);
    pesma_fft.iter().map(|x| x.norm()).collect()
}

pub fn open_sample_note(note: &Note) -> Vec<i16> {
    let mut filename = String::from("midi/");
    filename.push_str(&note.name[0..1]);
    filename.push_str("/");
    filename.push_str(&note.name);
    filename.push_str(".wav");

    let mut file = File::open(Path::new(&filename[..]))
        .expect(&format!("Can't open file named {filename}"));
    let (_, raw_data) = wav::read(&mut file)
        .expect(&format!("Can't read file, im retarded: ~{filename}~"));
    let data = raw_data.as_sixteen()
        .expect(&format!("Wav file : {filename} isn't 16 bit!"));

    data.to_vec() 
}

pub fn cross_corr_notes(first: &Note,
                        second: &Note,
                        window: &Vec<f32>,
                        nfft: usize,
                        padd: usize) -> Vec<f32>{

    let first_data = open_sample_note(&first);
    let second_data = open_sample_note(&second);

    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft + padd);
    let ifft = planner.plan_fft_inverse(nfft + padd);

    let mut first_fft_data = vec![Complex{ re: 0.0, im: 0.0};nfft + padd];
    let mut second_fft_data = vec![Complex{ re: 0.0, im: 0.0};nfft + padd];
    for i in 0..nfft{
        first_fft_data[i].re = (first_data[i] as f32) * window[i] / 65536.0;
        second_fft_data[i].re = (second_data[i] as f32) * window[i] / 65536.0;
    }

    fft.process(&mut first_fft_data);
    fft.process(&mut second_fft_data);
 
    let mut s_buffer = vec![Complex{ re: 0.0, im: 0.0}; nfft];
    s_buffer = first_fft_data.iter().zip(second_fft_data.iter())
                .map(|(x,y)| x*y.conj()).collect();
     
    ifft.process(&mut s_buffer);
    let mut out: Vec<f32> = s_buffer.iter().map(|a| a.norm()).collect();

    out 
}

pub fn average(data: &Vec<f32>) -> f32{
    let mut out = 0.0; 
    for bit in data{
        out += bit;
    }
    out
}

pub fn calculate_sample_note(note: &Note,
                             window: &Vec<f32>,
                             nfft: usize,
                             padd: usize) -> Vec<Complex<f32>>{
    
    let data = open_sample_note(&note); 
        
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(nfft + padd);

    let mut fft_data = vec![Complex{ re: 0.0, im: 0.0};nfft + padd];
    for i in 0..nfft{
        fft_data[i].re = (data[i] as f32) * window[i] / 65536.0;
    }

    fft.process(&mut fft_data);

    //println!("{} sample note calculated!", &note.name);
    fft_data 
}

//HANN FUNKCIJA MOŽE SE RAČUNATI KAO Y[k] = X[K] - 1/2(X[K-1] + X[K+1])
//ovim se izbegava operacija množenja jer se samo šiftuje i sabira
//ne pamti se ni prozorska funkcija

pub fn calculate_window_function(n: usize, wt: &str) -> Vec<f32>{
    let mut out: Vec<f32> = vec![];
    let a = 0.543478261;
    if wt == "blackman" {for i in 0..n{ out.push((a - (1.0-a)*(( 2.0 * PI * (i as f32) / n as f32) as f32).cos())); }}
    if wt == "hann" {for i in 0..n{ out.push((0.5 - 0.5*(( 2.0 * PI * (i as f32) / n as f32) as f32).cos())); }}
    if wt == "rect" {for i in 0..n{ out.push(1.0); }}
    if wt == "kaiser" {for i in 0..n{ out.push(1.0);}}
    if wt == "flattop"{for i in 0..n{ out.push(1.0);}}
    out
}

pub fn open_song(filename: &str, seconds: f32) -> Vec<i16>{
    let song_samples = (seconds * 44100.0) as usize;
 
    let mut file = File::open(Path::new(filename))
            .expect(&format!("Can't open file named {filename}"));
    let (_, raw_data) = wav::read(&mut file)
            .expect(&format!("Can't read file named {filename}, im retarded"));
    let data = raw_data.as_sixteen()
            .expect(&format!("Wav file : {filename} isn't 16 bit!"))[0..song_samples].to_vec();

    println!("song loaded!");
    data
}


