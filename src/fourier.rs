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
use crate::INTERPOL;
use crate::THREADS;
use crate::STRINGS;
use crate::post_processing::single_rolling_max_decemation;

pub fn threaded_dtft_and_conv(song: &Vec<i16>,
                              sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
                              window: &Vec<f32>
                              ) -> Vec<Vec<Vec<f32>>>{
    
    let mut chunks_of_the_song: Vec<Vec<i16>> = vec![];
    let mut chunk_lenght = song.len() / THREADS;

    for c in 0..THREADS{
        chunks_of_the_song.push(song[c*chunk_lenght..(c+1)*chunk_lenght].to_vec());
    }

    let mut handles = vec![]; 
    for i in 0..THREADS{
        let gas = chunks_of_the_song[i].clone();
        let win = window.clone();
        let sam = sample_ffts.clone();
        handles.push(thread::spawn(move || {
            dtft_and_conv(&gas, &sam, &win)
        }));
    }

    let mut joined_data: Vec<Vec<Vec<f32>>> = vec![vec![Vec::new();sample_ffts[0].len()]; 6];
    for handle in handles{ 
        let tmp = handle.join().unwrap(); 
        
        for string in 0..6{
        for note in 0..sample_ffts[0].len(){
            joined_data[string][note].extend(&tmp[string][note]);
        }
        }
    }

    println!("threaded dtft and conv done!");
    joined_data
}
/*
pub fn dtft_and_conv(input_chunk: &Vec<i16>,
                     sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
                     window: &Vec<f32>
                     ) -> Vec<Vec<Vec<f32>>>{ 
  
    let mut sf = sample_ffts.clone();
  let mut final_buffer = vec![vec![Vec::<f32>::new(); sample_ffts[0][0].len()];6];
  let mut planner = FftPlanner::<f32>::new();
  let fft = planner.plan_fft_forward(NFFT*2);
  let ifft = planner.plan_fft_inverse(NFFT*2);
   
  let out = sample_ffts.clone();
  let chunk_lenght = input_chunk.len();
  let num_of_chunks = chunk_lenght / SAMPLE * INTERPOL - INTERPOL;

    for string in 0..sample_ffts.len(){
        for note in 0..sample_ffts[0].len(){
            let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; NFFT*2];
            sf[string][note].extend(pesma_fft);
            
        }
    }
             
    for string in 0..sample_ffts.len(){
        for note in 0..sample_ffts[0].len(){
            let mut carry = vec![0.0; NFFT];

            for c in 0..num_of_chunks{
                let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; NFFT*2];
      
                for i in 0..NFFT{
                    pesma_fft[i].re = (input_chunk[c*SAMPLE/INTERPOL + i % SAMPLE] as f32) * window[i%SAMPLE];
                }
  
                fft.process(&mut pesma_fft);
    
                let mut s_buffer = vec![Complex{ re: 0.0, im: 0.0}; NFFT*2];

                s_buffer = pesma_fft.iter().zip(sf[string][note].iter())
                .map(|(x,y)| x*y.conj() / (SAMPLE * SAMPLE * 1000000) as f32).collect();
     
                ifft.process(&mut s_buffer);

                let mut out: Vec<f32> = s_buffer.iter().map(|a| a.norm()).collect();
                for i in 0..NFFT{
                    out[i] += carry[i];
                    carry[i] = out[i + NFFT];
                } 
                // decemate here; save RAM
                let decemated = single_rolling_max_decemation(&out[0..NFFT].to_vec(), AVG_LEN); 

                final_buffer[string][note].extend(decemated);
            }
        }
    }

  println!("dtft and conv done!");
  final_buffer 
}
*/
//THESE TWO DONT JIVE, E2 has an erronious bump 
pub fn dtft_and_conv(input_chunk: &Vec<i16>,
                     sample_ffts: &Vec<Vec<Vec<Complex<f32>>>>,
                     window: &Vec<f32>
                     ) -> Vec<Vec<Vec<f32>>>{ 
    
  let mut final_buffer = vec![vec![Vec::<f32>::new(); sample_ffts[0][0].len()];6];
  let mut planner = FftPlanner::<f32>::new();
  let fft = planner.plan_fft_forward(NFFT);
  let ifft = planner.plan_fft_inverse(SAMPLE);
    
  let out = sample_ffts.clone();
  let chunk_lenght = input_chunk.len();
  let num_of_chunks = chunk_lenght / SAMPLE * INTERPOL - INTERPOL;

  for c in 0..num_of_chunks{
    let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; NFFT];
      
    for i in 0..NFFT{
      pesma_fft[i].re = (input_chunk[c*SAMPLE/INTERPOL + i % SAMPLE] as f32) * window[i%SAMPLE];
    }
  
    fft.process(&mut pesma_fft);
    
    for string in 0..sample_ffts.len(){
    for note in 0..sample_ffts[0].len(){
      let mut s_buffer = vec![Complex{ re: 0.0, im: 0.0}; NFFT];

      s_buffer = pesma_fft.iter().zip(sample_ffts[string][note].iter())
      .map(|(x,y)| x*y.conj() / (SAMPLE * SAMPLE * 1000000) as f32).collect();
     
      ifft.process(&mut s_buffer);

      let out: Vec<f32> = s_buffer.iter().map(|a| a.norm()).collect();
      
      // decemate here; save RAM
      let out = single_rolling_max_decemation(&out, AVG_LEN); 

      final_buffer[string][note].extend(out);
    }
  }
    }

  println!("dtft and conv done!");
  final_buffer 
}


pub fn dtft(input_chunk: &Vec<i16>,
                     window: &Vec<f32>
                     ) -> Vec<f32>{ 
    
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(NFFT);
    
    let c: usize = 0;
    let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; NFFT];
      
    for i in 0..SAMPLE{
        pesma_fft[i].re = (input_chunk[c*SAMPLE + i] as f32) * window[i];
    }
  
    fft.process(&mut pesma_fft);
    pesma_fft[0..NFFT/16].iter().map(|x| x.norm()).collect()
}

pub fn threaded_conv(fft_data: &Vec<Vec<Complex<f32>>>,
                sample_ffts: &Vec<Vec<Complex<f32>>>,
                percent: f32)
                    -> Vec<Vec<Complex<f32>>>{
   
    let mut chunks_of_fft: Vec<Vec<Vec<Complex<f32>>>> = Vec::new();
    let mut chunk_lenght = ((fft_data.len() / THREADS) as f32 * percent) as usize;
    
    for c in 0..THREADS{
        chunks_of_fft.push(fft_data[c*chunk_lenght..(c+1)*chunk_lenght].to_vec());
    }

    let mut handles = vec![]; 
    for i in 0..THREADS{
        let gas = chunks_of_fft[i].clone();
        let sam = sample_ffts.clone(); 
        handles.push(thread::spawn(move || {
            conv_with_samples(&gas,&sam, &(i as i32))
        }));
    }

    let mut joined_data: Vec<Vec<Complex<f32>>> = vec![Vec::new();sample_ffts.len()];
    for handle in handles{ 
        let tmp = handle.join().unwrap(); 
        
        for note in 0..sample_ffts.len(){
            joined_data[note].extend(&tmp[note]);
        }
    }

    println!("threaded conv done!"); 
    joined_data 
}

fn conv_with_samples(fft_data: &Vec<Vec<Complex<f32>>>,
                    sample_ffts: &Vec<Vec<Complex<f32>>>,
                    thread_index: &i32)
                    -> Vec<Vec<Complex<f32>>>{
    //FFT DATA [time][0..SAMPLE]
    let mut final_buffer = vec![Vec::<Complex<f32>>::new();sample_ffts.len()];
    let mut planner = FftPlanner::<f32>::new();
    let ifft = planner.plan_fft_inverse(SAMPLE);
   
    for t in 0..fft_data.len(){
        for note in 0..sample_ffts.len(){
            //let mut time = Instant::now();
            let mut s_buffer = vec![Complex{ re: 0.0, im: 0.0}; SAMPLE];
            
            s_buffer = fft_data[t].iter().zip(sample_ffts[note].iter())
              .map(|(x,y)| x*y.conj() / (SAMPLE * SAMPLE * 1000000) as f32).collect();

            ifft.process(&mut s_buffer);
            //MAYBE SLOW; PROLLY NOT
            final_buffer[note].extend(s_buffer);
            //println!("{:?}",time.elapsed());
        }
    }

    final_buffer
}

pub fn threaded_fourier(filename: &str, window: &Vec<f32>) -> Vec<Vec<Complex<f32>>>{
    let mut fft_mem: Vec<Vec<Complex<f32>>> = vec![Vec::new();SAMPLE];
    let mut fajl_pesme = File::open(Path::new(filename)).unwrap();
    let (header, pesma) = wav::read(&mut fajl_pesme).unwrap();

    if !pesma.is_sixteen(){ panic!("Wav file : {filename} isn't 16 bit! "); }
    let data = pesma.as_sixteen().unwrap();
    
    let mut chunks_of_the_song: Vec<Vec<i16>> = vec![];
    let mut chunk_lenght = data.len() / THREADS;

    for c in 0..THREADS{
        chunks_of_the_song.push(data[c*chunk_lenght..(c+1)*chunk_lenght].to_vec());
    }

    let mut handles = vec![]; 
    for i in 0..THREADS{
        let gas = chunks_of_the_song[i].clone();
        let win = window.clone();
        handles.push(thread::spawn(move || {
            fourier(gas, &win, &(i as i32))
        }));
    }

    let mut joined_data: Vec<Vec<Complex<f32>>> = vec![];
    for handle in handles{ 
        let mut tmp = handle.join().unwrap(); 
       
        joined_data.extend(tmp);
    }

    println!("threaded fourier done!");
    joined_data
}

fn fourier(data: Vec<i16>, window: &Vec<f32>, thread_index: &i32) -> Vec<Vec<Complex<f32>>>{
    let mut fft_mem: Vec<Vec<Complex<f32>>> = vec![];
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(SAMPLE);
    
    for j in 0..data.len()/SAMPLE{
        let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; SAMPLE];
        
        for i in 0..SAMPLE{
            pesma_fft[i].re = (data[i + j*SAMPLE as usize] as f32) * window[i];
        }
        fft.process(&mut pesma_fft);
        fft_mem.push(pesma_fft);

        //let mut t = Instant::now();
        //println!("{:?}",t.elapsed());
    }
    fft_mem
}

pub fn calculate_sample_ffts(window: &Vec<f32>) -> Vec<Vec<Vec<Complex<f32>>>>{
    let start: usize = 0;
    let stop: usize = 22;
    let mut samples_fft = vec![vec![Vec::<Complex<f32>>::new();stop - start]; 6 ];
    
    for string in 0..STRINGS.len(){
    for note in start..stop{
        let mut filename = String::from("midi/");
        filename.push_str(STRINGS[string]);
        filename.push_str("/");
        filename.push_str(STRINGS[string]);
        filename.push_str(&note.to_string());
        filename.push_str(".wav");

        let mut file = File::open(Path::new(&filename[..]))
            .expect(&format!("Can't open file named {filename}"));
        let (_, raw_data) = wav::read(&mut file)
            .expect(&format!("Can't read file, im retarded: ~{filename}~"));
        if !raw_data.is_sixteen(){ panic!("Wav file : {filename} isn't 16 bit! "); }

        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(NFFT);
        let mut fft_data = vec![Complex{ re: 0.0, im: 0.0}; NFFT];

        let data = raw_data.as_sixteen().unwrap();
           
        let offset: usize = 0;
        for i in 0..NFFT{
            fft_data[i].re = (data[i%SAMPLE] as f32) * window[i%SAMPLE];
        }

        fft.process(&mut fft_data);
        samples_fft[string][note - start] = fft_data; //PODELI SA SAMPLES
    }
    }

    println!("samples calculated!");
    samples_fft
}

pub fn calculate_window_function(n: usize) -> Vec<f32>{
    let mut out: Vec<f32> = vec![];
    let a = 0.543478261;
    //for i in 0..n{ out.push((a - (1.0-a)*(( 2.0 * PI * (i as f32) / n as f32) as f32).cos())); }
    for i in 0..n{ out.push((0.5 - 0.5*(( 2.0 * PI * (i as f32) / n as f32) as f32).cos())); }
    //for i in 0..n{ out.push(1.0); }
    out
}

pub fn open_song(filename: &str, percent: f32) -> Vec<i16>{
    let mut fajl_pesme = File::open(Path::new(filename)).unwrap();
    let (header, pesma) = wav::read(&mut fajl_pesme).unwrap();

    if !pesma.is_sixteen(){ panic!("Wav file : {filename} isn't 16 bit! "); }
    let data = pesma.as_sixteen().unwrap().to_vec();
     
    println!("song loaded!");
    let tmp_len = ((data.len() as f32 * percent) as usize);
    let len = tmp_len - tmp_len % SAMPLE;
    data[0..len].to_vec()
}


