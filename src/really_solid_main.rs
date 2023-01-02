#[allow(dead_code)]
#[allow(unused_mut)]
#[allow(unused_parens)]

mod misc;

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

struct NotesToIndex{
    note: i32,
    index: i32,
}

struct Peaks{
    freq: f32,
    ampl: f32,
}

#[derive(Debug)]
struct NotePeak{
    time: f32,
    ampl: f32,
    index: usize,
}


const SAMPLE: usize = 4096 * 2;
const THREADS: usize = 8; 
const F_RES: f32 = 2.0 * 44100.0 / SAMPLE as f32;
const T_RES: f32 = 1.0 / (44100.0 * 2.0);

const BULLSHIT: &str = "pure_sine_samples/audiocheck.net_sin_";
const BS: &str = "Hz_-3dBFS_3s.wav";
const HERZ: [&str; 43] = [ "82.41", "87.31","92.5","98.0","103.83","110.0","116.54","123.47","130.81",
"138.59","146.83","155.56","164.81","174.61","185.0","196.0","207.65","220.0","233.08","246.94","261.63",
"277.18","293.66","311.13","329.63","349.23","369.99","392.0","415.3","440.0","466.16","493.88","523.25",
"554.37","587.33","622.25","659.26","698.46","739.99","830.61","880.0","932.33","987.77"];

//--------------------------------------------------------------
//
//SPEED OVERVIEW
//za N = 8192 i THREAD 8
//trenutno treba 40s za convoluciju cele pesme.
//treba 1.63s za foirier cele pesme
//duplo ubrzanje: koristi realfft biblioteku, ili radi fft nad 2 različita signala u isto
//vremeduplo ubrzanje: koristi realfft biblioteku, ili radi fft nad 2 različita signala u isto vreme
//
//za N = 4096 i THREAD 8
//fourier traje 1.25s
//conv oko 32s
//--------------------------------------------------------------

fn main(){
    //ADD THREAD DETECTION FOR INDIVIDUAL CPUs

    //let lookup = hertz_to_notes();
    let window = calculate_window_function(SAMPLE);
    let fft_data = threaded_fourier("januar.wav", &window);    
    //let conv_data = fourier(&bru);    
    //to_tabs(conv_data, lookup);
    //graph_by_index(&fft_data, &lookup);
    let sample_ffts = calculate_sample_ffts(&window);
    //conv_with_samples(&fft_data, &sample_ffts);
    let note_intensity = threaded_conv(&fft_data, &sample_ffts);
    parse_peaks(&note_intensity);
}

fn calculate_window_function(N: usize) -> Vec<f32>{
    let mut out: Vec<f32> = vec![]; 

    for i in 0..N{
        out.push((0.54 - 0.46 * (( 2.0 * PI * (i as f32) / SAMPLE as f32) as f32).cos()));
    }

    out
}
fn parse_peaks(data: &Vec<Vec<Complex<f32>>>){
    let note_timeline = data[20].clone(); 
    //for note_timeline in data{
        let list_of_peaks = find_peaks(&note_timeline, 0.3);
        for peak in list_of_peaks{
            println!("{:?}",peak);
        }
    //}
}

fn find_peaks(data: &Vec<Complex<f32>>, threshold_ratio: f32) -> Vec<NotePeak>{
    let mut local_peaks: Vec<NotePeak> = Vec::new();
    let threshold = data.iter().max_by(|a, b| a.norm().total_cmp(&b.norm())).unwrap().norm() * threshold_ratio;

    for t in 1..data.len()-1{
        if (data[t].norm() > data[t-1].norm())
        && (data[t].norm() > data[t+1].norm()){
            if(data[t].norm() > threshold){
                local_peaks.push(NotePeak{time: t as f32 * T_RES, ampl: data[t].norm(), index: t}); 
            }
        }
    }

    local_peaks 
}


fn conv_with_samples(fft_data: &Vec<Vec<Complex<f32>>>,
                    sample_ffts: &Vec<Vec<Complex<f32>>>,
                    thread_index: &i32)
                    -> Vec<Vec<Complex<f32>>>{

    let mut final_buffer = vec![Vec::<Complex<f32>>::new();43];
    let mut planner = FftPlanner::<f32>::new();
    let ifft = planner.plan_fft_inverse(SAMPLE);
   
    for t in 0..fft_data[0].len(){
        let mut s_buffer = vec![Complex{ re: 0.0, im: 0.0}; SAMPLE];
        
        // POTENCIAL SPEED UP; AVOID FOR; USE MAP
        for note in 0..43{
        //println!("time: {}s", (j as i32 + (thread_index * data.len() as i32 / SAMPLE as i32)) as f32 * SAMPLE as f32 * T_RES);
            //let mut time = Instant::now();
            //
            s_buffer = fft_data.iter().zip(sample_ffts[note].iter())
              .map(|(x,y)| x[t]*y.conj() / (SAMPLE * SAMPLE) as f32).collect();
            //
            /*
            s_buffer = fft_data[t].iter().zip(sample_ffts[note].iter())
              .map(|(x,y)| x*y.conj() / (SAMPLE * SAMPLE) as f32).collect();
            */
            //println!("{:?}",time.elapsed());

            ifft.process(&mut s_buffer);
            final_buffer[note].extend(s_buffer.clone());
        }

        // 14.92s
        // 14.75s
        
        //println!("time: {}s",((t as i32 + thread_index * fft_data[0].len() as i32) as f32 * SAMPLE as f32 * T_RES ));
    }

    /* 
    for note in 0..43{
        let mut filename = String::new();
        filename.push_str(HERZ[note]);
        filename.push_str(add_str);
        filename.push_str("_conv.png");
        draw_plot(&filename, final_buffer[note].iter()
                  .map(|a| a.norm() as f64).collect(),
                  T_RES as f32, 1).unwrap();
    }
    */ 

    final_buffer
}

fn threaded_conv(fft_data: &Vec<Vec<Complex<f32>>>,
                    sample_ffts: &Vec<Vec<Complex<f32>>>)
                    -> Vec<Vec<Complex<f32>>>{
   
    let mut chunks_of_fft: Vec<Vec<Vec<Complex<f32>>>> = Vec::new();
    let mut chunk_lenght = fft_data[0].len() / THREADS / 60;
    println!("{chunk_lenght}");

    for c in 0..THREADS{
        // /
        let mut tmp: Vec<Vec<Complex<f32>>> = Vec::new();
        for fsc in fft_data{ // fsc je fft sample channel, ima ih SAMPLE mnogo
            tmp.push(fsc[c*chunk_lenght..(c+1)*chunk_lenght].to_vec());
        }
        chunks_of_fft.push(tmp);
        // /
        //chunks_of_fft.push(fft_data[c*chunk_lenght..(c+1)*chunk_lenght].to_vec());
    }

    let mut handles = vec![]; 
    for i in 0..THREADS{
        let gas = chunks_of_fft[i].clone();
        let sam = sample_ffts.clone(); 
        handles.push(thread::spawn(move || {
            conv_with_samples(&gas,&sam, &(i as i32))
        }));
    }

    let mut joined_data: Vec<Vec<Complex<f32>>> = vec![Vec::new();43];
    for handle in handles{ 
        let tmp = handle.join().unwrap(); 
        
        for note in 0..43{
            joined_data[note].extend(&tmp[note]);
        }
    }
     
    for note in 24..43{
        let mut filename = String::new();
        filename.push_str(HERZ[note]);
        filename.push_str("_conv.png");
        draw_plot(&filename,joined_data[note].iter()
                  .map(|a| a.norm() as f64).collect(),
                  T_RES as f32, 1).unwrap();
    }
    
    joined_data 
}

fn graph_by_index(fft_data: &Vec<Vec<Complex<f32>>>, lookup: &Vec<String>){
    let mut index: i32 = 0;
    let mut ni_list: Vec<NotesToIndex> = Vec::new();
    for i in 0..43{
        let lup = HERZ[i as usize].parse::<f32>().unwrap();
        
        let index: i32 = (lup / F_RES).round() as i32;
        println!("{} {}", index, index as f32 * F_RES);
        ni_list.push(NotesToIndex{ note: i as i32, index: index });
    }

    for i in 0..43{
        let mut filename = String::new();
        filename.push_str(HERZ[i]);
        filename.push_str(".png");
            
        draw_plot(&filename, fft_data[ni_list[i].index as usize]
                  .iter().map(|a| a.norm() as f64).collect(), T_RES, 1).unwrap();
    }
}

fn to_tabs( conv: Vec<Vec<f32>>, lookup: Vec<String>){
    for t in 0..conv[0].len()/1024{
        let mut prisutne_freq = Vec::<usize>::new();
        print!("t: {}s :: ",1024.0 * T_RES * t as f32);
        for n in  0..43{
            if conv[n][t * 1024] > 30.0{
                print!("{} -:- ", lookup[n]);
                prisutne_freq.push(n);
            }
        }
        
        print!("\n");

        // PROVERI KOJI SU DEO ISTE GRUPE S
        for p in prisutne_freq.iter(){
            print!("{}-",p);
        }
        print!("\n");
       
        loop{
            if prisutne_freq.len() < 2 {
                break;
            }
            let mut broke = 0;
            //POLUPANA LOGIKA, NIŠTA NE VALJA
            for i in 0..prisutne_freq.len()-1{
                if(prisutne_freq[i] == prisutne_freq[i+1]-1){
                    if(conv[prisutne_freq[i]][t] > conv[prisutne_freq[i+1]][t] ){
                        prisutne_freq.remove(i+1);
                        broke = 1;
                        break;
                    }else{
                        prisutne_freq.remove(i);
                        broke = 1;
                        break;
                    }
                }    
            }
            if broke == 0{
                break;
            }
        }
        

        print!("\n");
        for p in prisutne_freq.iter(){
            print!("{}-",p);
        }
        print!("\n");
        
    }

}

fn threaded_fourier(filename: &str, window: &Vec<f32>) -> Vec<Vec<Complex<f32>>>{
    let mut fft_mem: Vec<Vec<Complex<f32>>> = vec![Vec::new();SAMPLE];
    let mut fajl_pesme = File::open(Path::new(filename)).unwrap();
    let (header, pesma) = wav::read(&mut fajl_pesme).unwrap();

    if !pesma.is_sixteen(){ panic!("Wav file : {filename} isn't 16 bit! "); }
    let data = pesma.as_sixteen().unwrap();
    
    let mut chunks_of_the_song: Vec<Vec<i16>> = vec![];
    let mut chunk_lenght = data.len() / THREADS;
    println!("{chunk_lenght}");

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

    //let mut joined_data: Vec<Vec<Complex<f32>>> = vec![Vec::new();SAMPLE];
    let mut joined_data: Vec<Vec<Complex<f32>>> = vec![];
    for handle in handles{ 
        let mut tmp = handle.join().unwrap(); 
       
        joined_data.extend(tmp);
        /*
        for note in 0..SAMPLE{
            joined_data[note].extend(&tmp[note]);
        }
        */
    }

    //SHIFT AROUND
    let mut out = vec![Vec::new();SAMPLE]; 
    for t in 0..joined_data.len(){
        for i in 0..SAMPLE{
            out[i].push(joined_data[t][i]);
        }

    }
    println!("threaded fourier done!");
    //joined_data
    out
}

fn fourier(data: Vec<i16>, window: &Vec<f32>, thread_index: &i32) -> Vec<Vec<Complex<f32>>>{
    let mut fft_mem: Vec<Vec<Complex<f32>>> = vec![];
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(SAMPLE);
    
    for j in 0..data.len()/SAMPLE{
        let mut t = Instant::now();
        let mut pesma_fft = vec![Complex{ re: 0.0, im: 0.0}; SAMPLE];
        
        for i in 0..SAMPLE{
            pesma_fft[i].re = (data[i + j*SAMPLE as usize] as f32) * window[i];
        }
        fft.process(&mut pesma_fft);

      
        fft_mem.push(pesma_fft);

        println!("{:?}",t.elapsed());

        //println!("time: {}s", (j as i32 + (thread_index * data.len() as i32 / SAMPLE as i32)) as f32 * SAMPLE as f32 * T_RES);
    }
    fft_mem
}

fn calculate_sample_ffts(window: &Vec<f32>) -> Vec<Vec<Complex<f32>>>{
    let mut samples_fft = vec![Vec::<Complex<f32>>::new();43];
    
    for note in 0..HERZ.len(){
        let mut filename = String::from(BULLSHIT);
        filename.push_str(HERZ[note]);
        filename.push_str(BS);

        let mut file = File::open(Path::new(&filename[..]))
            .expect(&format!("Can't open file named {filename}"));
        let (_, raw_data) = wav::read(&mut file)
            .expect("Can't read file, im retarded");
        if !raw_data.is_sixteen(){ panic!("Wav file : {filename} isn't 16 bit! "); }

        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(SAMPLE);
        let mut fft_data = vec![Complex{ re: 0.0, im: 0.0}; SAMPLE];

        let data = raw_data.as_sixteen().unwrap();
           
        let offset: usize = 0;
        for i in 0..SAMPLE{
            fft_data[i].re = (data[i + offset * SAMPLE] as f32) * window[i];
            fft_data[i].im = 0.0;
        }

        fft.process(&mut fft_data);
        samples_fft[note] = fft_data; //PODELI SA SAMPLES
    }

    println!("samples calculated!");
    samples_fft
} 

fn draw_plot(plot_name: &str, data: Vec<f64>,
             mul: f32, div: usize) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(plot_name, (1024, 768)).into_drawing_area();
    let mut freq: Vec<f32> = Vec::new();
    let max = data.iter().max_by(|a, b| a.total_cmp(b)).unwrap();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption("Area Chart Demo", ("sans-serif", 40))
        .build_cartesian_2d(0..(data.len()/div - 1), 0.0..*max)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_label_formatter(&|x| format!("{}", (*x as f32 * mul) as f32))
        .draw()?;

    chart.draw_series(
        AreaSeries::new(
            (0..).zip(data.iter()).map(|(x, y)| (x, *y)),
            0.0,
            &RED.mix(0.2),
        )
        .border_style(&RED),
    )?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", plot_name);
    println!("max is {max}");
    Ok(())
}
