use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
///Generator tablatura
pub struct Args {
    #[arg(short, long, default_value_t = 3072)]
    ///Lenght of fft sample size
    pub nfft: usize,

    #[arg(short, long, default_value_t = 0.1)]
    ///Filter cutoff frequency
    pub w: f32,

    #[arg(short, long, default_value_t = 100)]
    ///Lenght of fir filter
    pub lenght_fir: usize,

    #[arg(short, long, default_value_t = 1000)]
    ///Number of samples that get averaged into one
    pub decemation_len: usize,

    ///Circular convolution type: add, save
    #[arg(short, long, default_value = "add")]
    pub conv_type: String,

    #[arg(short, long)]
    ///Song file path
    pub file_name: String,

    #[arg(short, long, default_value_t = 10.0)]
    ///Number of seconds to analyze
    pub sec_to_run: f32,

    #[arg(short, long, default_value = "blackman")]
    ///Window function
    pub window_function: String,

    #[arg(short, long, default_value_t = 0.5)]
    ///Attenuation factor for harmonics
    pub attenuation_factor: f32,

    #[arg(short, long, default_value_t = 1.0)]
    ///Power for the harmonics ratio to be raised to
    pub power_of_harmonics: f32,

    #[arg(short, long, default_value_t = 0.1)]
    ///Outout filter cutoff frequency
    pub output_cutoff: f32,
}
