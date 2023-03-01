use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
///Generator tablatura
pub struct Args {
    #[arg(short, long, default_value_t = 3072)]
    ///Lenght of fft sample size
    pub nfft: usize,

    #[arg(short, long, default_value_t = 0.2)]
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

    #[arg(short, long, default_value = "songs/januar.wav")]
    ///Song file path
    pub file_name: String,

    #[arg(short, long, default_value_t = 20.0)]
    ///Number of seconds to analyze
    pub sec_to_run: f32,

    #[arg(short, long, default_value = "blackman")]
    ///Window function
    pub window_function: String,

    #[arg(short, long, default_value_t = 0.5)]
    ///Attenuation factor for harmonics
    pub attenuation_factor: f32,

    #[arg(short, long, default_value_t = 3.5)]
    ///Power for the harmonics ratio to be raised to
    pub power_of_harmonics: f32,

    #[arg(short, long, default_value_t = 0.1)]
    ///Outout filter cutoff frequency
    pub output_cutoff: f32,

    #[arg(short, long, default_value_t = 0.0)]
    ///Max plot value
    pub max_plot_value: f32,

    #[arg(short, long, default_value_t = 0.0)]
    ///Seek_offset_in_seconds
    pub seek_offset: f32,
}
