use std::f64::consts::TAU;

use real_time_fir_iir_filters::{iir::first::FirstOrderFilter, Filter};

const G_CLIP: f64 = 0.04;
const P_HARD_CLIP: f64 = 0.2;
const P_HARD_CLIP_FINAL: f64 = 0.8;
const FILTER_TUBE_CASCADE: usize = 9;
const FILTER_TUBE_FREQUENCIES: [f64; FILTER_TUBE_CASCADE] = [
    22.0,
    47.0,
    100.0,
    220.0,
    470.0,
    1000.0,
    2200.0,
    4700.0,
    10000.0
];
const G_OVERTONES: f64 = 1.0;
const P_SOFT_CLIP: f64 = 1.0;
const G_SOFT_CLIP: f64 = 1.13;
const P_SIGMOID_CLIP: f64 = 0.99;

#[derive(Clone, Copy)]
pub struct TubeStage
{
    filter_tube: [FirstOrderFilter<f64>; FILTER_TUBE_CASCADE],
}

impl TubeStage
{
    fn soft_clip(x: f64, min: f64, max: f64) -> f64
    {
        let x = x.max(min).min(max);
        x - 1.0/(max - x).exp() + 1.0/(x - min).exp()
    }
    

    pub fn new() -> Self
    {
        Self {
            filter_tube: FILTER_TUBE_FREQUENCIES.map(|f| FirstOrderFilter::new(TAU*f))
        }
    }

    pub fn next(&mut self, rate: f64, mut x: f64) -> f64
    {
        x = Self::soft_clip(x, -200.0, 140.0);
        x = Self::soft_clip(x, -20.0, 14.0)*P_HARD_CLIP_FINAL + (1.0 - P_HARD_CLIP_FINAL)*x;

        for (k, filter_tube) in self.filter_tube.iter_mut()
            .enumerate()
        {
            let [z1, z2] = filter_tube.filter(rate, x);

            x = z2*G_OVERTONES - z1*P_SOFT_CLIP;
            x = G_CLIP*x;
            x = P_SIGMOID_CLIP*((x*2.0).exp() - 1.0)/(1.0 + (x*2.0).exp()) + (1.0 - P_SIGMOID_CLIP)*x;
            x = x/G_CLIP;
            let div = 1.0/(k as f64 + 2.0).sqrt();
            x = Self::soft_clip(x, -200.0*div, 140.0*div)*P_HARD_CLIP + (1.0 - P_HARD_CLIP)*x;

            x = z1*(1.0 + P_SOFT_CLIP*G_SOFT_CLIP) + x*G_SOFT_CLIP;
        }
        
        x = G_CLIP*x;
        x = P_SIGMOID_CLIP*((x*2.0).exp() - 1.0)/(1.0 + (x*2.0).exp()) + (1.0 - P_SIGMOID_CLIP)*x;
        x = x/G_CLIP;

        x
    }

    pub fn reset(&mut self)
    {
        for filter in self.filter_tube.iter_mut()
        {
            filter.reset();
        }
    }
}