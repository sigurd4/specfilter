#![feature(adt_const_params)]
#![feature(split_array)]

use std::f32::EPSILON;
use std::f32::consts::TAU;
use std::sync::Arc;
use std::sync::atomic::{Ordering, AtomicU8};

use real_time_fir_iir_filters::iir::{SecondOrderFilter, IIRFilter};
use vst::{prelude::*, plugin_main};

use self::parameters::{BasicFilterParameters};

pub mod parameters;

const CHANGE: f32 = 0.2;

struct BasicFilterPlugin
{
    pub param: Arc<BasicFilterParameters>,
    filter: [SecondOrderFilter; CHANNEL_COUNT],
    rate: f32
}

const CHANNEL_COUNT: usize = 2;

impl BasicFilterPlugin
{

}

impl Plugin for BasicFilterPlugin
{
    fn new(_host: HostCallback) -> Self
    where
        Self: Sized
    {
        BasicFilterPlugin {
            param: Arc::new(BasicFilterParameters {
                filter: AtomicU8::from(0),
                frequency: AtomicFloat::from(880.0),
                resonance: AtomicFloat::from(0.5f32.sqrt()),
                mix: AtomicFloat::from(1.0)
            }),
            filter: array_init::array_init(|_| SecondOrderFilter::new(TAU*880.0, 1.0)),
            rate: 44100.0
        }
    }

    fn get_info(&self) -> Info
    {
        Info {
            name: "Basic Filter".to_string(),
            vendor: "Soma FX".to_string(),
            presets: 0,
            parameters: 4,
            inputs: CHANNEL_COUNT as i32,
            outputs: CHANNEL_COUNT as i32,
            midi_inputs: 0,
            midi_outputs: 0,
            unique_id: 6436354,
            version: 1,
            category: Category::Effect,
            initial_delay: 0,
            preset_chunks: false,
            f64_precision: true,
            silent_when_stopped: true,
            ..Default::default()
        }
    }

    fn set_sample_rate(&mut self, rate: f32)
    {
        self.rate = rate;
    }

    fn process(&mut self, buffer: &mut AudioBuffer<f32>)
    {
        let filter_type = self.param.filter.load(Ordering::Relaxed);
        let mix = self.param.mix.get();

        for (c, (input_channel, output_channel)) in buffer.zip().enumerate()
        {
            self.filter[c].omega = CHANGE*self.param.frequency.get()*TAU + (1.0 - CHANGE)*self.filter[c].omega;
            self.filter[c].zeta = CHANGE*0.5/(self.param.resonance.get() + EPSILON) + (1.0 - CHANGE)*self.filter[c].zeta;

            for (input_sample, output_sample) in input_channel.into_iter()
                .zip(output_channel.into_iter())
            {
                let x = *input_sample;
                let y = self.filter[c].filter(self.rate, x)[filter_type as usize];
                *output_sample = y*mix + x*(1.0 - mix);
            }
        }
    }

    fn get_parameter_object(&mut self) -> Arc<dyn PluginParameters>
    {
        self.param.clone()
    }
}

plugin_main!(BasicFilterPlugin);