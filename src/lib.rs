#![feature(adt_const_params)]
#![feature(split_array)]
#![feature(variant_count)]
#![feature(array_chunks)]
#![feature(array_methods)]
#![feature(generic_const_exprs)]

use std::f32::EPSILON;
use std::f32::consts::TAU;
use std::process::Command;
use std::sync::Arc;
use std::sync::atomic::{Ordering, AtomicU8};

use array_math::{ArrayOps, SliceMath};
use filter_type::FilterType;
use num_traits::float::TotalOrder;
use num_traits::{Float, FloatErrorKind, ParseFloatError, Zero};
use parameters::{SpecfilterParam, SpecfilterParamData};
use signal_processing::analysis::FiltOrd;
use signal_processing::gen::filter::{Butter, Cheby1, Cheby2, Ellip, FilterGenError, FilterGenPlane, IirDesign};
use signal_processing::operations::filtering::FilterMut;
use signal_processing::systems::{Rtf, Sos, Tf, Zpk};
use signal_processing::transforms::filter::Stabilize;
use signal_processing::transforms::system::ToSos;
use signal_processing::Plane;
use tube_stage::TubeStage;
use vst::{prelude::*, plugin_main};

use crate::filter_kind::FilterKind;

use self::parameters::{SpecfilterParameters};

pub mod parameters;
pub mod filter_type;
pub mod filter_kind;
pub mod tube_stage;

const CHANGE: f32 = 2000.0;
const MAX_ORDER: usize = 64;

struct SpecfilterPlugin
{
    pub param: Arc<SpecfilterParameters>,
    param_prev: Option<SpecfilterParamData>,
    filter_type: FilterType,
    filter: [Rtf<f64, Sos<f64, [f64; 3], [f64; 3], Vec<Tf<f64, [f64; 3], [f64; 3]>>>>; CHANNEL_COUNT],
    tubes: [TubeStage; CHANNEL_COUNT],
    rate: f64,
    host: HostCallback
}

const CHANNEL_COUNT: usize = 2;

#[test]
fn test()
{
    let mut plugin = SpecfilterPlugin::new(HostCallback::default());

    plugin.param.frequencies[0].set(1.18);
    plugin.param.frequencies[1].set(4242.703);

    plugin.generate_filter(usize::MAX).unwrap();

    println!("Ok!")
}

impl SpecfilterPlugin
{
    fn generate_filter(&mut self, buf_len: usize) -> Result<(), Box<dyn std::error::Error>>
    {
        self.param.rate.set(self.rate as f32);

        let param_next: SpecfilterParamData = (&*self.param).into();
        if let Some(param_prev) = &mut self.param_prev
        {
            if &param_next == param_prev
            {
                return Ok(())
            }
            if param_prev.filter_kind != param_next.filter_kind
            {
                for rtf in self.filter.iter_mut()
                {
                    rtf.w.clear()
                }
            }
            param_prev.change(param_next, 1.0 - (-CHANGE*buf_len as f32/(self.filter[0].sys.filtord() + 1) as f32/self.rate as f32).exp());
        }
        else
        {
            self.param_prev = Some(param_next)
        }
        let param_data = self.param_prev.as_ref()
            .unwrap();

        let mut freq = param_data.frequencies(self.rate as f32);
        let rp = param_data.passband_ripple as f64;
        let rs = param_data.stopband_attenuation as f64;

        if let Ok(Ok(freq)) = &mut freq
        {
            freq.0.sort_by(TotalOrder::total_cmp);
            freq.1.sort_by(TotalOrder::total_cmp);
        }

        let filter = match freq
        {
            Ok(freq) => match match self.param.filter_kind()
            {
                FilterKind::Butterworth => {
                    match freq
                    {
                        Ok((fp, fs)) => {
                            let (n, wp, ws, t) = signal_processing::gen::filter::buttord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::butter(n, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        },
                        Err((fp, fs)) => {
                            let (n, wp, ws, t) = signal_processing::gen::filter::buttord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::butter(n, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        }
                    }
                },
                FilterKind::Chebyshev1 => {
                    match freq
                    {
                        Ok((fp, fs)) => {
                            let (n, wp, ws, rp, t) = signal_processing::gen::filter::cheb1ord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::cheby1(n, rp, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        },
                        Err((fp, fs)) => {
                            let (n, wp, ws, rp, t) = signal_processing::gen::filter::cheb1ord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::cheby1(n, rp, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        }
                    }
                },
                FilterKind::Chebyshev2 => {
                    match freq
                    {
                        Ok((fp, fs)) => {
                            let (n, wp, ws, rs, t) = signal_processing::gen::filter::cheb2ord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::cheby2(n, rs, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        },
                        Err((fp, fs)) => {
                            let (n, wp, ws, rs, t) = signal_processing::gen::filter::cheb2ord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::cheby2(n, rs, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        }
                    }
                },
                FilterKind::Elliptic => {
                    match freq
                    {
                        Ok((fp, fs)) => {
                            let (n, wp, ws, rp, rs, t) = signal_processing::gen::filter::ellipord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::ellip(n, rp, rs, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        },
                        Err((fp, fs)) => {
                            let (n, wp, ws, rp, rs, t) = signal_processing::gen::filter::ellipord(
                                fp,
                                fs,
                                rp,
                                rs,
                                FilterGenPlane::Z { sampling_frequency: Some(self.rate) }
                            )?;
                            let n = n.min(MAX_ORDER);
                            let w = wp.comap(ws, |wp, ws| (wp + ws)*0.5);
                            Zpk::ellip(n, rp, rs, w, t, FilterGenPlane::Z { sampling_frequency: None })
                        }
                    }
                },
            }
            {
                Ok(h) => h.stabilize(Plane::Z).to_sos((), ()),
                Err(err) => match err
                {
                    FilterGenError::ZeroOrder => Sos::one(),
                    _ => Err(err)?
                }
            },
            Err(pass) => if pass
            {
                Sos::one()
            }
            else
            {
                Sos::new(vec![Tf::new([0.0; 3], [0.0, 0.0, 1.0])])
            }
        };

        /*if filter.is_stable(0.01, Plane::Z)
        {
            self.filter = filter;
        }*/

        if filter.sos.iter()
            .any(|sos| sos.a.iter().any(|a| !a.is_finite()) || sos.b.iter().any(|b| !b.is_finite()))
            || filter.sos.iter().any(|sos| sos.a.is_zero())
        {
            "a".parse::<f64>()?;
            return Ok(())
        }

        let gain = (filter.sos.iter()
                .map(|sos| sos.b.trim_zeros_front()
                        .iter()
                        .copied()
                        .map(|b| b*b)
                        .sum::<f64>()
                    /sos.a.trim_zeros_front()
                        .iter()
                        .copied()
                        .map(|a| a*a)
                        .sum::<f64>()
                ).product::<f64>()
            /self.filter[0].sys.sos.iter()
                .map(|sos| sos.b.trim_zeros_front()
                        .iter()
                        .copied()
                        .map(|b| b*b)
                        .sum::<f64>()
                    /sos.a.trim_zeros_front()
                        .iter()
                        .copied()
                        .map(|a| a*a)
                        .sum::<f64>()
                ).product::<f64>()).sqrt();
        let filter_type = param_data.filter_type(self.rate as f32);
        if gain.is_finite() && filter.filtord() == self.filter[0].sys.filtord() && filter_type != FilterType::NoPass && !(filter_type == FilterType::BandPass && self.filter_type == FilterType::BandStop)
        {
            for rtf in self.filter.iter_mut()
            {
                for w in rtf.w.iter_mut()
                {
                    *w *= gain
                }
            }
        }
        else
        {
            for rtf in self.filter.iter_mut()
            {
                rtf.w.clear()
            }
        }

        for rtf in self.filter[1..].iter_mut()
        {
            rtf.sys = filter.clone()
        }
        self.filter[0].sys = filter;
        self.filter_type = filter_type;

        Ok(())
    }

    fn process<T>(&mut self, buffer: &mut AudioBuffer<T>)
    where
        T: Float
    {
        let buf_len = buffer.samples();

        let prev_data = self.param_prev;
        if let Err(error) = self.generate_filter(buf_len)
        {
            self.param_prev = prev_data;
            if let Some(prev_data) = prev_data
            {
                self.param.reset_to(prev_data);
                if let Some(dispatcher) = self.host.raw_callback()
                {
                    // Calls dispatch in the same way that setParameterAutomated does
                    // see:
                    // opcode audioMasterAutomate is given in https://github.com/R-Tur/VST_SDK_2.4/blob/master/pluginterfaces/vst2.x/aeffect.h
                    // dispatch is given in setParameterAutomated in https://github.com/R-Tur/VST_SDK_2.4/blob/master/public.sdk/source/vst2.x/audioeffect.cpp
                    for param_id in [SpecfilterParam::Frequency1, SpecfilterParam::Frequency2, SpecfilterParam::Bandwidth1, SpecfilterParam::Bandwidth2]
                    {
                        dispatcher(self.host.raw_effect(), 0, param_id as i32, 0, std::ptr::null_mut(), self.param.get_parameter(param_id as i32));
                    }
                }
            }
            let estr = format!("{}", error);
            Command::new("cmd").args(
                &["start", "cmd.exe", "@cmd" , "/k", "echo", &estr]
            ).spawn()
            .expect("Failed to start cmd");
        }
        
        let mix = self.param.mix.get() as f64;

        for (((input_channel, output_channel), filter), tube) in buffer.zip()
            .zip(self.filter.iter_mut())
            .zip(self.tubes.iter_mut())
        {
            let x: Vec<_> = input_channel.into_iter()
                .map(|&x| x.to_f64().unwrap())
                .collect();

            let mut z = filter.filter_mut(x.as_slice());
            
            if z.iter()
                .any(|&z| z.is_nan())
            {
                filter.w.clear();
                z.fill(0.0);
                
                if let Some(prev_data) = prev_data
                {
                    self.param.reset_to(prev_data);
                    if let Some(dispatcher) = self.host.raw_callback()
                    {
                        // Calls dispatch in the same way that setParameterAutomated does
                        // see:
                        // opcode audioMasterAutomate is given in https://github.com/R-Tur/VST_SDK_2.4/blob/master/pluginterfaces/vst2.x/aeffect.h
                        // dispatch is given in setParameterAutomated in https://github.com/R-Tur/VST_SDK_2.4/blob/master/public.sdk/source/vst2.x/audioeffect.cpp
                        for param_id in [SpecfilterParam::Frequency1, SpecfilterParam::Frequency2, SpecfilterParam::Bandwidth1, SpecfilterParam::Bandwidth2]
                        {
                            dispatcher(self.host.raw_effect(), 0, param_id as i32, 0, std::ptr::null_mut(), self.param.get_parameter(param_id as i32));
                        }
                    }
                }
                self.param_prev = None;
            }

            for ((y, z), x) in output_channel.into_iter()
                .zip(z)
                .zip(x)
            {
                let mut z = tube.next(self.rate, z);
                if z.is_nan()
                {
                    z = 0.0;
                    tube.reset();
                }
                let z = z.max(-10.0).min(10.0);
                *y = T::from(z*mix + x*(1.0 - mix)).unwrap();
            }
        }
    }
}

impl Plugin for SpecfilterPlugin
{
    fn new(host: HostCallback) -> Self
    where
        Self: Sized
    {        
        SpecfilterPlugin {
            param: Arc::new(SpecfilterParameters::default()),
            param_prev: None,
            filter_type: FilterType::AllPass,
            filter: core::array::from_fn(|_| Rtf::new(Sos::one(), ())),
            tubes: core::array::from_fn(|_| TubeStage::new()),
            rate: 44100.0,
            host
        }
    }

    fn get_info(&self) -> Info
    {
        Info {
            name: "Specfilter".to_string(),
            vendor: "Soma FX".to_string(),
            presets: 0,
            parameters: SpecfilterParam::VARIANT_COUNT as i32,
            inputs: CHANNEL_COUNT as i32,
            outputs: CHANNEL_COUNT as i32,
            midi_inputs: 0,
            midi_outputs: 0,
            unique_id: 235925,
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
        self.rate = rate as f64;
    }

    fn process(&mut self, buffer: &mut AudioBuffer<f32>)
    {
        self.process(buffer)
    }
    
    fn process_f64(&mut self, buffer: &mut AudioBuffer<f64>)
    {
        self.process(buffer)
    }

    fn get_parameter_object(&mut self) -> Arc<dyn PluginParameters>
    {
        self.param.clone()
    }
}

plugin_main!(SpecfilterPlugin);