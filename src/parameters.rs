use core::f32::EPSILON;
use core::sync::atomic::AtomicBool;
use std::sync::atomic::{AtomicU8, Ordering};

use num_traits::float::TotalOrder;
use num_traits::{Float, Zero};
use vst::prelude::PluginParameters;
use vst::util::AtomicFloat;

use crate::filter_type::FilterType;
use crate::filter_kind::FilterKind;

const MIN_FREQ: f32 = 1.0;
const MAX_FREQ: f32 = 30000.0;
const MIN_TRANSITION_BAND: f32 = 0.01;
const MIN_RIPPLE: f32 = 1.0;
const MAX_RIPPLE: f32 = 100.0;
const BW_EPS: f32 = 0.00001;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum SpecfilterParam
{
    FilterKind,
    Mix,
    PassbandRipple,
    StopbandAttenuation,
    Frequency1,
    Frequency2,
    Bandwidth1,
    Bandwidth2,
}

impl SpecfilterParam
{
    pub const VARIANT_COUNT: usize = core::mem::variant_count::<Self>();
    pub const VARIANTS: [Self; Self::VARIANT_COUNT] = [
        Self::FilterKind,
        Self::Mix,
        Self::PassbandRipple,
        Self::StopbandAttenuation,
        Self::Frequency1,
        Self::Frequency2,
        Self::Bandwidth1,
        Self::Bandwidth2,
    ];
}

pub struct SpecfilterParameters
{
    pub filter_kind: AtomicU8,
    pub passband_ripple: AtomicFloat,
    pub stopband_attenuation: AtomicFloat,
    pub frequencies: [AtomicFloat; 2],
    pub bandwidths: [AtomicFloat; 2],
    pub mix: AtomicFloat,
    pub rate: AtomicFloat
}

impl From<&SpecfilterParameters> for SpecfilterParamData
{
    fn from(param: &SpecfilterParameters) -> Self
    {
        Self {
            filter_kind: param.filter_kind(),
            passband_ripple: param.passband_ripple.get(),
            stopband_attenuation: param.stopband_attenuation.get(),
            frequencies: param.frequencies.each_ref()
                .map(|f| f.get()),
            bandwidths: param.bandwidths.each_ref()
                .map(|f| f.get())
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
pub struct SpecfilterParamData
{
    pub filter_kind: FilterKind,
    pub passband_ripple: f32,
    pub stopband_attenuation: f32,
    pub frequencies: [f32; 2],
    pub bandwidths: [f32; 2]
}

impl SpecfilterParamData
{
    pub fn change(&mut self, new: SpecfilterParamData, change: f32)
    {
        self.filter_kind = new.filter_kind;
        self.passband_ripple = new.passband_ripple*change + self.passband_ripple*(1.0 - change);
        self.stopband_attenuation = new.stopband_attenuation*change + self.stopband_attenuation*(1.0 - change);
        for (f1, f2) in self.frequencies.iter_mut()
            .zip(new.frequencies)
        {
            *f1 = f2*change + *f1*(1.0 - change)
        }
        for (f1, f2) in self.bandwidths.iter_mut()
            .zip(new.bandwidths)
        {
            *f1 = f2*change + *f1*(1.0 - change)
        }
    }
}

impl SpecfilterParamData
{
    pub fn frequency_data(&self, rate: f32) -> ([f32; 4], bool, bool, bool)
    {
        let mut freq = self.frequencies;
        let mut bw = self.bandwidths;

        let max_freq = MAX_FREQ.min(rate/2.0);

        let stop = freq[0] >= freq[1];
        let mut nolb = freq[0].min(freq[1]) <= MIN_FREQ + EPSILON;
        let mut noub = freq[0].max(freq[1]) >= max_freq - EPSILON;

        if stop
        {
            freq.reverse();
            bw.reverse();
        }
        let bws = ((bw[0] >= 0.0) ^ stop, (bw[1] >= 0.0) ^ stop);
        for bw in bw.iter_mut()
        {
            *bw = MIN_TRANSITION_BAND + (*bw).abs()*(1.0 - MIN_TRANSITION_BAND)
        }
        let freq = loop {
            let freq = match (nolb, noub)
            {
                (false, false) => match bws
                {
                    (false, false) => {
                        let a = (freq[0].log2()*(1.0 - bw[0]) + bw[0]*freq[1].log2()).exp2();
                        let b = (freq[1].log2()*(1.0 - bw[1]) + bw[1]*freq[0].log2()).exp2();
                        [freq[0], a.min(b), b.max(a), freq[1]]
                    },
                    (true, false) => {
                        [(freq[0].log2()*(1.0 - bw[0]) + MIN_FREQ.log2()*bw[0]).exp2(), freq[0], (freq[1].log2()*(1.0 - bw[1]) + freq[0].log2()*bw[1]).exp2(), freq[1]]
                    },
                    (false, true) => {
                        [freq[0], (freq[0].log2()*(1.0 - bw[0]) + freq[1].log2()*bw[0]).exp2(), freq[1], (freq[1].log2()*(1.0 - bw[1]) + max_freq.log2()*bw[1]).exp2()]
                    }
                    (true, true) => {
                        [(freq[0].log2()*(1.0 - bw[0]) + MIN_FREQ.log2()*bw[0]).exp2(), freq[0], freq[1], (freq[1].log2()*(1.0 - bw[1]) + max_freq.log2()*bw[1]).exp2()]
                    }
                },
                (true, false) => match bws.1
                {
                    false => {
                        [MIN_FREQ, MIN_FREQ, (freq[1].log2()*(1.0 - bw[1]) + MIN_FREQ.log2()*bw[1]).exp2(), freq[1]]
                    },
                    true => {
                        [MIN_FREQ, MIN_FREQ, freq[1], (freq[1].log2()*(1.0 - bw[1]) + max_freq.log2()*bw[1]).exp2()]
                    }
                },
                (false, true) => match bws.0
                {
                    false => {
                        [freq[0], (freq[0].log2()*(1.0 - bw[0]) + max_freq.log2()*bw[0]).exp2(), max_freq, max_freq]
                    },
                    true => {
                        [(freq[0].log2()*(1.0 - bw[0]) + max_freq.log2()*bw[0]).exp2(), freq[0], max_freq, max_freq]
                    }
                },
                (true, true) => [MIN_FREQ, MIN_FREQ, max_freq, max_freq]
            };
            let mut changed = false;
            if !nolb && !(freq[1]/freq[0] > 1.0 + BW_EPS)
            {
                changed = true;
                nolb = true
            }
            if !noub && !(freq[3]/freq[2] > 1.0 + BW_EPS)
            {
                changed = true;
                noub = true
            }
            if !changed
            {
                break freq
            }
        };

        (freq, stop, nolb, noub)
    }

    pub fn frequencies<T>(&self, rate: f32) -> Result<Result<([T; 2], [T; 2]), ([T; 1], [T; 1])>, bool>
    where
        T: Float
    {
        let (freq, stop, nolb, noub) = self.frequency_data(rate);

        if nolb && noub
        {
            return Err(!stop)
        }
        if nolb
        {
            return Ok(Err(if stop
                {
                    ([T::from(freq[3]).unwrap()], [T::from(freq[2]).unwrap()])
                }
                else
                {
                    ([T::from(freq[2]).unwrap()], [T::from(freq[3]).unwrap()])
                }))
        }
        if noub
        {
            return Ok(Err(if stop
                {
                    ([T::from(freq[0]).unwrap()], [T::from(freq[1]).unwrap()])
                }
                else
                {
                    ([T::from(freq[1]).unwrap()], [T::from(freq[0]).unwrap()])
                }))
        }
        Ok(Ok(if stop
            {
                ([T::from(freq[0]).unwrap(), T::from(freq[3]).unwrap()], [T::from(freq[1]).unwrap(), T::from(freq[2]).unwrap()])
            }
            else
            {
                ([T::from(freq[1]).unwrap(), T::from(freq[2]).unwrap()], [T::from(freq[0]).unwrap(), T::from(freq[3]).unwrap()])
            }))
    }

    pub fn filter_type(&self, rate: f32) -> FilterType
    {
        let (_, stop, nolb, noub) = self.frequency_data(rate);

        if nolb && noub
        {
            return if stop
            {
                FilterType::NoPass
            }
            else
            {
                FilterType::AllPass
            }
        }
        if nolb
        {
            return if stop
            {
                FilterType::HighPass
            }
            else
            {
                FilterType::LowPass
            }
        }
        if noub
        {
            return if stop
            {
                FilterType::LowPass
            }
            else
            {
                FilterType::HighPass
            }
        }
        if stop
        {
            FilterType::BandStop
        }
        else
        {
            FilterType::BandPass
        }
    }
}

impl Default for SpecfilterParameters
{
    fn default() -> Self
    {
        let rate = 44100.0;
        let max_freq = MAX_FREQ.min(rate/2.0);

        Self {
            filter_kind: AtomicU8::new(FilterKind::Butterworth as u8),
            passband_ripple: AtomicFloat::new(3.0),
            stopband_attenuation: AtomicFloat::new(40.0),
            frequencies: [0.3, 0.7].map(|w| AtomicFloat::new((w*(max_freq.log2() - MIN_FREQ.log2()) + MIN_FREQ.log2()).exp2())),
            bandwidths: [0.5, 0.5].map(|w| AtomicFloat::new(w)),
            mix: AtomicFloat::new(1.0),
            rate: AtomicFloat::new(rate)
        }
    }
}

impl SpecfilterParameters
{
    pub fn reset_to(&self, to: SpecfilterParamData)
    {
        self.filter_kind.store(to.filter_kind as u8, Ordering::Relaxed);
        self.passband_ripple.set(to.passband_ripple);
        self.stopband_attenuation.set(to.stopband_attenuation);
        for (f1, f2) in self.frequencies.iter()
            .zip(to.frequencies)
        {
            f1.set(f2)
        }
        for (f1, f2) in self.bandwidths.iter()
            .zip(to.bandwidths)
        {
            f1.set(f2)
        }
    }

    pub fn filter_kind(&self) -> FilterKind
    {
        FilterKind::VARIANTS[self.filter_kind.load(Ordering::Relaxed) as usize]
    }

    pub fn frequency_data(&self) -> ([f32; 4], bool, bool, bool)
    {
        SpecfilterParamData::from(self)
            .frequency_data(self.rate.get())
    }

    pub fn bandwidths(&self) -> [f32; 2]
    {
        let (freq, stop, _, _) = self.frequency_data();
        
        if stop
        {
            [freq[3] - freq[2], freq[1] - freq[0]]
        }
        else
        {
            [freq[1] - freq[0], freq[3] - freq[2]]
        }
    }

    pub fn filter_type(&self) -> FilterType
    {
        SpecfilterParamData::from(self)
            .filter_type(self.rate.get())
    }
}

impl PluginParameters for SpecfilterParameters
{
    fn get_parameter_label(&self, index: i32) -> String
    {
        match SpecfilterParam::VARIANTS[index as usize]
        {
            SpecfilterParam::FilterKind => format!("{}", self.filter_type()),
            SpecfilterParam::PassbandRipple => "dB".to_string(),
            SpecfilterParam::StopbandAttenuation => "dB".to_string(),
            SpecfilterParam::Mix => "%".to_string(),
            SpecfilterParam::Frequency1 => "Hz".to_string(),
            SpecfilterParam::Frequency2 => "Hz".to_string(),
            SpecfilterParam::Bandwidth1 => "Hz".to_string(),
            SpecfilterParam::Bandwidth2 => "Hz".to_string(),
        }
    }

    fn get_parameter_text(&self, index: i32) -> String
    {
        match SpecfilterParam::VARIANTS[index as usize]
        {
            SpecfilterParam::FilterKind => format!("{}", self.filter_kind()),
            SpecfilterParam::PassbandRipple => format!("{:.3}", self.passband_ripple.get()),
            SpecfilterParam::StopbandAttenuation => format!("{:.3}", self.stopband_attenuation.get()),
            SpecfilterParam::Mix => format!("{:.3}", 100.0*self.mix.get()),
            SpecfilterParam::Frequency1 => format!("{:.3}", self.frequencies[0].get()),
            SpecfilterParam::Frequency2 => format!("{:.3}", self.frequencies[1].get()),
            SpecfilterParam::Bandwidth1 => format!("{:.3}", self.bandwidths()[0]),
            SpecfilterParam::Bandwidth2 => format!("{:.3}", self.bandwidths()[1]),
        }
    }

    fn get_parameter_name(&self, index: i32) -> String
    {
        match SpecfilterParam::VARIANTS[index as usize]
        {
            SpecfilterParam::FilterKind => "Filter type".to_string(),
            SpecfilterParam::PassbandRipple => "Passband ripple".to_string(),
            SpecfilterParam::StopbandAttenuation => "Stopband attenuation".to_string(),
            SpecfilterParam::Mix => "Mix".to_string(),
            SpecfilterParam::Frequency1 => "Frequency 1".to_string(),
            SpecfilterParam::Frequency2 => "Frequency 2".to_string(),
            SpecfilterParam::Bandwidth1 => "Bandwidth 1".to_string(),
            SpecfilterParam::Bandwidth2 => "Bandwidth 2".to_string(),
        }
    }

    /// Get the value of parameter at `index`. Should be value between 0.0 and 1.0.
    fn get_parameter(&self, index: i32) -> f32
    {
        match SpecfilterParam::VARIANTS[index as usize]
        {
            SpecfilterParam::FilterKind => self.filter_kind.load(Ordering::Relaxed) as f32/(FilterKind::VARIANT_COUNT - 1) as f32,
            SpecfilterParam::PassbandRipple => (self.passband_ripple.get() - MIN_RIPPLE)/(MAX_RIPPLE - MIN_RIPPLE),
            SpecfilterParam::StopbandAttenuation => (self.stopband_attenuation.get() - MIN_RIPPLE)/(MAX_RIPPLE - MIN_RIPPLE),
            SpecfilterParam::Mix => self.mix.get(),
            SpecfilterParam::Frequency1 => {
                let max_freq = MAX_FREQ.min(self.rate.get()/2.0);
                (self.frequencies[0].get().log2() - MIN_FREQ.log2())/(max_freq.log2() - MIN_FREQ.log2())
            },
            SpecfilterParam::Frequency2 => {
                let max_freq = MAX_FREQ.min(self.rate.get()/2.0);
                (self.frequencies[1].get().log2() - MIN_FREQ.log2())/(max_freq.log2() - MIN_FREQ.log2())
            },
            SpecfilterParam::Bandwidth1 => (self.bandwidths[0].get() + 1.0)*0.5,
            SpecfilterParam::Bandwidth2 => (self.bandwidths[1].get() + 1.0)*0.5,
        }.min(1.0).max(0.0)
    }
    
    fn set_parameter(&self, index: i32, value: f32)
    {
        match SpecfilterParam::VARIANTS[index as usize]
        {
            SpecfilterParam::FilterKind => self.filter_kind.store((value*(FilterKind::VARIANT_COUNT - 1) as f32).round() as u8, Ordering::Relaxed),
            SpecfilterParam::PassbandRipple => self.passband_ripple.set(value*(MAX_RIPPLE - MIN_RIPPLE) + MIN_RIPPLE),
            SpecfilterParam::StopbandAttenuation => self.stopband_attenuation.set(value*(MAX_RIPPLE - MIN_RIPPLE) + MIN_RIPPLE),
            SpecfilterParam::Mix => self.mix.set(value),
            SpecfilterParam::Frequency1 => {
                let max_freq = MAX_FREQ.min(self.rate.get()/2.0);
                self.frequencies[0].set((value*1.000001*(max_freq.log2() - MIN_FREQ.log2()) + MIN_FREQ.log2()).exp2().min(max_freq).max(MIN_FREQ))
            },
            SpecfilterParam::Frequency2 => {
                let max_freq = MAX_FREQ.min(self.rate.get()/2.0);
                self.frequencies[1].set((value*1.000001*(max_freq.log2() - MIN_FREQ.log2()) + MIN_FREQ.log2()).exp2().min(max_freq).max(MIN_FREQ))
            },
            SpecfilterParam::Bandwidth1 => self.bandwidths[0].set(value*2.0 - 1.0),
            SpecfilterParam::Bandwidth2 => self.bandwidths[1].set(value*2.0 - 1.0),
        }
    }

    fn change_preset(&self, preset: i32) {}

    fn get_preset_num(&self) -> i32
    {
        0
    }

    fn set_preset_name(&self, name: String) {}

    fn get_preset_name(&self, preset: i32) -> String
    {
        "".to_string()
    }

    fn can_be_automated(&self, index: i32) -> bool
    {
        index < SpecfilterParam::VARIANT_COUNT as i32
    }

    fn get_preset_data(&self) -> Vec<u8>
    {
        SpecfilterParam::VARIANTS.map(|v| self.get_parameter(v as i32).to_le_bytes())
            .concat()
    }

    fn get_bank_data(&self) -> Vec<u8>
    {
        SpecfilterParam::VARIANTS.map(|v| self.get_parameter(v as i32).to_le_bytes())
            .concat()
    }

    fn load_preset_data(&self, data: &[u8])
    {
        for (v, &b) in SpecfilterParam::VARIANTS.into_iter()
            .zip(data.array_chunks())
        {
            self.set_parameter(v as i32, f32::from_le_bytes(b));
        }
    }

    fn load_bank_data(&self, data: &[u8])
    {
        for (v, &b) in SpecfilterParam::VARIANTS.into_iter()
            .zip(data.array_chunks())
        {
            self.set_parameter(v as i32, f32::from_le_bytes(b));
        }
    }
}