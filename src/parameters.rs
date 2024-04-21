use std::sync::atomic::{AtomicU8, Ordering};

use vst::prelude::PluginParameters;
use vst::util::AtomicFloat;

const MIN_FREQ: f32 = 20.0;
const MAX_FREQ: f32 = 20000.0;
const MIN_RES: f32 = 0.0;
const RES_CURVE: f32 = 4.0;
const MAX_RES: f32 = 20.0;

pub struct BasicFilterParameters
{
    pub filter: AtomicU8,
    pub frequency: AtomicFloat,
    pub resonance: AtomicFloat,
    pub mix: AtomicFloat
}

impl PluginParameters for BasicFilterParameters
{
    fn get_parameter_label(&self, index: i32) -> String
    {
        match index
        {
            0 => "".to_string(),
            1 => "Hz".to_string(),
            2 => "".to_string(),
            3 => "%".to_string(),
            _ => "".to_string()
        }
    }

    fn get_parameter_text(&self, index: i32) -> String
    {
        match index
        {
            0 => match self.filter.load(Ordering::Relaxed)
            {
                0 => "Low-pass".to_string(),
                1 => "Notch".to_string(),
                2 => "High-pass".to_string(),
                _ => "ERR".to_string()
            }
            1 => format!("{:.3}", self.frequency.get()),
            2 => format!("{:.3}", self.resonance.get()),
            3 => format!("{:.3}", 100.0*self.mix.get()),
            _ => "".to_string()
        }
    }

    fn get_parameter_name(&self, index: i32) -> String
    {
        match index
        {
            0 => "Filter".to_string(),
            1 => "Frequency".to_string(),
            2 => "Resonance".to_string(),
            3 => "Mix".to_string(),
            _ => "".to_string()
        }
    }

    /// Get the value of parameter at `index`. Should be value between 0.0 and 1.0.
    fn get_parameter(&self, index: i32) -> f32
    {
        match index
        {
            0 => (self.filter.load(Ordering::Relaxed) as f32)/2.0,
            1 => (self.frequency.get().log2() - MIN_FREQ.log2())/(MAX_FREQ.log2() - MIN_FREQ.log2()),
            2 => ((self.resonance.get() - MIN_RES)/(MAX_RES - MIN_RES)).powf(1.0/RES_CURVE),
            3 => self.mix.get(),
            _ => 0.0
        }
    }
    
    fn set_parameter(&self, index: i32, value: f32)
    {
        match index
        {
            0 => self.filter.store((value*2.0).round() as u8, Ordering::Relaxed),
            1 => self.frequency.set((value*(MAX_FREQ.log2() - MIN_FREQ.log2()) + MIN_FREQ.log2()).exp2()),
            2 => self.resonance.set(MIN_RES + (MAX_RES - MIN_RES)*value.powf(RES_CURVE)),
            3 => self.mix.set(value),
            _ => ()
        }
    }

    fn change_preset(&self, preset: i32) {}

    fn get_preset_num(&self) -> i32 {
        0
    }

    fn set_preset_name(&self, name: String) {}

    fn get_preset_name(&self, preset: i32) -> String {
        "".to_string()
    }

    fn can_be_automated(&self, index: i32) -> bool {
        index < 4
    }

    fn get_preset_data(&self) -> Vec<u8> {
        [
            vec![self.filter.load(Ordering::Relaxed)],
            self.frequency.get().to_le_bytes().to_vec(),
            self.resonance.get().to_le_bytes().to_vec(),
            self.mix.get().to_le_bytes().to_vec()
        ].concat()
    }

    fn get_bank_data(&self) -> Vec<u8> {
        [
            vec![self.filter.load(Ordering::Relaxed)],
            self.frequency.get().to_le_bytes().to_vec(),
            self.resonance.get().to_le_bytes().to_vec(),
            self.mix.get().to_le_bytes().to_vec()
        ].concat()
    }

    fn load_preset_data(&self, data: &[u8])
    {
        self.filter.store(data[0], Ordering::Relaxed);
        self.frequency.set(f32::from_le_bytes(*data[(data.len() - 12)..].split_array_ref().0));
        self.resonance.set(f32::from_le_bytes(*data[(data.len() - 8)..].split_array_ref().0));
        self.mix.set(f32::from_le_bytes(*data[(data.len() - 4)..].split_array_ref().0));
    }

    fn load_bank_data(&self, data: &[u8])
    {
        self.filter.store(data[0], Ordering::Relaxed);
        self.frequency.set(f32::from_le_bytes(*data[(data.len() - 12)..].split_array_ref().0));
        self.resonance.set(f32::from_le_bytes(*data[(data.len() - 8)..].split_array_ref().0));
        self.mix.set(f32::from_le_bytes(*data[(data.len() - 4)..].split_array_ref().0));
    }
}