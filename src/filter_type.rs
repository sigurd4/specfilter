use core::fmt::Display;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum FilterType
{
    LowPass,
    HighPass,
    BandPass,
    BandStop,
    NoPass,
    AllPass
}

impl FilterType
{
    pub const VARIANT_COUNT: usize = core::mem::variant_count::<Self>();
    pub const VARIANTS: [Self; Self::VARIANT_COUNT] = [
        Self::LowPass,
        Self::HighPass,
        Self::BandPass,
        Self::BandStop,
        Self::NoPass,
        Self::AllPass
    ];
    pub const VARIANT_NAMES: [&'static str; Self::VARIANT_COUNT] = [
        "Low-pass",
        "High-pass",
        "Band-pass",
        "Band-stop",
        "No-pass",
        "All-pass"
    ];
}

impl Display for FilterType
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        write!(f, "{}", Self::VARIANT_NAMES[*self as usize])
    }
}