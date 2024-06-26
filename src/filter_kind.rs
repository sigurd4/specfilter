use core::fmt::Display;

use signal_processing::gen::filter::IirFilterType;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum FilterKind
{
    //Bessel,
    Butterworth,
    Chebyshev1,
    Chebyshev2,
    Elliptic
}

impl FilterKind
{
    pub const VARIANT_COUNT: usize = core::mem::variant_count::<Self>();
    pub const VARIANTS: [Self; Self::VARIANT_COUNT] = [
        Self::Butterworth,
        Self::Chebyshev1,
        Self::Chebyshev2,
        Self::Elliptic
    ];
    pub const VARIANT_NAMES: [&'static str; Self::VARIANT_COUNT] = [
        "Butterworth",
        "Chebyshev 1",
        "Chebyshev 2",
        "Elliptic"
    ];
}

impl Into<IirFilterType> for FilterKind
{
    fn into(self) -> IirFilterType
    {
        match self
        {
            FilterKind::Butterworth => IirFilterType::Butterworth,
            FilterKind::Chebyshev1 => IirFilterType::Chebyshev1,
            FilterKind::Chebyshev2 => IirFilterType::Chebyshev2,
            FilterKind::Elliptic => IirFilterType::Elliptic
        }
    }
}

impl Display for FilterKind
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        write!(f, "{}", Self::VARIANT_NAMES[*self as usize])
    }
}