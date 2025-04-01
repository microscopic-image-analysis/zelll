use crate::Angstrom;

#[derive(Clone, Copy, Debug, Default)]
pub enum Element {
    #[default]
    Carbon,
    Hydrogen,
    Oxygen,
    Nitrogen,
    Sulfur,
    Selenium,
}

impl Element {
    /// Returns the van-der-Waals radius of this `Element` variant (in `Å = 10⁻¹⁰m`).
    #[inline]
    pub fn radius(&self) -> Angstrom {
        use Element::*;
        match self {
            Carbon => 1.70,
            Hydrogen => 1.09, // 1.10, // 1.20,
            Oxygen => 1.52,
            Nitrogen => 1.55,
            Sulfur => 1.80,
            Selenium => 1.90,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Atom {
    pub(crate) element: Element,
    pub(crate) coords: [Angstrom; 3],
}

// The blanket implementation for Particle can use this
impl From<Atom> for [Angstrom; 3] {
    fn from(atom: Atom) -> Self {
        atom.coords
    }
}
