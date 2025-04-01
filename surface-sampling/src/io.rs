//! PDB input parsing and serialization of oriented pointcloud surfaces.
use crate::{
    Angstrom,
    atom::{Atom, Element},
};
use pdbtbx::{Atom as PdbAtom, Element as PdbElement};

impl TryFrom<PdbElement> for Element {
    type Error = &'static str;

    fn try_from(value: PdbElement) -> Result<Self, Self::Error> {
        match value {
            PdbElement::C => Ok(Element::Carbon),
            PdbElement::H => Ok(Element::Hydrogen),
            PdbElement::O => Ok(Element::Oxygen),
            PdbElement::N => Ok(Element::Nitrogen),
            PdbElement::S => Ok(Element::Sulfur),
            PdbElement::Se => Ok(Element::Selenium),
            _ => Err("This kind of element is not supported"),
        }
    }
}

impl TryFrom<&PdbAtom> for Atom {
    type Error = &'static str;

    fn try_from(value: &PdbAtom) -> Result<Self, Self::Error> {
        if let Some(elem) = value.element() {
            let element = Element::try_from(*elem)?;

            Ok(Atom {
                element,
                // TODO: `nuts-rs` wants f64 anyway and this would be cleaner
                // coords: value.pos().into(),
                coords: [
                    value.x() as Angstrom,
                    value.y() as Angstrom,
                    value.z() as Angstrom,
                ],
            })
        } else {
            Err("This atom has no associated element")
        }
    }
}

pub struct PointCloud {
    pub points: Vec<Atom>,
}

impl PointCloud {
    pub fn from_pdb_atoms<'a>(atoms: impl Iterator<Item = &'a PdbAtom>) -> PointCloud {
        PointCloud {
            points: atoms.filter_map(|a| Atom::try_from(a).ok()).collect(),
        }
    }
}
