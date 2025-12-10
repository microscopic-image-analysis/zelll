use ::zelll::*;
use bincode::serde::{decode_from_slice, encode_to_vec};
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::types::PyIterator;

// TODO: While having a borrowing iterator is nice, it's expensive on Python's side
// TODO: (every call to next() is an expensive python function call).
// TODO: So, we might prefer to exposing some specialized functionality
// TODO: where we accept an arbitrary `ParticlesIterable` as we currently do
// TODO: but also perhaps clone into some buffer so we don't have to hold the GIL for too long
#[derive(Clone)]
struct ParticlesIterable<'py> {
    inner: Bound<'py, PyAny>,
}

impl<'py> IntoIterator for ParticlesIterable<'py> {
    type Item = [f64; 3];
    type IntoIter = ParticlesIterator<'py>;

    fn into_iter(self) -> Self::IntoIter {
        ParticlesIterator {
            // PyO3 also just `unwrap()`s in their specific iterators
            iter: self.inner.try_iter().unwrap(),
        }
    }
}

#[derive(Clone)]
struct ParticlesIterator<'py> {
    iter: Bound<'py, PyIterator>,
}

impl<'py> Iterator for ParticlesIterator<'py> {
    type Item = [f64; 3];

    fn next(&mut self) -> Option<Self::Item> {
        // TODO: document behavior:
        // while next element is some Error
        // retry getting new next element
        // if it's None, break loop and return
        // if it's Some, attempt conversion and break loop if it was successful
        // otherwise, retry with next element
        loop {
            match self.iter.next().transpose() {
                Ok(Some(p)) => match <[f64; 3] as FromPyObject>::extract_bound(&p) {
                    Ok(p) => break Some(p),
                    Err(_) => (),
                },
                Ok(None) => break None,
                Err(_) => (),
            }
        }
    }
}

/// 3D cell grid
#[derive(Clone)]
#[pyclass(name = "CellGrid", module = "zelll")]
pub struct PyCellGrid {
    inner: CellGrid<[f64; 3]>,
}

#[pymethods]
impl PyCellGrid {
    #[new]
    #[pyo3(signature = (particles=None, /, cutoff=1.0))]
    fn new<'py>(py: Python<'py>, particles: Option<&Bound<'py, PyAny>>, cutoff: f64) -> Self {
        let inner = match particles {
            Some(p) => {
                // TODO: see if we can simplify ParticlesIterable, it wraps Bound<>, so holds the GIL
                // TODO: would like to call ::new() detached from the GIL if that's possible?
                let particles = ParticlesIterable { inner: p.clone() };
                CellGrid::new(particles, cutoff)
            }
            // nutpie+multithreading needs to serialize/deserialize PyCellGrid
            // and the latter calls ::new() (albeit with default arguments)
            // so in this case we can detach from the GIL
            // (although __setstate__() attaches again)
            _ => py.detach(|| CellGrid::default()),
        };

        Self { inner }
    }

    #[pyo3(signature = (particles, /, cutoff=None))]
    fn rebuild<'py>(
        mut slf: PyRefMut<'_, Self>,
        particles: &Bound<'py, PyAny>,
        cutoff: Option<f64>,
    ) -> () {
        let particles = ParticlesIterable {
            inner: particles.clone(),
        };

        slf.inner.rebuild_mut(particles, cutoff);
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyCellGridIter {
        PyCellGridIter::new(slf)
    }

    fn aabb(slf: PyRef<'_, Self>) -> ([f64; 3], [f64; 3]) {
        (
            slf.inner.info().bounding_box().inf(),
            slf.inner.info().bounding_box().sup(),
        )
    }

    fn cutoff(slf: PyRef<'_, Self>) -> f64 {
        slf.inner.info().cutoff()
    }
    #[pyo3(signature = (coordinates: "typing.Sequence[float]") -> "typing.Iterator[typing.Sequence[float]] | None")]
    fn query_neighbors(
        slf: PyRef<'_, Self>,
        coordinates: &Bound<'_, PyAny>,
    ) -> Option<PyCellQueryIter> {
        PyCellQueryIter::new(slf, coordinates)
    }

    // TODO: document that this filters by inner cutoff and returns Vec<> instead
    fn neighbors(slf: PyRef<'_, Self>, coordinates: [f64; 3]) -> Option<Vec<(usize, [f64; 3])>> {
        let cutoff_squared = slf.inner.info().cutoff().powi(2);
        let iter = slf.inner.query_neighbors(coordinates)?;
        slf.py().detach(|| {
            // filter by squared euclidean distance
            let out = iter.filter(|(_, other)| {
                let [x, y, z] =
                    std::array::from_fn(|i| coordinates[i] - other[i]).map(|diff| diff * diff);
                x + y + z <= cutoff_squared
            });
            Some(out.collect())
        })
    }

    fn __getstate__(&self, py: Python) -> PyResult<Py<PyAny>> {
        let bytes = encode_to_vec(&self.inner, bincode::config::standard()).unwrap();
        PyBytes::new(py, &bytes).into_py_any(py)
    }

    fn __setstate__(&mut self, state: &Bound<'_, PyAny>) -> PyResult<()> {
        // we can extract &[u8] directly because PyBytes transparently wraps PyAny
        // and that's how we produced `state` through `__getstate__()`
        let bytes = state.extract::<&[u8]>()?;
        self.inner = match decode_from_slice(bytes, bincode::config::standard()) {
            Ok((inner, _)) => Ok(inner),
            Err(_) => Err(PyErr::new::<PyTypeError, _>(
                "Could not unpickle CellGrid from this type",
            )),
        }?;
        Ok(())
    }
}

// PyCellGridIter is unsendable because we need to store a PyRef directly.
// cf. https://docs.rs/pyo3/latest/pyo3/attr.pyclass.html
#[pyclass(name = "CellGridIter", module = "zelll", unsendable)]
pub struct PyCellGridIter {
    // TODO: it looks like we probably don't need `_owner`
    // TODO: also haven't yet figured out how to migrate this from IntoPy to IntoPyObject
    // TODO: likely need more unsafe using Py::from_*_pointer()
    // _owner: PyObject,
    // TODO: `_keep_borrow` is enough to maintain correct drop order *and* prevents `PyCellGrid`
    // TODO: from being mutated while `PyCellGridIter` is still alive
    _keep_borrow: PyRef<'static, PyCellGrid>,
    iter: Box<dyn Iterator<Item = ((usize, [f64; 3]), (usize, [f64; 3]))>>,
}

impl PyCellGridIter {
    fn new(py_cellgrid: PyRef<'_, PyCellGrid>) -> Self {
        // let py = py_cellgrid.py();
        // let _owner = (&py_cellgrid).into_py(py);
        let iter = Box::new((&py_cellgrid).inner.particle_pairs());
        // SAFETY: unclear
        // SAFETY: (but the idea is that `_keep_borrow` makes sure that `iter`s lifetime can be extended)
        // replicating some ideas from
        // https://github.com/PyO3/pyo3/issues/1085 and
        // https://github.com/PyO3/pyo3/issues/1089
        let iter = unsafe {
            std::mem::transmute::<
                Box<dyn Iterator<Item = ((usize, [f64; 3]), (usize, [f64; 3]))> + '_>,
                Box<dyn Iterator<Item = ((usize, [f64; 3]), (usize, [f64; 3]))> + 'static>,
            >(iter)
        };

        // SAFETY: unclear
        // SAFETY: (the idea was that `_owner` ensures our static reference is valid as long as we hold `_owner`
        // SAFETY: but experiments show that it does not seem to be necessary.
        // SAFETY: Even if the owning `PyCellGrid` is dropped, ie. goes out of scope or using `del`, the GC
        // SAFETY: seems to keep it until the last `PyCellGridIter` is dropped but
        // SAFETY: I guess this should be profiled somehow)
        // SAFETY: (this might depend on `PyRef` being a wrapper around `Bound`, which might change to `Borrowed`,
        // SAFETY: could this cause problems?)
        let _keep_borrow: PyRef<'static, PyCellGrid> = unsafe { std::mem::transmute(py_cellgrid) };

        Self {
            // _owner,
            _keep_borrow,
            iter,
        }
    }
}

// impl Drop for PyCellGrid {
//     fn drop(&mut self) {
//         eprintln!("Dropping PyCellGrid");
//     }
// }

// impl Drop for PyCellGridIter {
//     fn drop(&mut self) {
//         eprintln!("Dropping PyCellGridIter");
//     }
// }

#[pymethods]
impl PyCellGridIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<((usize, [f64; 3]), (usize, [f64; 3]))> {
        slf.iter.next()
    }
}

// PyCellQueryIter is unsendable because we need to store a PyRef directly.
#[pyclass(name = "CellQueryIter", module = "zelll", unsendable)]
pub struct PyCellQueryIter {
    _keep_borrow: PyRef<'static, PyCellGrid>,
    iter: Box<dyn Iterator<Item = (usize, [f64; 3])>>,
}

impl PyCellQueryIter {
    fn new(py_cellgrid: PyRef<'_, PyCellGrid>, coordinates: &Bound<'_, PyAny>) -> Option<Self> {
        let coordinates = <[f64; 3] as FromPyObject>::extract_bound(coordinates).ok()?;
        let iter = Box::new((&py_cellgrid).inner.query_neighbors(coordinates)?);
        // SAFETY: see PyCellGridIter
        let iter = unsafe {
            std::mem::transmute::<
                Box<dyn Iterator<Item = (usize, [f64; 3])> + '_>,
                Box<dyn Iterator<Item = (usize, [f64; 3])> + 'static>,
            >(iter)
        };

        // SAFETY: see PyCellGridIter
        let _keep_borrow: PyRef<'static, PyCellGrid> = unsafe { std::mem::transmute(py_cellgrid) };

        Some(Self { _keep_borrow, iter })
    }
}

#[pymethods]
impl PyCellQueryIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(usize, [f64; 3])> {
        slf.iter.next()
    }
}

#[doc = include_str!("../README.md")]
#[pymodule]
pub mod zelll {
    #[pymodule_export]
    pub use super::{PyCellGrid, PyCellGridIter, PyCellQueryIter};
}
