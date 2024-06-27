use ::zelll::cellgrid::*;
use pyo3::prelude::*;
use pyo3::types::PyIterator;

#[derive(Clone)]
struct PointsIterable<'py> {
    inner: Bound<'py, PyAny>,
}

impl<'py> IntoIterator for PointsIterable<'py> {
    type Item = [f64; 3];
    type IntoIter = PointsIterator<'py>;

    fn into_iter(self) -> Self::IntoIter {
        PointsIterator {
            // panicking is probably not the best idea but I can't bubble up an exception to the rust-python boundary.
            // and PyO3 seems to do the same for their specific iterators
            iter: self.inner.iter().unwrap(),
        }
    }
}

#[derive(Clone)]
struct PointsIterator<'py> {
    iter: Bound<'py, PyIterator>,
}

impl<'py> Iterator for PointsIterator<'py> {
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
                Ok(Some(point)) => match <[f64; 3] as FromPyObject>::extract_bound(&point) {
                    Ok(point) => break Some(point),
                    Err(_) => (),
                },
                Ok(None) => break None,
                Err(_) => (),
            }
        }
    }
}

/// 3D cell grid
#[pyclass(name = "CellGrid", module = "zelll")]
pub struct PyCellGrid {
    inner: CellGrid<3>,
}

#[pymethods]
impl PyCellGrid {
    #[new]
    fn new<'py>(points: &Bound<'py, PyAny>, cutoff: f64) -> Self {
        let points = PointsIterable {
            inner: points.clone(),
        };

        Self {
            inner: CellGrid::new(points, cutoff),
        }
    }

    #[pyo3(signature = (points, /, cutoff=None))]
    fn rebuild<'py>(
        mut slf: PyRefMut<'_, Self>,
        points: &Bound<'py, PyAny>,
        cutoff: Option<f64>,
    ) -> () {
        let points = PointsIterable {
            inner: points.clone(),
        };

        slf.inner.rebuild_mut(points, cutoff);
        // slf.inner = (&slf.inner).clone().rebuild(points, cutoff);
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyCellGridIter {
        PyCellGridIter::new(slf)
    }
}

// PyCellGridIter is unsendable because we need to store a PyRef directly.
// cf. https://docs.rs/pyo3/latest/pyo3/attr.pyclass.html
#[pyclass(name = "CellGridIter", module = "zelll", unsendable)]
pub struct PyCellGridIter {
    #[pyo3(get)]
    pub owner: PyObject,
    _uphold_borrow: PyRef<'static, PyCellGrid>,
    iter: Box<dyn Iterator<Item = (usize, usize)>>,
}

impl PyCellGridIter {
    fn new(py_cellgrid: PyRef<'_, PyCellGrid>) -> Self {
        let py = py_cellgrid.py();
        let owner = (&py_cellgrid).into_py(py);
        let iter = Box::new((&py_cellgrid).inner.point_pairs());
        // SAFETY: lol
        // replicating some ideas from
        // https://github.com/PyO3/pyo3/issues/1085
        let iter = unsafe {
            std::mem::transmute::<
                Box<dyn Iterator<Item = (usize, usize)> + '_>,
                Box<dyn Iterator<Item = (usize, usize)> + 'static>,
            >(iter)
        };
        // SAFETY: lol
        let static_borrow: PyRef<'static, PyCellGrid> = unsafe { std::mem::transmute(py_cellgrid) };

        Self {
            owner,
            _uphold_borrow: static_borrow,
            iter,
        }
    }
}

#[pymethods]
impl PyCellGridIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<(usize, usize)> {
        slf.iter.next()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn zelll(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCellGrid>()?;
    m.add_class::<PyCellGridIter>()?;
    Ok(())
}
