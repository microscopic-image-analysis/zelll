use ::zelll::cellgrid::*;
use pyo3::prelude::*;

/// 3D cell grid
#[pyclass(name = "CellGrid")]
pub struct PyCellGrid {
    inner: CellGrid<3>,
}

#[pymethods]
impl PyCellGrid {
    #[new]
    fn new() -> Self {
        let points: Vec<[f64; 3]> = vec![
            [0.534, -0.491, -0.22],
            [0.108, 0.954, 0.221],
            [-0.629, -0.916, -0.009],
            [-0.9, -0.15, 0.017],
            [-0.01, -0.323, -0.173],
            [-0.382, 0.329, 0.912],
            [-0.609, -0.665, -0.004],
            [0.939, 0.129, 0.676],
            [0.544, 0.397, 0.248],
            [0.185, 0.609, 0.12],
        ];
        Self {
            inner: CellGrid::new(&points, 1.0),
        }
    }

    fn rebuild(mut slf: PyRefMut<'_, Self>) -> () {
        let points: Vec<[f64; 3]> = vec![
            [0.534, -0.491, -0.22],
            [0.108, 0.954, 0.221],
            [-0.629, -0.916, -0.009],
            [-0.9, -0.15, 0.017],
            [-0.01, -0.323, -0.173],
            [-0.382, 0.329, 0.912],
            [-0.609, -0.665, -0.004],
            [0.939, 0.129, 0.676],
            [0.544, 0.397, 0.248],
            [0.185, 0.609, 0.12],
        ];
        slf.inner.rebuild_mut(points.iter().rev(), None);
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyCellGridIter {
        PyCellGridIter::new(slf)
    }
}

// PyCellGridIter is unsendable because we need to store a PyRef directly.
// cf. https://docs.rs/pyo3/latest/pyo3/attr.pyclass.html
#[pyclass(name = "CellGridIter", unsendable)]
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
