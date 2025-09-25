// TODO: document motivation for this
// TODO: in principle we could solve this using slices instead of CellSliceMeta
// TODO: but then we'd have to distinguish mutability when storing the slices
// TODO: Also, bumpalo would be a nice approach too but there we'd have
// TODO: to use allocator-api2 and resetting bumpalo is not nice (again because of mutability)
// TODO: So instead let's just use a simple wrapper type around Vec with some specialized functionality
// TODO: This wrapper could have a field len/cursor s.t. resizing and stuff does not drop values
// TODO: but essentially we're doing the same things as Vec anyway and usize doesn't implement Drop
// TODO: The downside of this approach (compared to e.g. Bumpalo) is that we have less guarantees
// TODO: CellSliceMeta is decoupled from CellStorage (on purpose)
// TODO: we're basically doing "unsafe" stuff using indices instead of pointers
// TODO: (resulting in potential `panic!()`s instead of unsound stuff. Well, unsound stuff results now in logic errors)
// TODO: Basically my only issue with bumpalo is that I can't have a persistent HashMap<usize, Vec<usize, &Bump>>
// TODO: in CellGrid. I'd have to re-allocate this whenever CellGrid is rebuilt
// TODO: but maybe this is not an issue anyway. This CellStorage wrapper type only solves this issue
// TODO: Also hashbrown supports bumpalo, so I could have a separate Bump instance if I'd care about this
// TODO: or can I do this:
// TODO: https://users.rust-lang.org/t/reuse-a-vec-t-s-allocation-beyond-the-scope-of-its-content/66398/4
// TODO: with hashbrown::HashMap? So far I didn't manage to do it
// TODO: I can:
// TODO: https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&code=use+hashbrown%3A%3AHashMap%3B%0A%0Afn+main%28%29+%7B%0A++++let+mut+empty_hm%3A+HashMap%3Cusize%2C+usize%3E+%3D+HashMap%3A%3Anew%28%29%3B%0A%0A++++for+i+in+0..10+%7B%0A++++++++let+mut+hm+%3D+empty_hm%3B%0A%0A++++++++%2F%2F+Fill+the+hashmap+with+some+data+bound+to+the+scope+of+the+loop.%0A++++++++for+j+in+0..i+%7B%0A++++++++hm.insert%28j%2C+i%29%3B++++%0A++++++++%7D%0A++++++++%0A++++%0A++++++++%2F%2F+sanity+check%3A+address+stays+the+same%0A++++++++println%21%28%22%7B%3Ap%7D%22%2C+%26hm%29%3B%0A++++%0A++++++++%2F%2F+Do+things+with+the+data.%0A++++++++for+s+in+%26hm+%7B+%0A++++++++++++std%3A%3Ahint%3A%3Ablack_box%28s%29%3B%0A++++++++%7D%0A++++%0A++++++++%2F%2F+Clear+the+vector%2C+ensuring+there+is+no+references+left+in+the+vector.%0A++++++++hm.clear%28%29%3B%0A++++++++empty_hm+%3D+hm.into_iter%28%29.map%28%7C_%7C+unreachable%21%28%29%29.collect%28%29%3B%0A++++%7D%0A%7D
// TODO: but I need to check if this is still possible with values of type Vec<_, &Bump>
// TODO: and see how this works without without this scoping
use core::ops::Range;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
// TODO: learn from https://davidlattimore.github.io/posts/2025/09/02/rustforge-wild-performance-tricks.html

// TODO: impl<T> Deref<Target = Vec<T>> for CellStorage<T>
// TODO: impl<T> DerefMut<Target = Vec<T>> for CellStorage<T>
// FIXME: remove methods overlapping with Vec<T>/&[T]
// FIXME: rename push() so it does not collide with Vec<T>
// TODO: expose CellStorage field in CellGrid, so users can use std::mem::swap to cheaply obtain stored particles?
// TODO: (for simulations, P: Particle should then also contain (references to) velocities/accelerations, position indices
// TODO: because those should be associated with P after reordering/swapping)
// TODO: (also tricky for HMC/NUTS, where we shouldn't shuffle parameters around)
// TODO: problem with swapping: invalidates every `CellSliceMeta`
// TODO: at least document that CellStorage allows to Deref for cheap swapping if the order of the data does not matter to the user
// TODO: because they don't do e.g. HMC or just store references/std::cell::* types
// TODO: (but the whole point of storing P: Particle in CellStorage is to avoid cache misses due to references...)
// TODO: could also make CellStorage<T> generic over buffer and require Index/IndexMut/SliceIndex traits
// TODO: https://doc.rust-lang.org/std/boxed/struct.Box.html#method.new_uninit_slice would be quite nice
// TODO: but then we couldn't add more particles to storage/grid
// FIXME: in the future we might change this to one Vec<> per GridCell to address
// FIXME: performance issues in CellGrid::new()/::rebuild_mut()
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub(crate) struct CellStorage<T> {
    buffer: Vec<T>,
}

impl<T> CellStorage<T> {
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
        }
    }

    // TODO: provide fallible version of this
    // FIXME: not taking `&CellSliceMeta` because we'd have to `clone()` the internal range
    // FIXME: now the caller has to do it themselves, so at least they are aware of it
    // FIXME: waiting for new copyable range types
    pub(crate) fn cell_slice(&self, metadata: CellSliceMeta) -> &[T] {
        &self.buffer[metadata.range]
    }

    pub(crate) fn try_cell_slice(&self, metadata: CellSliceMeta) -> Option<&[T]> {
        self.buffer.get(metadata.range)
    }

    // TODO: choose appropriate Error type
    pub(crate) fn try_push(&mut self, _value: T, _metadata: &mut CellSliceMeta) {
        todo!()
    }

    // TODO: `panic!()`s if OOB
    pub(crate) fn push(&mut self, value: T, metadata: &mut CellSliceMeta) {
        let slice = &mut self.buffer[metadata.range.clone()];
        slice[metadata.cursor] = value;
        metadata.move_cursor(1);
    }

    // TODO: potentially makes existing slice metadata unsound
    // TODO: does not shrink capacity
    // TODO: generational indices/ranges for metadata would make sense here but
    // TODO: I'm not trying to reinvent ECS etc. I just want a simple wrapped Vec
    // TODO: actually might store a &'s CellStorage inside of CellSliceMeta<'s>?
    // TODO: I think this is no problem (well I think it is though...) even when I'm handling &mut CellSliceMeta?
    // TODO: cf. https://github.com/CAD97/generativity but I'd like to avoid lifetime trickery
    // TODO: Then CellSliceMeta would be tied to specific storage
    // TODO: So, maybe CellStorage/CellSliceMeta should not really be exposed to public-facing API
    // TODO: in some sense, real `unsafe` for CellStorage would be more honest
    pub(crate) fn truncate(&mut self, len: usize) {
        self.buffer.truncate(len);
    }

    // TODO: this does not overwrite any memory
    // TODO: document the behaviour and intention clearly
    pub(crate) fn clear(&mut self) {
        self.buffer.clear()
    }
}

impl<T: Default> CellStorage<T> {
    // TODO: this resizes dynamically (but this only happens if we add particles to the point cloud)
    pub(crate) fn reserve_cell(&mut self, capacity: usize) -> CellSliceMeta {
        let range = self.buffer.len()..(self.buffer.len() + capacity);
        self.buffer.resize_with(range.end, Default::default);

        CellSliceMeta::new(range)
    }
}

// TODO: this type does not check bounds, this is responsibility of CellStorage
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub(crate) struct CellSliceMeta {
    cursor: usize,
    // TODO: Range is not Copy
    // TODO: see https://github.com/rust-lang/rust/pull/27186
    range: Range<usize>,
}
// TODO: note that if I make this Either<i32, CellSliceMeta>, I could also remove CellSliceMeta entirely
// TODO: if I create a  custom Enum for us:
// enum CellSlice {
//     Capacity(i32),
//     Slice(&[usize]),
//     MutSlice(&mut [usize]),
// }
// TODO: would then need to convert from Capacity(i32) to MutSlice by repeatedly calling split_at_mut()

impl CellSliceMeta {
    pub(crate) fn new(range: Range<usize>) -> Self {
        Self { cursor: 0, range }
    }
    // TODO: probably won't need this but in principle, we just reset the cursor to clear
    pub(crate) fn clear(&mut self) {
        self.cursor = 0;
    }

    // TODO: document OOB behavior (based on debug_assert/assert and CellStorage::push() behavior)
    pub(crate) fn move_cursor(&mut self, steps: usize) {
        // FIXME: this will be triggered on CellGrid::rebuild(_mut)()
        // FIXME: removing this debug_assert! for now since CellStorage is not part of the public API
        // FIXME: but CellStorage::push() does perform bounds checks on the actual slice anyway
        // debug_assert!(self.cursor + steps < self.range.len());
        self.cursor += steps;
    }
    // panics if moving back causes OOB (i.e. if cursor is currently 0)
    pub(crate) fn move_cursor_back(&mut self, steps: usize) {
        // FIXME: this will be triggered on CellGrid::rebuild(_mut)()
        // FIXME: removing this debug_assert! for now since CellStorage is not part of the public API
        // FIXME: but CellStorage::push() does perform bounds checks on the actual slice anyway
        // debug_assert!(self.cursor + steps < self.range.len());
        self.cursor -= steps;
    }

    //FIXME: technically a len() but might be misleading because cursor might exceed range.end...
    pub(crate) fn cursor(&self) -> usize {
        self.cursor
    }

    // TODO: proper error type?
    pub(crate) fn try_move_cursor(&mut self, _steps: usize) {
        todo!()
    }
}

use core::iter::FusedIterator;

#[derive(Debug, Default, Clone)]
pub(crate) struct DenseMap<K, V> {
    // TODO: should make this a Vec<Option<(i32, CellSliceMeta)>>
    inner: Vec<Option<(K, V)>>,
}

impl<K, V> DenseMap<K, V> {
    pub fn new() -> Self {
        Self { inner: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            inner: Vec::with_capacity(capacity),
        }
    }
}

impl<V: Default + Clone> GridStorage<V> for DenseMap<i32, V> {
    fn clear(&mut self) {
        self.inner.clear();
    }

    fn get(&self, k: &i32) -> Option<&V> {
        self.inner
            .get(*k as usize)
            .and_then(|c| c.as_ref().map(|(_, v)| v))
    }

    fn get_mut(&mut self, k: &i32) -> Option<&mut V> {
        self.inner
            .get_mut(*k as usize)
            .and_then(|c| c.as_mut().map(|(_, v)| v))
    }

    fn shrink_to_fit(&mut self) {
        self.inner.shrink_to_fit();
    }

    // FIXME: need entry API here instead
    fn insert(&mut self, k: i32, v: V) {
        if k as usize >= self.inner.len() {
            self.inner.resize(k as usize + 1, None);
        }
        self.inner[k as usize] = Some((k, v));
    }

    // TODO: typical vecmap: https://crates.io/crates/vecmap-rs
    // TODO: impl GridStorage trait for VecMap, hashbrown::HashMap, std HashMap, std BTreeMap
    // TODO: and sprs::CsVecBase
    // TODO: note that vecmap-rs could be used as a sparse storage with costly insert(), similar to CsVecBase
    // TODO: would have to prepulate vecmap-rs to make it a dense storage...
    // TODO: or simply keep my wrapper DenseMap
    // TODO: https://docs.rs/ordered-vecmap/latest/ordered_vecmap/ would also be interesting
    fn iter<'a>(&'a self) -> impl FusedIterator<Item = (&'a i32, &'a V)> + Clone + 'a
    where
        V: 'a,
    {
        self.inner
            .iter()
            .filter_map(|c| c.as_ref().map(|(k, v)| (k, v)))
    }

    fn keys<'a>(&'a self) -> impl FusedIterator<Item = &'a i32> + Clone + 'a
    where
        V: 'a,
    {
        self.inner.iter().filter_map(|c| c.as_ref().map(|(k, _)| k))
    }
}

use std::iter::FromIterator;
use std::ops::Index;

impl<V: Clone + Default> Index<&i32> for DenseMap<i32, V> {
    type Output = V;

    fn index(&self, index: &i32) -> &Self::Output {
        self.get(index).expect("index should not be out of bounds.")
    }
}

impl<V: Clone + Default> FromIterator<(i32, V)> for DenseMap<i32, V> {
    fn from_iter<T: IntoIterator<Item = (i32, V)>>(iter: T) -> Self {
        let iter = iter.into_iter();
        let mut dense_map = DenseMap::with_capacity(iter.size_hint().0);
        iter.for_each(|(k, v)| dense_map.insert(k, v));
        dense_map
    }
}

// FIXME: decide whether keep the trait methods' names as-is
// FIXME: or find more unique names
// FIXME: might be a bit inconvenient to use because
// FIXME: we might have to use fully qualified syntax?
// FIXME: actually that might not be a problem?
// FIXME: it's rather users of `zelll` if they decide to change the storage type
// FIXME: at the very least, we shouldn't re-export GridStorage?
// TODO: might need GridStorage: FromIterator<(K, V)> + Default
pub(crate) trait GridStorage<V = CellSliceMeta>: Default {
    // TODO: make Entry a generic parameter or hardcode it to CellSliceMeta?
    // default implementation is a no-op
    // but overriding this helps to replicate with_capacity() (since Default is super trait)
    fn reserve(&mut self, _additional: usize) {}

    fn shrink_to_fit(&mut self) {}

    // TODO: documentation:
    // TODO: clear this `impl GridStorage`
    // TODO: this is required by implementors
    // TODO: so that a grid storage type does not accumulate stale keys/cells
    fn clear(&mut self);

    // FIXME: actually would rather like ::entry() API here
    fn insert(&mut self, k: i32, v: V);

    fn get(&self, k: &i32) -> Option<&V>;

    fn get_mut(&mut self, k: &i32) -> Option<&mut V>;

    // FIXME: need also iter_mut()
    fn iter<'a>(&'a self) -> impl FusedIterator<Item = (&'a i32, &'a V)> + Clone + 'a
    where
        V: 'a;

    fn keys<'a>(&'a self) -> impl FusedIterator<Item = &'a i32> + Clone + 'a
    where
        V: 'a,
    {
        GridStorage::iter(self).map(|(k, _)| k)
    }
}

use hashbrown::HashMap;
use std::collections::{BTreeMap, HashMap as StdHashMap};

impl<V> GridStorage<V> for HashMap<i32, V> {
    fn reserve(&mut self, additional: usize) {
        self.reserve(additional);
    }

    fn shrink_to_fit(&mut self) {
        self.shrink_to_fit();
    }

    fn clear(&mut self) {
        self.clear()
    }

    fn insert(&mut self, k: i32, v: V) {
        // SAFETY:
        // This is safe because CellGrid only inserts unique keys *once*
        // or clears the whole impl GridStorage before inserting the same key again.
        unsafe {
            self.insert_unique_unchecked(k, v);
        }
    }

    fn get(&self, k: &i32) -> Option<&V> {
        self.get(k)
    }

    fn get_mut(&mut self, k: &i32) -> Option<&mut V> {
        self.get_mut(k)
    }

    fn iter<'a>(&'a self) -> impl FusedIterator<Item = (&'a i32, &'a V)> + Clone + 'a
    where
        V: 'a,
    {
        self.iter()
    }

    fn keys<'a>(&'a self) -> impl FusedIterator<Item = &'a i32> + Clone + 'a
    where
        V: 'a,
    {
        self.keys()
    }
}

impl<V> GridStorage<V> for StdHashMap<i32, V> {
    fn reserve(&mut self, additional: usize) {
        self.reserve(additional);
    }

    fn shrink_to_fit(&mut self) {
        self.shrink_to_fit();
    }

    fn clear(&mut self) {
        self.clear()
    }

    fn insert(&mut self, k: i32, v: V) {
        self.insert(k, v);
    }

    fn get(&self, k: &i32) -> Option<&V> {
        self.get(k)
    }

    fn get_mut(&mut self, k: &i32) -> Option<&mut V> {
        self.get_mut(k)
    }

    fn iter<'a>(&'a self) -> impl FusedIterator<Item = (&'a i32, &'a V)> + Clone + 'a
    where
        V: 'a,
    {
        self.iter()
    }

    fn keys<'a>(&'a self) -> impl FusedIterator<Item = &'a i32> + Clone + 'a
    where
        V: 'a,
    {
        self.keys()
    }
}

impl<V> GridStorage<V> for BTreeMap<i32, V> {
    fn clear(&mut self) {
        self.clear()
    }

    fn insert(&mut self, k: i32, v: V) {
        self.insert(k, v);
    }

    fn get(&self, k: &i32) -> Option<&V> {
        self.get(k)
    }

    fn get_mut(&mut self, k: &i32) -> Option<&mut V> {
        self.get_mut(k)
    }

    fn iter<'a>(&'a self) -> impl FusedIterator<Item = (&'a i32, &'a V)> + Clone + 'a
    where
        V: 'a,
    {
        self.iter()
    }

    fn keys<'a>(&'a self) -> impl FusedIterator<Item = &'a i32> + Clone + 'a
    where
        V: 'a,
    {
        self.keys()
    }
}
