//TODO: document motivation for this
//TODO: in principle we could solve this using slices instead of CellSliceMeta
//TODO: but then we'd have to distingquish mutability when storing the slices
//TODO: Also, bumpalo would be a nice approach too but there we'd have
//TODO: to use allocator-api2 and resetting bumpalo is not nice (again because of mutability)
//TODO: So instead let's just use a simple wrapper type around Vec with some specialized functionality
//TODO: This wrapper could have a field len/cursor s.t. resizing and stuff does not drop values
//TODO: but essentially we're doing the same things as Vec anyway and usize doesn't implement Drop
//TODO: The downside of this approach (compared to e.g. Bumpalo) is that we have less guarantees
//TODO: CellSliceMeta is decoupled from CellStorage (on purpose)
//TODO: we're basically doing "unsafe" stuff using indices instead of pointers
//TODO: (resulting in potential `panic!()`s instead of unsound stuff. Well, unsound stuff results now in logic errors)
//TODO: Basically my only issue with bumpalo is that I can't have a persistent HashMap<usize, Vec<usize, &Bump>>
//TODO: in CellGrid. I'd have to re-allocate this whenever CellGrid is rebuilt
//TODO: but maybe this is not an issue anyway. This CellStorage wrapper type only solves this issue
//TODO: Also hashbrown supports bumpalo, so I could have a separate Bump instance if I'd care about this
//TODO: or can I do this:
//TODO: https://users.rust-lang.org/t/reuse-a-vec-t-s-allocation-beyond-the-scope-of-its-content/66398/4
//TODO: with hashbrown::HashMap? So far I didn't manage to do it
//TODO: I can:
//TODO: https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&code=use+hashbrown%3A%3AHashMap%3B%0A%0Afn+main%28%29+%7B%0A++++let+mut+empty_hm%3A+HashMap%3Cusize%2C+usize%3E+%3D+HashMap%3A%3Anew%28%29%3B%0A%0A++++for+i+in+0..10+%7B%0A++++++++let+mut+hm+%3D+empty_hm%3B%0A%0A++++++++%2F%2F+Fill+the+hashmap+with+some+data+bound+to+the+scope+of+the+loop.%0A++++++++for+j+in+0..i+%7B%0A++++++++hm.insert%28j%2C+i%29%3B++++%0A++++++++%7D%0A++++++++%0A++++%0A++++++++%2F%2F+sanity+check%3A+address+stays+the+same%0A++++++++println%21%28%22%7B%3Ap%7D%22%2C+%26hm%29%3B%0A++++%0A++++++++%2F%2F+Do+things+with+the+data.%0A++++++++for+s+in+%26hm+%7B+%0A++++++++++++std%3A%3Ahint%3A%3Ablack_box%28s%29%3B%0A++++++++%7D%0A++++%0A++++++++%2F%2F+Clear+the+vector%2C+ensuring+there+is+no+references+left+in+the+vector.%0A++++++++hm.clear%28%29%3B%0A++++++++empty_hm+%3D+hm.into_iter%28%29.map%28%7C_%7C+unreachable%21%28%29%29.collect%28%29%3B%0A++++%7D%0A%7D
//TODO: but I need to check if this is still possible with values of type Vec<_, &Bump>
//TODO: and see how this works without without this scoping
use core::ops::Range;

#[derive(Debug, Default, Clone)]
pub struct CellStorage<T> {
    buffer: Vec<T>,
}

impl<T> CellStorage<T> {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
        }
    }

    //TODO: provide fallible version of this
    pub fn cell_slice(&self, metadata: &CellSliceMeta) -> &[T] {
        &self.buffer[metadata.range.clone()]
    }

    //TODO: choose appropriate Error type
    pub fn try_push(&mut self, value: T, metadata: &mut CellSliceMeta) {
        todo!()
    }

    //TODO: `panic!()`s if OOB
    pub fn push(&mut self, value: T, metadata: &mut CellSliceMeta) {
        let slice = &mut self.buffer[metadata.range.clone()];
        slice[metadata.cursor] = value;
        //TODO: use ::move_cursor(1)?
        metadata.cursor += 1;
    }

    //TODO: potentially makes existing slice metadata unsound
    //TODO: does not shrink capacity
    //TODO: generational indices/ranges for metadata would make sense here but
    //TODO: I'm not trying to reinvent ECS etc. I just want a simple wrapped Vec
    //TODO: actually might store a &'s CellStorage inside of CellSliceMeta<'s>?
    //TODO: I think this is no problem (well I think it is though...) even when I'm handling &mut CellSliceMeta?
    //TODO: Then CellSliceMeta would be tied to specific storage
    pub fn truncate(&mut self, len: usize) {
        self.buffer.truncate(len);
    }

    //TODO: this does not overwrite any memory
    //TODO: document the behaviour and intention clearly
    pub fn clear(&mut self) {
        self.buffer.clear()
    }
}

impl<T: Default> CellStorage<T> {
    //TODO: this resizes dynamically (but this only happens if we add particles to the point cloud)
    pub fn reserve_cell(&mut self, capacity: usize) -> CellSliceMeta {
        let range = self.buffer.len()..(self.buffer.len() + capacity);
        self.buffer.resize_with(range.end, Default::default);

        CellSliceMeta::new(range)
    }
}

//TODO: this type does not check bounds, this is responsibility of CellStorage
#[derive(Debug, Default, Clone)]
pub struct CellSliceMeta {
    cursor: usize,
    //TODO: Range is not Copy
    //TODO: see https://github.com/rust-lang/rust/pull/27186
    range: Range<usize>,
}

impl CellSliceMeta {
    fn new(range: Range<usize>) -> Self {
        Self { cursor: 0, range }
    }
    //TODO: probably won't need this but in principle, we just reset the cursor to clear
    pub fn clear(&mut self) {
        self.cursor = 0;
    }

    //TODO: panics on OOB
    pub fn move_cursor(&mut self, steps: usize) {
        //TODO: could use Range::contains() here but I think that would be actually more complicated?
        //Range<usize> impl ExactSizeIterator (and TrustedLen)
        assert!(self.cursor + steps < self.range.len());
        self.cursor += steps;
    }

    //TODO: proper error type?
    pub fn try_move_cursor(&mut self, steps: usize) {
        todo!()
    }
}

