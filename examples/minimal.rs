use allocator_api2::boxed::Box as BumpBox;
use allocator_api2::vec::Vec as BumpVec;
use bumpalo::Bump;
use hashbrown::{HashMap, HashSet};
use nalgebra::distance_squared;
use nohash_hasher::BuildNoHashHasher;
use std::hint::black_box;
use std::iter::FromIterator;
use zelll::cellgrid::{generate_points_random, CellGrid};

fn main() {
    for size in (2..=6).map(|exp| 10usize.pow(exp)) {
        let cutoff: f64 = 10.0;
        let conc = 10.0 / cutoff.powi(3); //i.e. 100mol per 10^3 volume units
        let a = 3.0 * cutoff;
        let b = 3.0 * cutoff;
        let c = (size as f64 / conc) / a / b;
        let vol_edges = (size as f64 / conc).cbrt();
        let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0]);

        let cg = CellGrid::new(pointcloud.iter(), cutoff);
        println!("{:?}", cg.shape());
        let cutoff_squared = cutoff.powi(2);

        //cg.for_each_point_pair(|_, _| black_box(()));
        let mut count: usize = 0;
        cg.filter_point_pairs(
            |_, _| {
                count += 1;
            },
            |i, j| true, //distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared,
        );
        println!("{}", count);
    }

    /*let vals = 0..100usize;
    let vals = vec![
        443, 262, 52, 843, 917, 777, 323, 801, 73, 961, 780, 212, 220, 852, 619, 649, 297, 689,
        588, 617, 146, 105, 560, 481, 108, 442, 796, 299, 261, 257, 740, 541, 337, 435, 47, 197,
        545, 770, 564, 502, 753, 203, 524, 982, 29, 946, 600, 253, 558, 907, 886, 394, 491, 951,
        272, 151, 31, 904, 206, 863, 464, 772, 799, 610, 125, 501, 938, 399, 701, 44, 335, 543,
        631, 940, 759, 252, 351, 413, 978, 602, 972, 254, 193, 949, 553, 80, 324, 711, 851, 969,
        721, 308, 540, 576, 923, 875, 751, 14, 412, 27,
    ];
    let hs: HashSet<usize, BuildNoHashHasher<usize>> = HashSet::from_iter(vals);
    let collected: Vec<_> = hs.iter().collect();
    println!("{collected:?}");
    let collected: Vec<_> = hs
        .iter()
        .map(|i| (i, i % hs.raw_table().buckets()))
        .collect();
    println!("{collected:?}");
    println!(
        "{} {} {}",
        hs.len(),
        hs.capacity(),
        hs.raw_table().buckets()
    );

    let v = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let sub_lengths = vec![2, 3, 2, 1, 1, 2];
    let mut boxed = v.into_boxed_slice();

    let mut sliced_mut: Vec<&mut [usize]> = vec![];

    let mut rest = boxed.as_mut();

    for l in &sub_lengths {
        let (left, right) = rest.split_at_mut(*l);
        sliced_mut.push(left);

        rest = right;
    }

    //println!("{:?}", &boxed);

    for slice in sliced_mut {
        for elem in slice {
            *elem += 1;
        }
    }
    //TODO: https://doc.rust-lang.org/nomicon/transmutes.html
    //TODO: DON'T transmute &[usize] to &mut [usize] (other way around is still unsafe but okay)
    //TODO: when transmuting references we have to use life time annotations! (otherwise we'll get 'static)
    //TODO: we'll want to transmute HashMap<SomeIndexType, &mut [usize]> to HashMap<SomeIndexType, &[usize]>
    //TODO: can we rely on the layout of hashbrown even in repr(Rust)? (Probably not)
    //TODO: safe alternatives:
    //TODO: HashMap<SomeIndexType, Range<usize>> constructed from tuple windows of cumulative cell list particle counts:
    //TODO: vec![2,6,1,3,2] -> vec![2,8,9,12,14] -> (0,2), (2,8), (8,9), (9,12), (12,14)
    //TODO: Range<usize> does impl Iterator but this mutates the lower bound of the range
    //TODO: I.e. we temporarily need a "cursor" hash map (in any case, even with transmutes)
    //TODO: (see https://doc.rust-lang.org/std/io/struct.Cursor.html but we'd want this not for bytes but usize)
    //TODO: makes sense to make these cursors a HashMap<SomeIndexType, std::slice::IterMut<'_, usize>>
    //TODO: or store the cursor permanently, allowing to shrink individual lists
    //TODO: maybe this would allow us to implement some better rebuild heuristics (don't think so)
    //TODO: enum CellSclice to distinguish/convert between mut and normal slices
    //TODO: can use some arena allocator, like bumpalo, id-arena or typed-arena
    //TODO: maybe should have some enum like CellSlice { UnInit(usize), Init(&[usize]) }
    //TODO: could use bumpalo::collections::vec::Vec instead of slice
    //TODO: this way I wouldn't have to store a cursor and just grow the Vec inside of bumpalo (within its capacity)
    //TODO: BUT: bumpalos Vec is !Sync and !Send...
    //TODO: but I could convert it to a slice after filling (dunno if 'bump lifetime is an issue here though)
    //TODO: but enums containing a !Sync and !Send variant also are thread-unsafe...
    //TODO: (or bumpalos Box<[usize]> which is Send & Sync too)
    //TODO: enum CellSlice { UnInit(usize), Init(Vec<'bump, uize>), InitDone(&'bump [usize] or Box<'bump, [usize]>)}
    //TODO: Box<'bump, [usize]> has the advantage that I still could obtain mut slices
    //TODO: or just start with a Vec-hashmap and then use from_iter() to construct a Box<[usize]> hashmap
    //TODO: and hope that it's fast enough

    /*let mut sliced: Vec<&[usize]> = vec![];
    for slice in sliced_mut {
        sliced.push(slice as &[usize]);
    }*/

    println!("{:?}", boxed);

    let bump = Bump::with_capacity(8 * 11);
    let b = Bump::new();
    bump.set_allocation_limit(Some(8 * 11));
    //let mut sliced: Vec<&mut [usize]> = vec![];
    let mut hm: HashMap<usize, BumpVec<usize, _>> = HashMap::with_capacity(sub_lengths.len());
    for (i, l) in sub_lengths.iter().enumerate() {
        //sliced.push(bump.alloc_slice_fill_default(*l));
        hm.insert(i, BumpVec::with_capacity_in(*l, &bump));
    }

    for (i, l) in sub_lengths.iter().enumerate() {
        for j in 0..*l {
            if let Some(v) = hm.get_mut(&i) {
                v.push(j);
            }
        }
    }
    //TODO: instead of converting vec hashmap to box hashmap
    //TODO: could just create box hashmap and use mut sub_lengths vec/hashmap
    //TODO: to 1. allocate Box with given capacity, then use sub_lengths as cursor (decreasing l: &Box<[usize]>.len() - l)
    //TODO: creating the Box<[usize]> zero-initialized? or unsafe with MaybeUninit?
    //TODO: this is currently not nice since the allocator API is not yet stable
    //TODO: so for now I'll just create the temporary vec hashmap
    //TODO: (although hashbrown is bumpalo compatible; we could reuse the same memory if I knew the exact size of it in bytes)
    //TODO: but the total number of cells changes anyway between simulation steps
    //TODO: could also use allocator_api2 feature in bumpalo and use the provided Vec which is Send & Sync
    //TODO: then I could just use Vec instead of Box anyway  (or do the enum stuff again)
    //TODO: need to check if performance gain from computing each cell vec capacity beforehnd
    //TODO: is worth it or if I should just let bumpalo handle re-allocating growing vecs
    //TODO: Also maybe amortized re-allocation of std Vec's is efficient enough?
    //TODO: Also I could just let MultiIndex count the occurrences of each index (for compartmentalization, doesn't make a difference though)

    println!("{:p}", &hm);
    let mut hm: HashMap<usize, BumpBox<[usize], _>> =
        HashMap::from_iter(hm.into_iter().map(|(k, v)| (k, v.into_boxed_slice())));
    println!("{:p}", &hm);
    let test: &mut [usize] = hm.get_mut(&0).unwrap();

    let sl: &mut [usize] = bump.alloc_slice_fill_default(10);
    let test = BumpBox::new_in(sl, &bump);
    let test2: BumpVec<usize, _> = BumpVec::with_capacity_in(10, &bump);

    println!("{:?}", hm);*/
}

enum CellSlice<'s> {
    Mutable(&'s mut [usize]),
    Immutable(&'s [usize]),
}

impl<'s> CellSlice<'s> {
    fn as_immutable(&'s self) -> CellSlice<'s> {
        match self {
            CellSlice::Mutable(slice) => CellSlice::Immutable(slice as &'s [usize]),
            // can't use wildcard here since only one of the variants implements Copy
            CellSlice::Immutable(slice) => CellSlice::Immutable(slice),
        }
    }
}

