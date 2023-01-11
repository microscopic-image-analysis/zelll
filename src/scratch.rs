#![allow(dead_code)]
/// TODO: cell grid related types, i.e. flat and multi indices of cubical cells.
// TODO: I'll make this more general at a later point.
use nalgebra::Point3;
use ndarray::{Array2, Array3};

//TODO: use ndarray::Array<usize> as MultiIndex
//TODO: Neighborhood Type should represent/emulate neighbor lists (i.e. a representative particle
//TODO: for that cell + backing storage for all other particles in that cell)
//TODO: the ndarray ansatz has some limitations though: no negative indices
//TODO: could do manual computations by tracking the origin of the point cloud?
//TODO: the tricky part here is updating the Multi/GridIndex.
//TODO: Naively we'd just re-allocate a new underlying Array.
//TODO: While I'm going to do this first, we should keep in mind that maybe there's a more clever way?
//TODO: Also: MultiIndex should be sparse (no need to store empty cells explicitly)

//TODO: this is a bit of a leap but maybe keep https://www.sciencedirect.com/science/article/pii/S0010465512000057
//TODO: in mind. Can I do a nested MultiIndex using ndarray::Array<N>?
//TODO: (i.e. octrees in 3D, orthtrees/hyperoctrees in nD. Also see KD-tree-pyramids or just KD trees as an alternative)
//TODO: Ok, I don't think it makes sense to nest ndarray::Array's. Just directly use a KD-tree for the cells instead.
//TODO: Can we simplify KD-trees for that? We don't need a lot, essentially just the hierarchy. nearest neighbors are easy on a cell level, I assume
//TODO: (i.e. every cell at the same depth of the KD-tree is near to the others?)
//TODO: this would call for a Multi/GridIndex trait and implementors take care of the backing storage.
//TODO: see https://en.wikipedia.org/wiki/Spatial_database#Spatial_index
//TODO: OR: have an Order/FlatIndex/FlatOrder trait that defines the, well, order of the MultiIndex, or rather how MultiIndexes are flattened
//TODO: We'd only have one MultiIndex type and Implementors would take care of the order (could internally use kd-trees for nested ordering?)

//TODO: A General Multi/GridIndex should be independent of the precise Cell geometry as long as cells fill/tile the space
//TODO: see https://en.wikipedia.org/wiki/Crystal_system
//TODO: Therefore it might make sense to have a "CrystalSystem"/Lattice? trait + implementing types defining the mapping between MultiIndex and particle coordinates
//TODO: Q: Does this translate to higher dimension though? A: Not generally, I guess
//TODO: see https://en.wikipedia.org/wiki/Lattice_(group)

//TODO: All this lattice stuff is only important for periodic boundary conditions, I guess?
//TODO: cutoff radius always defines a sphere. The less the lattice looks like cubes, the more spurious distance computations?
//TODO: Or this depends on the definition of distance. (we should allow more general distance definitions)

//TODO: this is less part of the core algorithm/problem:
//TODO: types of particles (this influences how the distance-dependent potential will be computed)

//TODO: This will probably not be a real "list". So we should rename it at some point
//TODO: the interplay between this and `CellList` is important.
//TODO: In the end we want to iterate over pairs of potential neighbors
/*
struct CellGrid {
    //TODO: issue: do we want to store empty cells (might make updating `CellList`s and indexing easier)
    //TODO: if not, we'd probably have to always rebuild a `CellGrid` on updating the pointcloud/`MultiIndex`
    cells: Vec<CellList>, // Vec<CellList>? Or can I do this more efficiently using ndarray or some BTreeMap while keeping it ergonomic?
    index: MultiIndex,    //Do I need to store this or should MultiIndex be a CellList builder type?
}

//TODO: Here's the question do we only operate on indices of the original point cloud vec?
//TODO: The Good: only need to store indices, don't need to copy the actual positions
//TODO: The Bad: We rely on invariance of the point order in the cloud vec, computing potentials always requires indexing the original point cloud
//TODO: which makes building ergonomic iterators cumbersome, I guess
//TODO: also: by storing points in the same cell contiguously probably benefits from cache locality
//TODO: BUT: updating the point cloud means we need to do more work here than with indices? Or don't we?
//TODO: The Ugly: if `CellGrid` owns the point cloud data, I probably could/should use references instead of indices? -> probably will have lifetime issues then
//TODO: Also, owning point cloud data means we can sort it to increase cache locality (but then I need to keep track of correct indices/references in `Cell`)
//TODO: I really see now why it's supposed to be a cell *list*. But it would be nice to get the same behaviour with more contiguous memory for cache locality
//TODO: but maybe I'm trying to optimize prematurely
struct CellList {
    point_indices: Vec<usize>, // alternatives: LinkedList, VecDeque?
}

//TODO: do I want to store a flat index or is it enough to build an iterator?
/*struct MultiIndex {
    //storage: ndarray::Array, //tbd.
    grid_info: usize, //GridInfo,
}*/

//TODO: should this include the lattice type (mapping coordinates -> lattice points)
/*struct GridInfo {
    todo!(),
}*/
*/
struct PointCloud(Vec<Point3<f64>>);
struct AABB {
    min: Point3<f64>,
    max: Point3<f64>,
}
struct GridInfo {
    cutoff: f64,
    //shape: (usize, usize, usize), //TODO: make sure I provide constructors but do I want to store both for convenience?
    aabb: AABB,
}
struct MultiIndex(Vec<[usize; 3]>);

impl MultiIndex {
    from_points(points: &[[f64; 3]], cutoff: f64)
}

/// This is a naive structure of cell lists following Allen and Tildesley
/// We're not yet even using type aliases.
//TODO: would be nice to not have to know the number of cells (or bounding box dimensions and cutoff) beforehand
struct NaiveCellGrid {
    point_cloud: Vec<Point3<f64>>, // this is array `r` in Allen and Tildesley
    //TODO: Make sure we construct multi_index column-major (F order)
    //TODO: or store it transposed
    //TODO: or just use Array1/Vec<[usize; 3]>? It's the triples that we want anyway. These are fixed size, i.e. are stored contigouosly
    multi_index: Array2<usize>, // this is array `c` in Allen and Tildesley
    // We could make `cells: Array1<usize>` if we construct a flattened index from `multi_index`
    // but ndarray already makes sure Array3 is contiguous so we don't need manual strides.
    // Instead, we just map from `multi_index` to `cells` to `cell_list`
    // `None` means this cell is empty
    cells: Array3<Option<usize>>, // this is array `head` in Allen and Tildesley
    // `None` means the end of a cell list is reached
    //TODO: perhaps we could use specialized sum types instead of `Option<>` to store info about next non-empty cell?
    cell_list: Vec<Option<usize>>, // this is array `list` in Allen and Tildesley
}

impl NaiveCellGrid {
    fn new(points: &[[f64; 3]], cutoff: f64) -> Self {
        let point_cloud: Vec<Point3<f64>> = points.iter().map(|p| Point3::from(*p)).collect();
        let n = point_cloud.len();

        //TODO: would be nice to get the smallest bounding box around point cloud for multi-index computation?
        //TODO: leverage nalgebra's point type to compute inf/sup to get the AABB for the multi index
        //TODO: right now I'm mostly concerned with a sensible API (i.e. what type contains which behaviour and how should I plug them together)
        let cells = Array3::default(shape)

        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naive_cellgrid() {
        let points = vec![
            [72.46735284035672, 46.44550378090424, 66.63901205139712],
            [57.69674557034732, 54.53844814921267, 69.0970244304932],
            [34.67813071756326, 53.16888452609834, 61.616297234615],
            [40.13882919260389, 28.975922417310834, 55.253885301067996],
            [42.61828804027062, 35.20552057284047, 62.1781605010535],
            [63.14255785611127, 46.75634193050543, 17.216733004395778],
            [64.58535798025872, 43.924718999414885, 44.359472356623854],
            [56.508264823856265, 51.74535180326326, 88.07572874002203],
            [51.40591695755581, 75.62211010510869, 47.55107384210809],
            [44.08156272328777, 51.16721452838641, 52.49586501282618],
            [34.38034961466383, 35.35377362779141, 54.70795827372833],
            [58.079935039320354, 55.00772926184105, 54.68504356973868],
            [52.87188519819982, 46.67167992644182, 67.34078941111173],
            [38.17111850725933, 28.113313500884335, 53.987122308925414],
            [59.037413467130236, 50.55283549816238, 64.76182846629645],
            [43.12807740901327, 57.642793461405276, 57.89444992914932],
            [47.75499500945548, 54.99327943953817, 38.231867403361356],
            [25.188112076540705, 57.80348478947813, 41.94776503979341],
            [61.26791732496439, 50.3328612531578, 67.28236273479564],
            [46.31639982669126, 34.809947204968545, 55.41687957872044],
            [68.24489572989566, 64.36252649586478, 45.59532638768112],
            [39.642539405028884, 46.55482138885022, 20.66386840274636],
            [46.42323385809075, 44.84961716169285, 62.03938038350073],
            [63.83198454312149, 23.08193780980907, 60.85876845877909],
            [47.49313247204821, 40.0130740435201, 39.48034890954136],
            [63.70990546387425, 49.6834719629841, 38.76832383907555],
            [73.22976567347294, 34.38492771863741, 63.10497561401529],
            [31.56303787081657, 37.2810019672349, 49.186319948269805],
            [48.891914222109676, 45.908605995610905, 57.85909914298353],
            [42.19940343638644, 62.05045024256016, 46.04615587474219],
            [61.36836663969868, 63.610202615133026, 65.50585643407233],
            [53.22281550411799, 38.915573426336636, 56.142309179608084],
            [97.22990511601748, 54.21499759713862, 56.8923471572009],
            [44.45100277301918, 79.25799136076189, 63.8605867569631],
            [65.48588740173099, 33.618136876932795, 58.813801883008665],
            [59.84754006539127, 39.91580402346271, 27.31424273259812],
            [10.220666496210917, 67.15814360954151, 33.83605442858338],
            [49.86390666339886, 49.83180624543126, 66.48642667266306],
            [30.658738700695373, 30.13233834816357, 29.77580361060121],
            [40.55874707553164, 65.78441682847034, 65.329704096891],
            [78.34600396932015, 57.89344459966851, 63.147527973823806],
            [32.80493512781791, 83.82047561851411, 26.347874951756527],
            [64.21155331024228, 51.55767210116186, 44.488493919235765],
            [29.774253608613687, 60.36683123467222, 48.98986903681975],
            [61.9463615693191, 46.70151655296417, 40.60905723084488],
            [57.25175222727173, 58.267578504791096, 45.627601887896226],
            [48.990436706721155, 39.28691192908618, 65.99996069071857],
            [34.9551652667681, 33.993535834480646, 56.857243146927296],
            [31.552421500401945, 44.69629084263187, 36.978637735334694],
            [37.969456344079475, 56.73920865114978, 43.723353504599494],
            [100.7362316586383, 29.325475738352132, 39.2961631010229],
            [37.59217400693751, 54.6377588907565, 64.115798853594],
            [56.737265860640804, 39.04376659614116, 52.94883384952191],
            [61.138496819067534, 63.2680559138494, 51.68343302237743],
            [36.178791209183004, 37.07869560366909, 89.81510004537876],
            [25.367500334103767, 37.20404158565769, 71.20194533668426],
            [39.5658971262939, 37.71015668463188, 69.69175059087647],
            [74.59494233788367, 58.9386803091584, 49.89875712671498],
            [36.027919768941864, 26.709842713120487, 55.55531436350241],
            [60.19951310819511, 50.43685316547253, 48.84434124863185],
            [47.08445998269379, 51.42256580514277, 34.77265928236812],
            [16.239633312593277, 68.03538296971537, 45.690378394764245],
            [44.146625228728404, 71.5587709356264, 59.59311627038948],
            [67.73847715205741, 69.47789890147548, 72.04517453522703],
            [46.43867689855991, 65.93355978193483, 72.17843732119677],
            [64.76470527538888, 38.097835284089044, 23.65964644922338],
            [22.00615978111906, 76.67598552439202, 29.24252006289259],
            [64.08791191993704, 53.09100419914532, 48.754861682238875],
            [55.46117771890817, 24.53825031550556, 45.18370209343301],
            [86.07382920040648, 43.237350300817866, 55.591644162868114],
            [54.83014612842436, 46.91667248552966, 76.70889773887879],
            [61.67492326946063, 30.54211202203747, 41.95249605247814],
            [21.611653936078874, 40.489755209396264, 78.65322598488636],
            [27.253251225145146, 66.46042980543669, 62.14726515202773],
            [29.5965367421379, 58.48553446136442, 19.986911598303685],
            [47.108958973879545, 42.01956331048258, 34.129778053404785],
            [73.17808321663728, 57.65297654476593, 30.211840054197957],
            [13.667308188070699, 43.56974240978837, 61.20601784396419],
            [31.84382087562678, 47.226355580202274, 34.865970681666624],
            [31.547392451104017, 49.33602923976242, 33.12473753947354],
            [39.68776738777069, 23.089387957833956, 62.45456486865129],
            [58.86878521076362, 52.42993012186303, 30.29545978336863],
            [37.47610734759132, 29.434093661037558, 43.32423412815852],
            [58.57365118198644, 61.58027783599883, 61.411711315063194],
            [54.15337178381058, 59.34935299268743, 61.84440491486255],
            [62.19631196240582, 46.14088751172073, 57.496210851682065],
            [57.163914610126795, 35.61326734970044, 42.95106224503242],
            [24.73493184392569, 40.653528344244116, 44.673829547322676],
            [53.33466064076672, 2.1362788844113467, 58.997439366269994],
            [36.47450202576266, 56.04090352957959, 25.22421829800911],
            [52.9409340230885, 18.71357687358318, 66.62607832586342],
            [83.3291701136857, 54.40985511139716, 9.32770492411823],
            [42.89771938036293, 43.96081062358374, 44.69028713355895],
            [24.970927379898757, 39.22858982829084, 58.22650711255713],
            [84.76794322987249, 23.223616574963835, 50.20197916194471],
            [29.242059940317862, 46.68577310455115, 62.99451765905497],
            [40.39908569084808, 34.73316890811478, 45.49875863568221],
            [54.50362129284737, 39.2272313514261, 65.1884692955773],
            [44.07093095603651, 33.768156172870206, 54.373898951277155],
            [39.929742543788265, 62.62236865041245, 51.624880489282006],
        ];
        todo!()
    }
}

