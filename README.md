# Rustbustools

Rust library to interact with the kallisto/bus format of scRNAseq data (see [bustools](https://github.com/BUStools/bustools)).
At this point, the package is pretty mature, but there might be some minor features missing compared to the original bustools.

There's also a CLI mimicking [bustools](https://github.com/BUStools/bustools), see [bustools_cli](https://github.com/redst4r/bustools_cli-rs)
## Examples
For more examples, see the rust-docs.

### Iterating a bus file
```rust
use bustools::io::{BusReader};
let bus = BusReader::new("/tmp/some.bus");
for record in bus {
    // record.CB, record.UMI ...
}
```

### Iterating a bus file by cell
```rust
use bustools::io::BusReader;
use bustools::iterators::CellGroupIterator; //need to bring that trait into scope

let breader = BusReader::new("/path/to/some.bus");
for (cb, vector_of_records) in breader.groupby_cb() {
    // Example: the number of records in that cell
    let n_molecules: usize = vector_of_records.len();
}
```
