# Rustbustools

Rust version of [bustools](https://github.com/BUStools/bustools) command line interface.
At this point, it's **far from complete and correct**, but rather a project to learn rust.

This is heavily built on [rustbustools](https://github.com/redst4r/rustbustools), which handles all the basic interactions with busfiles.
## Examples
For more examples, see the rust-docs.

### Iterating a bus file
```rust
# use rustbustools::io::{BusReader};
let bus = BusReader::new("/tmp/some.bus");
for record in bus {
    // record.CB, record.UMI ...
}
```

### Iterating a bus file by cell
```rust
use rustbustools::io::BusReader;
use rustbustools::iterators::CellGroupIterator; //need to bring that trait into scope

let breader = BusReader::new("/path/to/some.bus");
for (cb, vector_of_records) in breader.groupby_cb() {
    // Example: the number of records in that cell
    let n_molecules: usize = vector_of_records.len();
}
```