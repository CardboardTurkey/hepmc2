# hepmc2

Read and write event files in the `hepmc2` format, also known as
`IO_GenEvent`.

## Caveats

This crate is inspired by the code for the `ReaderAsciiHepMC2` in the
[HepMC3 library](https://gitlab.cern.ch/hepmc/HepMC3), version
3.2.0. The aim is to be fully compatible, but be aware that the
current tests are not exhaustive.

## Example

```rust
// Read events from `events_in.hepmc2` and write them to `events_out.hepmc2`
use hepmc2::{Reader, Writer};

use std::io::BufReader;
use std::fs::File;

let input = BufReader::new(File::open("events_in.hepmc2")?);
let in_events = Reader::from(input);

let output = File::create("events_out.hepmc2")?;
let mut writer = Writer::try_from(output)?;

for event in in_events {
   let event = event?;
   println!("Current cross section: {}",  event.xs);
   writer.write(&event)?
}
writer.finish()?;
```

## Async API

By default this crate enables the `sync` feature which exposes a sync API. You
can however switch to using a `tokio`-backed async API by disabling the `sync`
feature and enabling the `tokio` feature.

Either run the following in the root of your crate:

```sh
cargo add hepmc2 --no-default-features -F tokio
```

or make sure a line like the following is present in your `Cargo.toml`:

```toml
hepmc2 = { version = "0.6.0", default-features = false, features = ["tokio"] }
```

The async API is exactly the same as the sync one but IO operations will return
futures that you will, as usual, need to call `await` on. For examples, generate
the async API documentation in the root of this project:

```sh
cargo doc --open --no-default-features -F tokio
```

License: GPL-3.0-or-later
