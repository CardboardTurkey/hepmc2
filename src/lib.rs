//! Read and write event files in the `hepmc2` format, also known as
//! `IO_GenEvent`.
//!
//! # Caveats
//!
//! This crate is inspired by the code for the `ReaderAsciiHepMC2` in the
//! [HepMC3 library](https://gitlab.cern.ch/hepmc/HepMC3), version
//! 3.2.0. The aim is to be fully compatible, but be aware that the
//! current tests are not exhaustive.
//!
//! # Example
//!
//! ```rust,no_run
//! // Read events from `events_in.hepmc2` and write them to `events_out.hepmc2`
//! use hepmc2::{Reader, Writer};
//!
//! use std::io::BufReader;
//! use std::fs::File;
//!
//! let input = BufReader::new(File::open("events_in.hepmc2")?);
//! let in_events = Reader::from(input);
//!
//! let output = File::create("events_out.hepmc2")?;
//! let mut writer = Writer::try_from(output)?;
//!
//! for event in in_events {
//!    let event = event?;
//!    println!("Current cross section: {}",  event.xs);
//!    writer.write(&event)?
//! }
//! writer.finish()?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
pub mod event;
pub mod reader;
pub mod writer;

pub use crate::event::Event;
pub use crate::reader::Reader;
pub use crate::writer::Writer;

#[cfg(all(feature = "sync", feature = "tokio"))]
compile_error!("One and only one sync/async feature must be enabled");
#[cfg(not(any(feature = "sync", feature = "tokio")))]
compile_error!("One and only one sync/async feature must be enabled");

#[cfg(test)]
mod tests {
    use super::*;

    #[maybe_async::test(
        feature = "sync",
        async(
            all(not(feature = "sync"), feature = "tokio"),
            tokio::test(flavor = "multi_thread")
        )
    )]
    async fn tst_read() {
        let mut reader = reader::Reader::from(EVENT_TXT);
        let next_line = reader.next().await.unwrap();
        assert!(next_line.is_ok());
        let next_line = reader.next().await;
        assert!(next_line.is_none());
    }

    #[maybe_async::test(
        feature = "sync",
        async(
            all(not(feature = "sync"), feature = "tokio"),
            tokio::test(flavor = "multi_thread")
        )
    )]
    async fn tst_read_write() {
        #[cfg(feature = "sync")]
        use std::io::BufReader;
        #[cfg(feature = "tokio")]
        use tokio::io::BufReader;

        let mut reader = reader::Reader::from(EVENT_TXT);
        let mut buf = Vec::<u8>::new();
        let event = reader.next().await.unwrap().unwrap();
        {
            let mut writer = writer::Writer::try_from(&mut buf).await.unwrap();
            writer.write(&event).await.unwrap();
        }
        let mut reader = reader::Reader::from(BufReader::new(buf.as_slice()));
        let event = reader.next().await.unwrap().unwrap();
        let mut buf2 = Vec::<u8>::new();
        let mut writer = writer::Writer::try_from(&mut buf2).await.unwrap();
        writer.write(&event).await.unwrap();
        writer.finish().await.unwrap();
        use std::str::from_utf8;
        assert_eq!(from_utf8(&buf), from_utf8(&buf2));
    }

    const EVENT_TXT: &[u8] = br#"
HepMC::Version 2.06.09
HepMC::IO_GenEvent-START_EVENT_LISTING
E 0 -1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 0 0 23 1 2 0 16 5.5606031127834702e-09 5.5606031127834702e-09 0 -1.0000000000000000e+00 2.6206847299999998e+00 5.8670509300000004e+00 6.1711295699999997e+00 -1.0000000000000000e+00 5.9781127600000001e+00 3.9987750200000001e+00 2.8889888799999999e+00 5.5606031100000000e-09 -1.0000000000000000e+00 8.6938853699999996e-01 3.5507575000000002e-01 5.6137609399999999e-01
N 16 "0" "Weight" "eventNumber" "phi" "phi1" "phi2" "phi3" "t" "t1" "t2" "t3" "unlops_weight" "z" "z1" "z2" "z3"
U GEV MM
C 5.5606031127834701e+00 5.3451183000000005e+04
V -1 0 0 0 0 0 0 1 0
P 3 -1 0 0 5.2051533588697652e+01 5.2051533588697652e+01 0 21 0 0 -3 1 2 501
V -2 0 0 0 0 0 0 2 0
P 4 2 0 0 -2.9818545945620773e+01 2.9818545945620773e+01 0 21 0 0 -3 1 1 501
P 9 21 5.7076480059072443e+00 3.2748659206013784e+00 -2.4249351062062976e+00 7.0130095413142710e+00 0 43 0 0 -8 2 1 502 2 501
V -3 0 0 0 0 0 0 1 0
P 5 24 -4.6339895476421589e-08 -1.0408832373798305e-07 2.2232987710669342e+01 8.1870080240182716e+01 7.8793428000000006e+01 22 0 0 -6 0
V -4 0 0 0 0 0 1 1 0
P 1 2212 0 0 6.9999999371178146e+03 7.0000000000000000e+03 9.3827000000000005e-01 4 0 0 -4 0
P 6 -1 0 0 5.2051533588697652e+01 5.2051533588697652e+01 0 42 0 0 -1 1 2 501
V -5 0 0 0 0 0 0 2 0
P 7 2 0 0 -3.6129848231502301e+01 3.6129848231502301e+01 0 41 0 0 -2 1 1 502
P 12 -2 2.5237885213344935e+00 3.3116569526639856e+00 -6.9143141497758080e+01 6.9268395365069210e+01 0 1 0 0 0 1 2 503
V -6 0 0 0 0 0 0 1 0
P 8 24 -5.7076481044982819e+00 -3.2748660546697370e+00 1.8346620499347921e+01 8.1168372992804763e+01 7.8793428000000006e+01 44 0 0 -23 0
V -7 0 0 0 0 0 0 2 0
P 10 21 0 0 -1.0175250837889556e+02 1.0175250837889556e+02 0 41 0 0 -5 2 1 502 2 503
P 15 21 2.9123097878176076e+00 2.8941914776608224e-01 -8.0591874171438874e+01 8.0644996709085106e+01 0 43 0 0 -18 2 1 503 2 504
V -8 0 0 0 0 0 0 1 0
P 11 21 3.1838594845727508e+00 -3.6791032062607254e-02 1.0955462441585269e+00 3.3672743236383127e+00 0 48 0 0 -10 2 1 502 2 501
V -9 0 0 0 0 0 0 2 0
P 13 21 0 0 -1.8125550661114170e+02 1.8125550661114170e+02 0 41 0 0 -7 2 1 502 2 504
P 18 2 -6.0739028889597069e-01 2.5046574652089002e+00 -2.2999284645202547e+03 2.2999299085281878e+03 0 43 0 0 -17 1 1 504
V -10 0 0 0 0 0 0 1 0
P 14 21 2.7154969675514318e-01 -3.2621017982868949e-01 2.1844221833512592e+00 2.2252758467993488e+00 0 48 0 0 -12 2 1 502 2 501
V -11 0 0 0 0 0 0 2 0
P 16 2 0 0 -2.4821608453435269e+03 2.4821608453435274e+03 0 41 0 0 -9 1 1 502
P 21 21 -1.5176954196303214e+00 -1.1581464938267174e+00 -2.1424758388954547e+02 2.1425608954052109e+02 0 43 0 0 -16 2 1 505 2 502
V -12 0 0 0 0 0 0 1 0
P 17 21 8.7893998565111386e-01 -2.8308676450375896e+00 1.2075479712207198e+00 3.2007060509972689e+00 0 48 0 0 -14 2 1 502 2 501
V -13 0 0 0 0 0 1 1 0
P 2 2212 0 0 -6.9999999371178146e+03 7.0000000000000000e+03 9.3827000000000005e-01 4 0 0 -13 0
P 19 2 0 0 -2.6963868241367072e+03 2.6963868241367081e+03 0 41 0 0 -11 1 1 505
V -14 0 0 0 0 0 0 1 0
P 20 21 2.3966354052814354e+00 -1.6727211512108722e+00 1.2291530675859121e+00 3.1705953036569099e+00 0 48 0 0 -15 2 1 502 2 501
V -15 0 0 0 0 0 0 2 0
P 22 21 1.6780099262880444e+00 -1.2585715889494420e-02 1.5618589221168571e+00 2.2924395315359480e+00 0 1 0 0 0 2 1 506 2 501
P 23 21 7.0154368522484822e-01 -1.6731704747007090e+00 -2.7440810198086751e+00 3.2896266692386575e+00 0 1 0 0 0 2 1 502 2 506
V -16 0 0 0 0 0 0 1 0
P 24 21 -1.5006136258617786e+00 -1.1451114544473862e+00 -2.1183620872426775e+02 2.1184461864340338e+02 0 1 0 0 0 2 1 505 2 502
V -17 0 0 0 0 0 0 2 0
P 25 2 -1.6033704065363041e+00 4.3863302327067144e+00 -2.1547456957401910e+03 2.1547507568163073e+03 0 51 0 0 -20 1 1 507
P 26 21 2.6011093609332741e+00 -1.7221584398720926e+00 -1.8960124875750688e+02 1.8962691035295006e+02 0 51 0 0 -19 2 1 504 2 507
V -18 0 0 0 0 0 0 1 0
P 27 21 1.3071805445246669e+00 1.2990482014036053e-01 -3.6173394193995875e+01 3.6197238068015622e+01 0 1 0 0 0 2 1 503 2 504
V -19 0 0 0 0 0 0 2 0
P 28 21 3.9678447216374000e+00 -1.0579393283155296e+00 -3.5312021479758994e+02 3.5314409116683936e+02 0 51 0 0 -22 2 1 504 2 508
P 29 21 -1.9866544379963598e+00 1.0316895689562262e+00 -6.6958108247799419e+02 6.6958482449322673e+02 0 51 0 0 -21 2 1 508 2 507
V -20 0 0 0 0 0 0 1 0
P 30 2 -9.8345132924407030e-01 2.6904215521939250e+00 -1.3216456472221134e+03 1.3216487515091915e+03 0 1 0 0 0 1 1 507
V -21 0 0 0 0 0 0 2 0
P 31 21 -9.8001674572977349e-01 1.3322582893034975e+00 -6.1464412971562535e+02 6.1464635485688007e+02 0 1 0 0 0 2 1 509 2 507
P 32 21 -6.7009385001849431e-01 -3.9030080197012751e-01 -8.4887830929292363e+01 8.4891372944274849e+01 0 1 0 0 0 2 1 508 2 509
V -22 0 0 0 0 0 0 1 0
P 33 21 3.6313008793893080e+00 -9.6820724669267344e-01 -3.2316933663066641e+02 3.2319118785891123e+02 0 1 0 0 0 2 1 504 2 508
V -23 0 0 0 0 0 0 2 0
P 34 -13 6.7280635852793118e+00 2.6135691981490424e+01 3.5757112587329956e+01 4.4798699749357091e+01 1.0565837000000000e-01 1 0 0 0 0
P 35 14 -1.2435711898718612e+01 -2.9410557911101773e+01 -1.7410492121416357e+01 3.6369672571802241e+01 0 1 0 0 0 0
"#;
}
