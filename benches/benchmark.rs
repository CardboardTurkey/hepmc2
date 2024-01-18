#![cfg(feature = "sync")]

use std::convert::{AsRef, From};
use std::default::Default;
use std::f64::consts::PI;
use std::io::BufReader;

use criterion::{criterion_group, criterion_main, Criterion};
use hepmc2::{Reader, Writer};
use rand::distributions::{Alphanumeric, Distribution, Standard};
use rand::{Rng, SeedableRng};

struct Event(hepmc2::event::Event);

impl AsRef<hepmc2::event::Event> for Event {
    fn as_ref(&self) -> &hepmc2::event::Event {
        &self.0
    }
}

impl Distribution<Event> for Standard {
    fn sample<R: Rng + ?Sized>(&self, mut rng: &mut R) -> Event {
        Event(hepmc2::event::Event {
            alpha_qcd: rng.gen_range(0.1..0.12),
            alpha_qed: 1. / 137.,
            energy_unit: Default::default(),
            length_unit: Default::default(),
            mpi: rng.gen(),
            number: rng.gen(),
            pdf_info: rng.gen::<PdfInfo>().into(),
            random_states: {
                let len = rng.gen_range(0..4);
                (0..len).map(|_| rng.gen()).collect()
            },
            scale: rng.gen(),
            signal_process_id: rng.gen(),
            signal_process_vertex: rng.gen(),
            vertices: {
                let len = rng.gen_range(1..6);
                (0..len).map(|_| rng.gen::<Vertex>().into()).collect()
            },
            weight_names: {
                let len = rng.gen_range(0..11);
                (0..len).map(|_| gen_random_name(&mut rng)).collect()
            },
            weights: {
                let len = rng.gen_range(0..11);
                (0..len).map(|_| rng.gen()).collect()
            },
            xs: rng.gen::<CrossSection>().into(),
            heavy_ion_info: None,
        })
    }
}

struct PdfInfo(hepmc2::event::PdfInfo);

impl From<PdfInfo> for hepmc2::event::PdfInfo {
    fn from(p: PdfInfo) -> Self {
        p.0
    }
}

impl Distribution<PdfInfo> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> PdfInfo {
        PdfInfo(hepmc2::event::PdfInfo {
            parton_id: rng.gen(),
            pdf_id: rng.gen(),
            scale: rng.gen(),
            x: rng.gen(),
            xf: rng.gen(),
        })
    }
}

struct Vertex(hepmc2::event::Vertex);

impl From<Vertex> for hepmc2::event::Vertex {
    fn from(v: Vertex) -> Self {
        v.0
    }
}

impl Distribution<Vertex> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vertex {
        Vertex(hepmc2::event::Vertex {
            barcode: 0,
            particles_in: {
                let len = rng.gen_range(0..=2);
                (0..len).map(|_| rng.gen::<Particle>().into()).collect()
            },
            particles_out: {
                let len = rng.gen_range(0..=5);
                (0..len).map(|_| rng.gen::<Particle>().into()).collect()
            },
            status: rng.gen(),
            t: rng.gen(),
            weights: {
                let len = rng.gen_range(0..3);
                (0..len).map(|_| rng.gen()).collect()
            },
            x: rng.gen(),
            y: rng.gen(),
            z: rng.gen(),
        })
    }
}

struct Particle(hepmc2::event::Particle);

impl From<Particle> for hepmc2::event::Particle {
    fn from(v: Particle) -> Self {
        v.0
    }
}

impl Distribution<Particle> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Particle {
        Particle(hepmc2::event::Particle {
            end_vtx: 0,
            flows: Default::default(),
            id: rng.gen_range(-30..30),
            m: rng.gen_range(0.0..175.),
            p: rng.gen::<FourVector>().into(),
            phi: rng.gen_range(-PI..PI),
            status: rng.gen_range(-1..3),
            theta: rng.gen_range(0.0..PI),
        })
    }
}

struct FourVector(hepmc2::event::FourVector);

impl From<FourVector> for hepmc2::event::FourVector {
    fn from(v: FourVector) -> Self {
        v.0
    }
}

impl Distribution<FourVector> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> FourVector {
        FourVector(hepmc2::event::FourVector(rng.gen()))
    }
}

struct CrossSection(hepmc2::event::CrossSection);

impl From<CrossSection> for hepmc2::event::CrossSection {
    fn from(xs: CrossSection) -> Self {
        xs.0
    }
}

impl Distribution<CrossSection> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> CrossSection {
        CrossSection(hepmc2::event::CrossSection {
            cross_section: rng.gen(),
            cross_section_error: rng.gen(),
        })
    }
}

fn gen_random_name<R: Rng>(rng: &mut R) -> String {
    let len = rng.gen_range(0..6);
    rng.sample_iter(&Alphanumeric)
        .take(len)
        .map(char::from)
        .collect()
}

const NEVENTS: usize = 3_000;

fn criterion_benchmark(c: &mut Criterion) {
    let mut buf: Vec<u8> = Vec::new();

    {
        let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(0);
        let events: Vec<Event> = (0..NEVENTS).map(|_| rng.gen()).collect();
        c.bench_function("write", |b| {
            b.iter(|| {
                buf.clear();
                let mut writer = Writer::new(&mut buf).unwrap();
                for event in &events {
                    writer.write(event.as_ref()).unwrap()
                }
                writer.finish().unwrap()
            })
        });
    }

    c.bench_function("read", |b| {
        b.iter(|| {
            let mut count = 0;
            let buf = BufReader::new(buf.as_slice());
            let reader = Reader::new(buf);
            for _event in reader {
                count += 1
            }
            assert_eq!(count, NEVENTS)
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
