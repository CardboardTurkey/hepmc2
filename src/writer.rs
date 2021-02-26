use std::default::Default;
use std::fmt::Display;
use std::io::{self, Write};
use std::mem::take;

use crate::event::*;

use log::error;

const DEFAULT_HEADER: &str = "HepMC::Version 2.06.09
HepMC::IO_GenEvent-START_EVENT_LISTING
";

const DEFAULT_FOOTER: &[u8] = b"HepMC::IO_GenEvent-END_EVENT_LISTING\n";

/// Writer for the HepMC2 format
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
pub struct Writer<T: Write> {
    stream: T,
    finished: bool,
}

impl<T: Write + Default> Writer<T> {
    pub fn new() -> Self {
        Self::default()
    }

    /// Retrieve the underlying writer
    pub fn into_inner(mut self) -> T {
        // ensure that the destructor doesn't do anything
        self.finished = true;
        take(&mut self.stream)
    }
}

impl<T: Write> Writer<T> {
    pub fn finish(mut self) -> Result<(), std::io::Error> {
        self.ref_finish()
    }

    fn ref_finish(&mut self) -> Result<(), std::io::Error> {
        self.stream.write_all(DEFAULT_FOOTER)?;
        self.finished = true;
        Ok(())
    }

    fn write_header<U: Display>(&mut self, header: U) -> Result<(), io::Error> {
        write!(self.stream, "{}", header)
    }

    pub fn try_from(stream: T) -> Result<Self, io::Error> {
        Self::with_header(stream, DEFAULT_HEADER)
    }

    pub fn with_header<U: Display>(
        stream: T,
        header: U,
    ) -> Result<Self, io::Error> {
        let mut writer = Self { stream, finished: false };
        writer.write_header(header)?;
        Ok(writer)
    }

    pub fn write(&mut self, event: &Event) -> Result<(), io::Error> {
        self.write_event_line(&event)?;
        if !event.weight_names.is_empty() {
            self.write_weight_names_line(&event.weight_names)?;
        }
        self.write_unit_line(&event)?;
        self.write_cross_section_line(&event.xs)?;
        self.write_pdf_info_line(&event.pdf_info)?;
        for vertex in &event.vertices {
            self.write_vertex_line(vertex)?;
            let particles = vertex
                .particles_in
                .iter()
                .chain(vertex.particles_out.iter());
            for particle in particles {
                self.write_particle_line(particle)?;
            }
        }
        Ok(())
    }

    fn write_event_line(&mut self, event: &Event) -> Result<(), io::Error> {
        write!(
            self.stream,
            "E {} {} {} {} {} {} {} {} 0 0 {}",
            event.number,
            event.mpi,
            ryu::Buffer::new().format(event.scale),
            ryu::Buffer::new().format(event.alpha_qcd),
            ryu::Buffer::new().format(event.alpha_qed),
            event.signal_process_id,
            event.signal_process_vertex,
            event.vertices.len(),
            event.random_states.len()
        )?;
        for state in &event.random_states {
            write!(self.stream, " {}", state)?;
        }
        write!(self.stream, " {}", event.weights.len())?;
        let mut buffer = ryu::Buffer::new();
        for weight in &event.weights {
            write!(self.stream, " {}", buffer.format(*weight))?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_vertex_line(&mut self, vertex: &Vertex) -> Result<(), io::Error> {
        write!(
            self.stream,
            "V {} {} {} {} {} {} 0 {} {}",
            vertex.barcode,
            vertex.status,
            ryu::Buffer::new().format(vertex.x),
            ryu::Buffer::new().format(vertex.y),
            ryu::Buffer::new().format(vertex.z),
            ryu::Buffer::new().format(vertex.t),
            vertex.particles_in.len() + vertex.particles_out.len(),
            vertex.weights.len()
        )?;
        for weight in &vertex.weights {
            write!(self.stream, " {}", weight)?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_particle_line(
        &mut self,
        particle: &Particle,
    ) -> Result<(), io::Error> {
        write!(
            self.stream,
            "P 0 {} {} {} {} {} {} {} {} {} {} {}",
            particle.id,
            ryu::Buffer::new().format(particle.p[1]),
            ryu::Buffer::new().format(particle.p[2]),
            ryu::Buffer::new().format(particle.p[3]),
            ryu::Buffer::new().format(particle.p[0]),
            ryu::Buffer::new().format(particle.m),
            particle.status,
            ryu::Buffer::new().format(particle.theta),
            ryu::Buffer::new().format(particle.phi),
            particle.end_vtx,
            particle.flows.len()
        )?;
        for (idx, val) in &particle.flows {
            write!(self.stream, " {} {}", idx, val)?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_weight_names_line(
        &mut self,
        names: &[String],
    ) -> Result<(), io::Error> {
        write!(self.stream, "N {}", names.len())?;
        for name in names {
            write!(self.stream, r#" "{}""#, name)?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_unit_line(&mut self, event: &Event) -> Result<(), io::Error> {
        if event.energy_unit.is_empty() || event.length_unit.is_empty() {
            return Ok(());
        }
        writeln!(self.stream, "U {} {}", event.energy_unit, event.length_unit)
    }

    fn write_cross_section_line(
        &mut self,
        xs: &CrossSection,
    ) -> Result<(), io::Error> {
        writeln!(
            self.stream,
            "C {} {}",
            ryu::Buffer::new().format(xs.cross_section),
            ryu::Buffer::new().format(xs.cross_section_error)
        )
    }

    fn write_pdf_info_line(&mut self, pdf: &PdfInfo) -> Result<(), io::Error> {
        writeln!(
            self.stream,
            "F {} {} {} {} {} {} {} {} {}",
            pdf.parton_id[0],
            pdf.parton_id[1],
            ryu::Buffer::new().format(pdf.x[0]),
            ryu::Buffer::new().format(pdf.x[1]),
            ryu::Buffer::new().format(pdf.scale),
            ryu::Buffer::new().format(pdf.xf[0]),
            ryu::Buffer::new().format(pdf.xf[1]),
            pdf.pdf_id[0],
            pdf.pdf_id[1],
        )
    }
}

impl<T: Write> Drop for Writer<T> {
    fn drop(&mut self) {
        if !self.finished {
            error!("Hepmc2 writer dropped before finished.");
            error!("Call finish() manually to fix this error.");
            self.ref_finish().unwrap()
        }
    }
}
