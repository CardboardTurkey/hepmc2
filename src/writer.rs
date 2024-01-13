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
    /// Retrieve the underlying writer
    pub fn into_inner(mut self) -> T {
        // ensure that the destructor doesn't do anything
        self.finished = true;
        take(&mut self.stream)
    }
}

impl<T: Write> Writer<T> {
    /// Construct new `Writer`
    ///
    /// This automatically tries to write the mandatory HepMC header,
    /// which may fail.
    ///
    /// # Example
    ///
    /// ```rust
    /// use hepmc2::writer::Writer;
    ///
    /// let mut output = Vec::new();
    /// let mut writer = Writer::new(&mut output)?;
    /// // always call finish at the end
    /// writer.finish()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(stream: T) -> Result<Self, io::Error> {
        Self::with_header(stream, DEFAULT_HEADER)
    }

    /// Construct new `Writer`, trying to write a custom header
    ///
    /// `hepmc2` ignores headers, but other implementations of the
    /// format may be less lenient.
    ///
    /// # Example
    ///
    /// ```rust
    /// use hepmc2::writer::Writer;
    ///
    /// let mut output = Vec::new();
    /// let mut writer = Writer::with_header(output, "")?;
    /// // always call finish at the end
    /// writer.finish()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn with_header<U: Display>(stream: T, header: U) -> Result<Self, io::Error> {
        let mut writer = Self {
            stream,
            finished: false,
        };
        writer.write_header(header)?;
        Ok(writer)
    }

    /// Finish writing, consuming the `Writer`
    ///
    /// This tries to write the mandatory HepMC footer, which may fail.
    ///
    /// # Example
    ///
    /// ```rust
    /// use hepmc2::writer::Writer;
    ///
    /// let mut output = Vec::new();
    /// let mut writer = Writer::new(&mut output)?;
    /// // always call finish at the end
    /// writer.finish()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn finish(mut self) -> Result<(), std::io::Error> {
        self.ref_finish()
    }

    /// Write an event
    ///
    /// # Example
    ///
    /// ```rust
    /// use hepmc2::writer::Writer;
    /// use hepmc2::event::Event;
    ///
    /// let mut output = Vec::new();
    /// let mut writer = Writer::new(&mut output)?;
    /// let event = Event::default();
    /// writer.write(&event)?;
    /// // always call finish at the end
    /// writer.finish()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write(&mut self, event: &Event) -> Result<(), io::Error> {
        self.write_event_line(event)?;
        if !event.weight_names.is_empty() {
            self.write_weight_names_line(&event.weight_names)?;
        }
        self.write_unit_line(event)?;
        self.write_cross_section_line(&event.xs)?;
        self.write_pdf_info_line(&event.pdf_info)?;
        if let Some(hi) = event.heavy_ion_info {
            self.write_heavy_ion_info_line(&hi)?;
        }
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

    pub fn try_from(stream: T) -> Result<Self, io::Error> {
        Self::with_header(stream, DEFAULT_HEADER)
    }

    fn ref_finish(&mut self) -> Result<(), std::io::Error> {
        self.stream.write_all(DEFAULT_FOOTER)?;
        self.finished = true;
        Ok(())
    }

    fn write_header<U: Display>(&mut self, header: U) -> Result<(), io::Error> {
        self.stream.write_all(header.to_string().as_bytes())
    }

    fn write_event_line(&mut self, event: &Event) -> Result<(), io::Error> {
        self.stream.write_all(
            format!(
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
            )
            .as_bytes(),
        )?;
        for state in &event.random_states {
            self.stream.write_all(format!(" {state}").as_bytes())?
        }
        self.stream
            .write_all(format!(" {}", event.weights.len()).as_bytes())?;
        let mut buffer = ryu::Buffer::new();
        for weight in &event.weights {
            self.stream
                .write_all(format!(" {}", buffer.format(*weight)).as_bytes())?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_vertex_line(&mut self, vertex: &Vertex) -> Result<(), io::Error> {
        self.stream.write_all(
            format!(
                "V {} {} {} {} {} {} 0 {} {}",
                vertex.barcode,
                vertex.status,
                ryu::Buffer::new().format(vertex.x),
                ryu::Buffer::new().format(vertex.y),
                ryu::Buffer::new().format(vertex.z),
                ryu::Buffer::new().format(vertex.t),
                vertex.particles_in.len() + vertex.particles_out.len(),
                vertex.weights.len()
            )
            .as_bytes(),
        )?;
        for weight in &vertex.weights {
            self.stream.write_all(format!(" {weight}").as_bytes())?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_particle_line(&mut self, particle: &Particle) -> Result<(), io::Error> {
        self.stream.write_all(
            format!(
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
            )
            .as_bytes(),
        )?;
        for (idx, val) in &particle.flows {
            self.stream
                .write_all(format!(" {} {}", idx, val).as_bytes())?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_weight_names_line(&mut self, names: &[String]) -> Result<(), io::Error> {
        self.stream
            .write_all(format!("N {}", names.len()).as_bytes())?;
        for name in names {
            self.stream
                .write_all(format!(r#" "{}""#, name).as_bytes())?;
        }
        self.stream.write_all(b"\n")
    }

    fn write_unit_line(&mut self, event: &Event) -> Result<(), io::Error> {
        self.stream
            .write_all(format!("U {:?} {:?}\n", event.energy_unit, event.length_unit).as_bytes())
    }

    fn write_cross_section_line(&mut self, xs: &CrossSection) -> Result<(), io::Error> {
        self.stream.write_all(
            format!(
                "C {} {}\n",
                ryu::Buffer::new().format(xs.cross_section),
                ryu::Buffer::new().format(xs.cross_section_error)
            )
            .as_bytes(),
        )
    }

    fn write_pdf_info_line(&mut self, pdf: &PdfInfo) -> Result<(), io::Error> {
        self.stream.write_all(
            format!(
                "F {} {} {} {} {} {} {} {} {}\n",
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
            .as_bytes(),
        )
    }

    fn write_heavy_ion_info_line(&mut self, hi: &HeavyIonInfo) -> Result<(), io::Error> {
        self.stream.write_all(
            format!(
                "H {} {} {} {} {} {} {} {} {} {} {} {} {}\n",
                hi.ncoll_hard,
                hi.npart_proj,
                hi.npart_targ,
                hi.ncoll,
                hi.spectator_neutrons,
                hi.spectator_protons,
                hi.n_nwounded_collisions,
                hi.nwounded_n_collisions,
                hi.nwounded_nwounded_collisions,
                ryu::Buffer::new().format(hi.impact_parameter),
                ryu::Buffer::new().format(hi.event_plane_angle),
                ryu::Buffer::new().format(hi.eccentricity),
                ryu::Buffer::new().format(hi.sigma_inel_nn),
            )
            .as_bytes(),
        )
    }
}

impl<T: Write> Drop for Writer<T> {
    fn drop(&mut self) {
        if !self.finished {
            error!("Hepmc2 writer dropped before finished.");
            error!("Call finish() manually to fix this error.");
            if let Err(err) = self.ref_finish() {
                error!("Error writing footer: {}", err);
            }
        }
    }
}
