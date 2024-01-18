use std::collections::BTreeMap;
use std::default::Default;
use std::fmt::{self, Display};
use std::io;
use std::num::{ParseFloatError, TryFromIntError};

use crate::event::*;

use hepmc2_macros::read_bound;
use nom::{
    bytes::complete::{take_until, take_while1},
    character::complete::{char, i32, space1, u64},
    combinator::opt,
    number::complete::double,
    sequence::{delimited, preceded, tuple},
    IResult,
};
use thiserror::Error;

const BUF_SIZE: usize = 256;

/// Reader for the HepMC2 format
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
pub struct Reader<T> {
    stream: T,
    line: String,
    line_nr: usize,
}

impl<T> Reader<T> {
    /// Retrieve the underlying reader
    pub fn into_inner(self) -> T {
        self.stream
    }
}

#[read_bound]
impl<T> Reader<T> {
    /// Construct a new Reader
    pub fn new(stream: T) -> Self {
        stream.into()
    }
}

#[read_bound]
impl<T> From<T> for Reader<T> {
    fn from(stream: T) -> Self {
        Self {
            stream,
            line: String::with_capacity(BUF_SIZE),
            line_nr: 0,
        }
    }
}

#[read_bound]
impl<T> Reader<T> {
    #[maybe_async::maybe_async]
    async fn skip_headers(&mut self) -> Result<(), io::Error> {
        while self.line.trim().is_empty() || self.line.starts_with("HepMC") {
            self.line.clear();
            if self.stream.read_line(&mut self.line).await? == 0 {
                break;
            }
            self.line_nr += 1;
        }
        Ok(())
    }

    #[maybe_async::maybe_async]
    async fn parse_event_inner(&mut self) -> Result<Event, ParseError> {
        let mut event = parse_event_line(&self.line)?;
        loop {
            self.line.clear();
            if self.stream.read_line(&mut self.line).await? == 0 {
                break;
            };
            self.line_nr += 1;
            match self.line.as_bytes().first() {
                Some(b'E') => break,
                Some(b'V') => parse_vertex_line(&self.line, &mut event)?,
                Some(b'P') => parse_particle_line(&self.line, &mut event)?,
                Some(b'U') => parse_units_line(&self.line, &mut event)?,
                Some(b'F') => parse_pdf_info_line(&self.line, &mut event)?,
                Some(b'H') => {
                    if self.line.starts_with("HepMC") {
                        continue;
                    }
                    parse_heavy_ion_line(&self.line, &mut event)?
                }
                Some(b'N') => parse_weight_names_line(&self.line, &mut event)?,
                Some(b'C') => parse_xs_info_line(&self.line, &mut event)?,
                _ => {
                    if self.line.trim().is_empty() {
                        continue;
                    } else {
                        return Err(ParseError::BadPrefix);
                    }
                }
            };
        }
        Ok(event)
    }

    #[maybe_async::maybe_async]
    async fn parse_event(&mut self) -> Result<Event, LineParseError> {
        self.parse_event_inner()
            .await
            .map_err(|err| LineParseError {
                err,
                line: self.line.clone(),
                line_nr: self.line_nr,
            })
    }

    #[maybe_async::async_impl]
    /// Read the next event from the stream
    // NOTE: This function is a tempory measure! The proper way to replicate an iterator in async
    // context is to implement `futures::stream::Stream`. This should be done in a later
    // contribution.
    pub async fn next(&mut self) -> Option<std::result::Result<Event, LineParseError>> {
        if let Err(err) = self.skip_headers().await {
            return Some(Err(LineParseError {
                err: err.into(),
                line: self.line.clone(),
                line_nr: self.line_nr,
            }));
        }
        if self.line.is_empty() {
            return None;
        }
        Some(self.parse_event().await)
    }
}

fn whitespace(line: &str) -> IResult<&str, &str> {
    space1(line)
}

fn non_whitespace(line: &str) -> IResult<&str, &str> {
    take_while1(|c: char| !c.is_ascii_whitespace())(line)
}

fn ws_nonws(line: &str) -> IResult<&str, &str> {
    preceded(whitespace, non_whitespace)(line)
}

fn ws_i32(line: &str) -> IResult<&str, i32> {
    preceded(whitespace, i32)(line)
}

fn ws_u64(line: &str) -> IResult<&str, u64> {
    preceded(whitespace, u64)(line)
}

fn ws_double(line: &str) -> IResult<&str, f64> {
    preceded(whitespace, double)(line)
}

fn string(line: &str) -> IResult<&str, &str> {
    delimited(char('"'), take_until("\""), char('"'))(line)
}

fn parse_event_line(line: &str) -> Result<Event, ParseError> {
    let rest = &line[1..];

    let (rest, event_number) = ws_i32(rest)?;
    let (rest, mpi) = ws_i32(rest)?;
    let (rest, event_scale) = ws_double(rest)?;
    let (rest, alpha_qcd) = ws_double(rest)?;
    let (rest, alpha_qed) = ws_double(rest)?;
    let (rest, signal_process_id) = ws_i32(rest)?;
    let (rest, signal_process_vertex) = ws_i32(rest)?;
    let (rest, num_vertices) = ws_u64(rest)?;
    let num_vertices = num_vertices.try_into()?;
    let (rest, _beam1) = ws_nonws(rest)?;
    let (rest, _beam2) = ws_nonws(rest)?;
    let (mut rest, nrandom_states) = ws_u64(rest)?;

    let nrandom_states = nrandom_states.try_into()?;
    let mut random_states = Vec::with_capacity(nrandom_states);
    for _ in 0..nrandom_states {
        let (rem, random_state) = ws_i32(rest)?;
        rest = rem;
        random_states.push(random_state);
    }
    let (mut rest, nweights) = ws_u64(rest)?;
    let nweights = nweights.try_into()?;
    let mut weights = Vec::with_capacity(nweights);
    for _ in 0..nweights {
        let (rem, weight) = ws_double(rest)?;
        rest = rem;
        weights.push(weight);
    }
    let event = Event {
        number: event_number,
        mpi,
        scale: event_scale,
        alpha_qcd,
        alpha_qed,
        signal_process_id,
        signal_process_vertex,
        random_states,
        weights,
        vertices: Vec::with_capacity(num_vertices),
        weight_names: Default::default(),
        xs: Default::default(),
        energy_unit: Default::default(),
        length_unit: Default::default(),
        pdf_info: Default::default(),
        heavy_ion_info: None,
    };
    Ok(event)
}

fn parse_vertex_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let rest = &line[1..];
    let (rest, barcode) = ws_i32(rest)?;
    let (rest, status) = ws_i32(rest)?;
    let (rest, x) = ws_double(rest)?;
    let (rest, y) = ws_double(rest)?;
    let (rest, z) = ws_double(rest)?;
    let (rest, t) = ws_double(rest)?;
    let (rest, _num_orphans_int) = ws_i32(rest)?;
    let (rest, num_particles_out) = ws_u64(rest)?;
    let num_particles_out = num_particles_out.try_into()?;
    let (mut rest, num_weights) = ws_u64(rest)?;
    let num_weights = num_weights.try_into()?;
    let mut weights = Vec::with_capacity(num_weights);
    for _ in 0..num_weights {
        let (rem, weight) = ws_double(rest)?;
        rest = rem;
        weights.push(weight);
    }
    let vertex = Vertex {
        barcode,
        status,
        x,
        y,
        z,
        t,
        weights,
        particles_in: Vec::new(),
        particles_out: Vec::with_capacity(num_particles_out),
    };
    event.vertices.push(vertex);
    Ok(())
}

fn parse_particle_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let rest = &line[1..];
    let (rest, _barcode) = ws_i32(rest)?;
    let (rest, id) = ws_i32(rest)?;
    let (rest, px) = ws_double(rest)?;
    let (rest, py) = ws_double(rest)?;
    let (rest, pz) = ws_double(rest)?;
    let (rest, e) = ws_double(rest)?;
    let (rest, m) = ws_double(rest)?;
    let (rest, status) = ws_i32(rest)?;
    let (rest, theta) = ws_double(rest)?;
    let (rest, phi) = ws_double(rest)?;
    let (rest, end_vtx_code) = ws_i32(rest)?;
    let (mut rest, flowsize) = ws_i32(rest)?;
    let mut flows = BTreeMap::new();
    for _ in 0..flowsize {
        let (rem, flowidx) = ws_i32(rest)?;
        let (rem, flowval) = ws_i32(rem)?;
        rest = rem;
        flows.insert(flowidx, flowval);
    }
    let particle = Particle {
        id,
        p: FourVector::txyz(e, px, py, pz),
        m,
        status,
        theta,
        phi,
        flows,
        end_vtx: end_vtx_code,
    };
    // TODO: handling of end_vtx is ReaderAsciiHepMC2.cc is obscure and undocumented
    if let Some(vertex) = event.vertices.last_mut() {
        if particle.end_vtx == vertex.barcode {
            vertex.particles_in.push(particle);
        } else {
            vertex.particles_out.push(particle);
        }
    } else {
        return Err(ParseError::NoVertex);
    }
    Ok(())
}

fn parse_units_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let rest = &line[1..];

    let (rest, energy) = ws_nonws(rest)?;
    let (_rest, length) = ws_nonws(rest)?;
    event.energy_unit = energy.parse()?;
    event.length_unit = length.parse()?;
    Ok(())
}

fn parse_pdf_info_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let rest = &line[1..];

    let (rest, id0) = ws_i32(rest)?;
    let (rest, id1) = ws_i32(rest)?;
    let (rest, x0) = ws_double(rest)?;
    let (rest, x1) = ws_double(rest)?;
    let (rest, scale) = ws_double(rest)?;
    let (rest, xf0) = ws_double(rest)?;
    let (rest, xf1) = ws_double(rest)?;
    let (_rest, parsed) = tuple((
        whitespace,
        opt(i32), // pdf_id0
        whitespace,
        opt(i32), // pdf_id1
    ))(rest)?;
    let (_, pdf_id0, _, pdf_id1) = parsed;
    let pdf_info = PdfInfo {
        parton_id: [id0, id1],
        x: [x0, x1],
        scale,
        xf: [xf0, xf1],
        pdf_id: [pdf_id0.unwrap_or(0), pdf_id1.unwrap_or(0)],
    };
    event.pdf_info = pdf_info;
    Ok(())
}

fn parse_heavy_ion_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let rest = &line[1..];

    let (rest, ncoll_hard) = ws_i32(rest)?;
    let (rest, npart_proj) = ws_i32(rest)?;
    let (rest, npart_targ) = ws_i32(rest)?;
    let (rest, ncoll) = ws_i32(rest)?;
    let (rest, spectator_neutrons) = ws_i32(rest)?;
    let (rest, spectator_protons) = ws_i32(rest)?;
    let (rest, n_nwounded_collisions) = ws_i32(rest)?;
    let (rest, nwounded_n_collisions) = ws_i32(rest)?;
    let (rest, nwounded_nwounded_collisions) = ws_i32(rest)?;
    let (rest, impact_parameter) = ws_double(rest)?;
    let (rest, event_plane_angle) = ws_double(rest)?;
    let (rest, eccentricity) = ws_double(rest)?;
    let (_rest, sigma_inel_nn) = ws_double(rest)?;
    event.heavy_ion_info = Some(HeavyIonInfo {
        ncoll_hard,
        npart_proj,
        npart_targ,
        ncoll,
        spectator_neutrons,
        spectator_protons,
        n_nwounded_collisions,
        nwounded_n_collisions,
        nwounded_nwounded_collisions,
        impact_parameter,
        event_plane_angle,
        eccentricity,
        sigma_inel_nn,
    });
    Ok(())
}

fn parse_weight_names_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let rest = &line[1..];
    let (mut rest, nnames) = ws_u64(rest)?;
    let nnames = nnames.try_into()?;
    let mut weight_names = Vec::with_capacity(nnames);
    for _ in 0..nnames {
        let (rem, (_, name)) = tuple((whitespace, string))(rest)?;
        weight_names.push(name.to_owned());
        rest = rem;
    }
    event.weight_names = weight_names;
    Ok(())
}

fn parse_xs_info_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let rest = &line[1..];

    let (rest, cross_section) = ws_double(rest)?;
    let (_rest, cross_section_error) = ws_double(rest)?;
    event.xs = CrossSection {
        cross_section,
        cross_section_error,
    };
    Ok(())
}

#[maybe_async::sync_impl]
impl<T: std::io::BufRead> Iterator for Reader<T> {
    type Item = Result<Event, LineParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Err(err) = self.skip_headers() {
            return Some(Err(LineParseError {
                err: err.into(),
                line: self.line.clone(),
                line_nr: self.line_nr,
            }));
        }
        if self.line.is_empty() {
            return None;
        }
        Some(self.parse_event())
    }
}

/// Error when parsing a line
#[derive(Debug)]
pub struct LineParseError {
    /// The actual error
    pub err: ParseError,
    /// The line where the error occurred
    pub line: String,
    /// The line number where the error occurred
    pub line_nr: usize,
}

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("I/O error")]
    Io(#[from] io::Error),
    #[error("Parsing error: {0}")]
    Parse(String),
    #[error("Integer conversion error")]
    ConvertInt(#[from] TryFromIntError),
    #[error("Float conversion error")]
    ConvertFloat(#[from] ParseFloatError),
    #[error("Enum parsing error")]
    StrumErr(#[from] strum::ParseError),
    #[error("Unrecognized prefix")]
    BadPrefix,
    #[error("Tried to add particle without vertex")]
    NoVertex,
}

impl<T: Display> From<nom::Err<T>> for ParseError {
    fn from(err: nom::Err<T>) -> Self {
        match err {
            nom::Err::Failure(err) => ParseError::Parse(err.to_string()),
            nom::Err::Error(err) => ParseError::Parse(err.to_string()),
            _ => unreachable!(),
        }
    }
}

impl Display for LineParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\n in line {}:\n{}", self.err, self.line_nr, self.line)
    }
}

impl std::error::Error for LineParseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&(self.err))
    }
}
