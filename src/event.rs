use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use strum::EnumString;

/// Scattering event
#[derive(Debug, PartialEq, Default, Clone, Serialize, Deserialize)]
pub struct Event {
    pub number: i32,
    pub mpi: i32,
    pub scale: f64,
    pub alpha_qcd: f64,
    pub alpha_qed: f64,
    pub signal_process_id: i32,
    pub signal_process_vertex: i32,
    pub random_states: Vec<i32>,
    pub weights: Vec<f64>,
    pub weight_names: Vec<String>,
    pub vertices: Vec<Vertex>,
    pub xs: CrossSection,
    pub pdf_info: PdfInfo,
    pub energy_unit: EnergyUnit,
    pub length_unit: LengthUnit,
    pub heavy_ion_info: Option<HeavyIonInfo>,
}

/// Interaction vertex
#[derive(Debug, PartialEq, Default, Clone, Serialize, Deserialize)]
pub struct Vertex {
    pub barcode: i32,
    pub status: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub t: f64,
    pub weights: Vec<f64>,
    pub particles_in: Vec<Particle>,
    pub particles_out: Vec<Particle>,
}

/// Particle
#[derive(Debug, PartialEq, Default, Clone, Serialize, Deserialize)]
pub struct Particle {
    pub id: i32,
    pub p: FourVector,
    pub m: f64,
    pub status: i32,
    pub theta: f64,
    pub phi: f64,
    pub flows: BTreeMap<i32, i32>,
    pub end_vtx: i32,
}

/// Simple Lorentz vector with components (t, x, y, z)
#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct FourVector(pub [f64; 4]);

impl FourVector {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn txyz(t: f64, x: f64, y: f64, z: f64) -> Self {
        FourVector([t, x, y, z])
    }
}

impl std::ops::Index<usize> for FourVector {
    type Output = f64;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl std::ops::IndexMut<usize> for FourVector {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.0[idx]
    }
}

/// Cross section with error
#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct CrossSection {
    pub cross_section: f64,
    pub cross_section_error: f64,
}

impl std::fmt::Display for CrossSection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ± {}", self.cross_section, self.cross_section_error)
    }
}

/// PDF information
#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct PdfInfo {
    pub parton_id: [i32; 2],
    pub x: [f64; 2],
    pub scale: f64,
    pub xf: [f64; 2],
    pub pdf_id: [i32; 2],
}

/// Information for heavy ion collisions
#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct HeavyIonInfo {
    pub ncoll_hard: i32,
    pub npart_proj: i32,
    pub npart_targ: i32,
    pub ncoll: i32,
    pub spectator_neutrons: i32,
    pub spectator_protons: i32,
    pub n_nwounded_collisions: i32,
    pub nwounded_n_collisions: i32,
    pub nwounded_nwounded_collisions: i32,
    pub impact_parameter: f64,
    pub event_plane_angle: f64,
    pub eccentricity: f64,
    pub sigma_inel_nn: f64,
}

/// Energy units
#[derive(
    EnumString, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone, Serialize, Deserialize,
)]
pub enum EnergyUnit {
    MEV,
    GEV,
}

impl std::default::Default for EnergyUnit {
    fn default() -> Self {
        Self::GEV
    }
}

/// Length units
#[derive(
    EnumString, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone, Serialize, Deserialize,
)]
pub enum LengthUnit {
    MM,
    CM,
}

impl std::default::Default for LengthUnit {
    fn default() -> Self {
        Self::CM
    }
}
