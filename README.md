# XRDSimulator: Interactive Single-Crystal X-ray Diffraction Demo

**Developed by Dr. Yugang Zhang**  
**Date**: 2025-08-10

This project provides an interactive simulation of single-crystal X-ray diffraction from the (111) plane of an FCC crystal. It visualizes the diffraction geometry, detector signal, and rocking curve as a function of crystal rotation angles, lattice constant, and X-ray wavelength.

## ✨ Features

- Interactive 3D visualization of beam-crystal-detector geometry  
- Simulates Bragg diffraction and spot detection  
- Plots rocking curve (θ-scan) in real-time  
- Adjustable parameters: χ (X-rot), θ (Y-rot), φ (Z-rot), lattice constant (a), and X-ray wavelength (λ)


## 🧩 Notes on Crystal Rotations

The rotation order is **X → Y → Z**, applied as matrix multiplications:

\[
R = R_X(\chi) \cdot R_Y(\theta) \cdot R_Z(\phi)
\]

These correspond to:

- **χ (chi)**: Rotation about the **beam axis (X-axis)**  
- **θ (theta)**: Rocking scan (**Y-axis rotation** – critical for Bragg condition)  
- **φ (phi)**: Azimuthal twist about the **Z-axis (initial normal direction)**

---

## 📊 Intensity Model

The intensity is modeled using a **Gaussian distribution** centered at the Bragg condition:

\[
I = \exp\left( -\left( \frac{\Delta \theta}{1^\circ} \right)^2 \right)
\]

Where:
- \( \Delta \theta \) is the deviation (in radians) from the ideal Bragg angle.
- The function creates a sharp, symmetric peak around the Bragg angle.
- Peak width is determined by the 1° standard deviation (adjustable for sharpness).

---

## 📍 Detector Projection

If the **Bragg condition** is satisfied **and** the scattered ray is **forward-directed** (i.e., positive X component), then the detector projection is calculated:

The outgoing beam vector is:

\[
\vec{k}_{\text{out}} = \vec{k}_{\text{in}} - 2(\vec{k}_{\text{in}} \cdot \vec{n}) \vec{n}
\]

Where:
- \( \vec{k}_{\text{in}} \) is the normalized incoming X-ray beam.
- \( \vec{n} \) is the rotated (111) plane normal vector.

The detector coordinates \( (y, z) \) are determined by intersecting this vector with a plane perpendicular to X at distance `L = detector_distance`.

---

## 📚 Integration Example

```python
from xrdsimulator import XRDSimulator

sim = XRDSimulator(a=4.08, wavelength=1.54)

# Rotate the crystal
normal = sim.rotate_crystal(chi_deg=10, theta_deg=19.5, phi_deg=0)

# Calculate diffraction
y, z, intensity = sim.calculate_diffraction(normal)

# Visualize geometry and rocking curve
fig = sim.visualize_setup(chi_deg=10, theta_deg=19.5, phi_deg=0)

🧠 Limitations & Assumptions
Only simulates (111) diffraction of an FCC crystal.

Assumes monochromatic X-ray beam (no polychromatic support).

Uses a flat detector in the YZ plane, fixed at X = +L.

Assumes perfect crystal with no:

Mosaic spread

Orientation distribution

Beam divergence

Intensity model is Gaussian, not based on full dynamical theory.

🛠 Future Improvements (Optional Ideas)
Add support for multiple crystallographic planes: (100), (110), etc.

Allow user-defined reciprocal lattice vectors.

Simulate multiple wavelengths or continuous spectra.

Integrate 3D reciprocal space mapping and orientation pole figures.

Implement Ewald sphere visualization and intersection logic.

Export interactive or static snapshots as image/video/GIF.

Add sample mosaicity or beam divergence effects.

📬 Contact
For suggestions, bug reports, or feature requests:

Dr. Yugang Zhang
Brookhaven National Lab
📧 yuzhang@bnl.gov

---

## 🧪 Installation

Create a virtual environment (optional but recommended), then install dependencies:

```bash
pip install -r requirements.txt

# 📘 XRDSimulator Class Documentation

This document provides a comprehensive technical overview of the `XRDSimulator` class used to simulate and visualize single-crystal X-ray diffraction based on Bragg’s law and crystal rotation.

---

## 📦 Class: `XRDSimulator`

Simulates diffraction from the (111) plane of an FCC crystal. Includes 3D geometry visualization, detector spot calculation, and rocking curve intensity analysis.

### 🔧 Constructor

```python
XRDSimulator(a=4.08, wavelength=1.54)


