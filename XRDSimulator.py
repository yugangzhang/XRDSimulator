# Dr. Yugang Zhang Developed at 2025/8/10 for the demo of single crystal diffraction pattern as a function of rotation angles of samples, lattice constant, and X-ray wavelength
# Please contact Dr. Zhang if you have questions.
# --------------------------------------------------------------
#  IMPORTS
# --------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D      # noqa: F401 (keeps the 3-D projection import)
import ipywidgets as widgets
from IPython.display import display, clear_output

# --------------------------------------------------------------
#  XRDSimulator  
# --------------------------------------------------------------
class XRDSimulator:
    """
    Simple X-ray-diffraction visualiser for a single (111) plane of an FCC crystal.
    Angles are supplied **in degrees** (the UI works with degrees for readability).
    """
    def __init__(self, a=4.08, wavelength=1.54):
        '''
        Parameters:
        a (float): Lattice constant in angstroms (default: 4.08 Å for FCC metals like Ni or Cu).
        
        wavelength (float): X-ray wavelength in angstroms (default: 1.54 Å for Cu Kα radiation).
        
        Attributes:
        self.a: Lattice constant.
        
        self.wavelength: X-ray wavelength.
        
        self.d111: Interplanar spacing for the (111) plane.
        
        self.theta_bragg: Bragg angle for the (111) reflection (in radians).
        
        self.initial_normal: Unit vector for the (111) plane in the initial (unrotated) crystal frame.
        
        self.k_in: Incoming X-ray beam vector (assumed along +X).
        
        self.detector_distance: Distance from sample to detector (in mm).
        
        self.detector_size: Width/height of square detector (in mm).

        '''
        # crystal
        self.a = a                              # lattice constant (Å)
        self.d111 = self.a / np.sqrt(3)         # d-spacing for (111)

        # X-ray
        self.wavelength = wavelength            # Å
        self.theta_bragg = np.arcsin(self.wavelength / (2 * self.d111))

        # geometry
        self.initial_normal = np.array([0, 0, 1])   # (111) normal in crystal frame
        self.k_in = np.array([1, 0, 0])             # incident beam (along +X)

        # detector
        self.detector_distance = 100.0   # mm (distance from crystal to detector plane)
        self.detector_size = 200.0       # mm (square detector)

    # ------------------------------------------------------------------
    #  Rotation matrices (angles **in radians**)
    # ------------------------------------------------------------------
    @staticmethod
    def _rot_x(angle):
        '''
        Returns a rotation matrix for rotation about the X-axis.

        Input: angle (float) – angle in radians

        Returns: 3x3 np.ndarray (rotation matrix)
        '''
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[1, 0, 0], [0, c, -s], [0, s,  c]])

    @staticmethod
    def _rot_y(angle):
        '''
        Returns a rotation matrix for rotation about the Y-axis.

        Input: angle (float) – angle in radians

        Returns: 3x3 np.ndarray
        '''
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[ c, 0, s], [ 0, 1, 0], [-s, 0, c]])

    @staticmethod
    def _rot_z(angle):
        '''
        Returns a rotation matrix for rotation about the Z-axis.

        Input: angle (float) – angle in radians

        Returns: 3x3 np.ndarray

        '''
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, -s, 0], [s,  c, 0], [0,  0, 1]])

    # ------------------------------------------------------------------
    #  Public API – rotation supplied in **degrees**
    # ------------------------------------------------------------------
    def rotate_crystal(self, chi_deg=0.0, theta_deg=0.0, phi_deg=0.0):
        """Return the (111) normal after X-χ, Y-θ, Z-φ rotations (degrees).
        Applies X, Y, Z rotations (in that order) to the initial (111) plane normal.

        Inputs:
        
        chi_deg (float) – Rotation about X-axis (χ), in degrees
        
        theta_deg (float) – Rotation about Y-axis (θ), in degrees
        
        phi_deg (float) – Rotation about Z-axis (φ), in degrees
        
        Returns:
        
        rotated_normal (np.ndarray, shape=(3,)) – Unit vector of the rotated (111) plane normal.

        """
        chi   = np.radians(chi_deg)   # about X
        theta = np.radians(theta_deg) # about Y
        phi   = np.radians(phi_deg)   # about Z

        R = (self._rot_x(chi) @ self._rot_y(theta) @ self._rot_z(phi))
        n = R @ self.initial_normal
        return n / np.linalg.norm(n)

    # ------------------------------------------------------------------
    #  Diffraction calculation -- ❗ FIX #1 IS HERE ❗
    # ------------------------------------------------------------------
    def calculate_diffraction(self, plane_normal):
        """Return detector (y,z) position (mm) and intensity (0-1).
        Calculates whether the provided plane normal satisfies the Bragg condition with the incoming beam. If satisfied, computes the detector coordinates and intensity.
        
        Input:
        
        plane_normal (np.ndarray) – A unit vector representing the normal of a crystal plane in lab coordinates.
        
        Returns:
        
        det_y (float or None) – Y-coordinate on detector (mm)
        
        det_z (float or None) – Z-coordinate on detector (mm)
        
        intensity (float) – Value between 0 and 1 depending on angular deviation from the Bragg condition. Peak is Gaussian-shaped.
        
        Bragg Condition:
        
        The angle between k_in and n must be close to      𝜋/2 - 𝜃_B         
        
        A tolerance of ±5° (in radians) is used.
        
        Forward Diffraction Check:
        
        Rejects backward reflections (k_out[0] <= 0.1)        
        
        """
        n = plane_normal / np.linalg.norm(plane_normal)
        k_in_norm = self.k_in / np.linalg.norm(self.k_in)

        # This is the angle between the beam and the PLANE NORMAL
        angle_to_normal = np.arccos(np.clip(np.dot(k_in_norm, n), -1, 1))

        # The Bragg condition requires the angle to the NORMAL to be 90° - θ_bragg
        # This is the crucial correction.
        bragg_condition_angle = np.pi/2 - self.theta_bragg

        if np.abs(angle_to_normal - bragg_condition_angle) < np.radians(5.0): # Use a tolerance
            k_out = k_in_norm - 2 * np.dot(k_in_norm, n) * n
            # Intensity is a Gaussian based on how close we are to the exact angle
            intensity = np.exp(-((angle_to_normal - bragg_condition_angle) / np.radians(1.0)) ** 2)

            if k_out[0] > 0.1:              # forward-going ray
                t = self.detector_distance / k_out[0]
                det_y = k_out[1] * t
                det_z = k_out[2] * t
                return det_y, det_z, intensity

        return None, None, 0.0

    # ------------------------------------------------------------------
    #  Visualiser – ❗ FIX #2 IS HERE ❗
    # ------------------------------------------------------------------
    def visualize_setup(self, chi_deg=0.0, theta_deg=0.0, phi_deg=0.0):
        """Draw the four panels for the supplied rotation angles (degrees).
        Creates a 2×2 figure panel for the specified rotation state:
        
        3D Geometry View: Lab axes and (111) normal.
        
        Detector View: Simulated diffraction spot (if present).
        
        Rocking Curve: Intensity vs θ (Y-rotation).
        
        Info Panel: Lattice constants, Bragg angle, and spot status.
        
        Inputs:
        
        chi_deg, theta_deg, phi_deg – Rotation angles in degrees.
        
        Returns:
        
        matplotlib.Figure object for display.

        """
        from matplotlib.gridspec import GridSpec
        #fig = plt.figure(figsize=(14, 9))
        fig = plt.figure(figsize=(14, 9), constrained_layout=True)

        fig.suptitle(f"χ={chi_deg:.1f}°, θ={theta_deg:.1f}°, φ={phi_deg:.1f}°", fontsize=14)
        gs = GridSpec(nrows=2, ncols=3, figure=fig)
        ax3d    = fig.add_subplot(gs[:, :2], projection='3d')  # big 3D plot spans 2 rows, 2 cols
        ax_det  = fig.add_subplot(gs[0, 2])                    # top right
        ax_bragg = fig.add_subplot(gs[1, 2])                   # bottom right    
        # ------------------------------------------------------------
        # 1. 3D View – Beam, Plane Normal, Detector Plane
        # ------------------------------------------------------------
        #ax3d = fig.add_subplot(221, projection='3d')
        ax3d.set_box_aspect([1, 1, 1])
    
        # Plot lab axes
        axes_length = 1.0
        ax3d.quiver(0, 0, 0, 1, 0, 0, color='k', linewidth=1, arrow_length_ratio=0.1)
        #ax3d.text(1.05, 0, 0, '[100]', color='k')
        ax3d.text(1.25, 0, 0, 'x', color='k')
        ax3d.quiver(0, 0, 0, 0, 1, 0, color='k', linewidth=1, arrow_length_ratio=0.1)
        #ax3d.text(0, 1.05, 0, '[010]', color='k')
        ax3d.text(0, 1.25, 0, 'y', color='k')
        ax3d.quiver(0, 0, 0, 0, 0, 1, color='k', linewidth=1, arrow_length_ratio=0.1)
        #ax3d.text(0, 0, 1.05, '[001]', color='k')
        ax3d.text(0, 0, 1.25, 'z', color='k')
    
        # Beam direction (incident)
        ax3d.quiver(0, 0, 0, 1, 0, 0, color='blue', linewidth=2, arrow_length_ratio=0.1)
        ax3d.text(0.9, 0.1, 0.1, 'k_in →', color='blue')
    
        # Rotated (111) normal
        n = self.rotate_crystal(chi_deg, theta_deg, phi_deg)
        ax3d.quiver(0, 0, 0, n[0], n[1], n[2], color='purple', linewidth=3, arrow_length_ratio=0.1)
        ax3d.text(n[0]*1.1, n[1]*1.1, n[2]*1.1, '(111)', color='purple')
    
        # Optional: Bragg-reflected ray
        det_y, det_z, intensity = self.calculate_diffraction(n)
        if intensity > 0.01:
            k_in_norm = self.k_in / np.linalg.norm(self.k_in)
            k_out = k_in_norm - 2 * np.dot(k_in_norm, n) * n
            ax3d.quiver(0, 0, 0, k_out[0], k_out[1], k_out[2], color='red', linewidth=2, arrow_length_ratio=0.1)
            ax3d.text(k_out[0], k_out[1], k_out[2], 'k_out →', color='red')
    
        # Detector plane
        d = self.detector_distance
        # Detector surface grid (at x = detector_distance)
        plane_size = self.detector_size / 100  # scale to match axes (~2x2 units)
        x_det = np.full((2, 2), self.detector_distance / 100)
        y_det = np.array([[-1, 1], [-1, 1]]) * (plane_size / 2)
        z_det = np.array([[-1, -1], [1, 1]]) * (plane_size / 2)        
        ax3d.plot_surface(x_det, y_det, z_det, color='gray', alpha=0.2, edgecolor='none')    
        ax3d.set_xlim(-1, 1)
        ax3d.set_ylim(-1, 1)
        ax3d.set_zlim(-1, 1)
        ax3d.set_xlabel('X')
        ax3d.set_ylabel('Y')
        ax3d.set_zlabel('Z')
        ax3d.set_title('3D Diffraction Geometry')
    
        # ------------------------------------------------------------
        # 2. Detector View
        # ------------------------------------------------------------
        #ax_det = fig.add_subplot(222)
        if intensity > 0.01:
            circ = Circle((det_y, det_z), radius=2 + 15 * intensity, color='red', alpha=0.7)
            ax_det.add_patch(circ)
        ax_det.set_xlim(-self.detector_size/2, self.detector_size/2)
        ax_det.set_ylim(-self.detector_size/2, self.detector_size/2)
        ax_det.set_xlabel('Y (mm)')
        ax_det.set_ylabel('Z (mm)')
        ax_det.set_title('Detector View')
        ax_det.grid(True, ls='--', alpha=0.5)
        ax_det.set_aspect('equal')
    
        # ------------------------------------------------------------
        # 3. Rocking Curve (θ-scan)
        # ------------------------------------------------------------
        #ax_bragg = fig.add_subplot(223)
        theta_scan_range = np.linspace(theta_deg - 15, theta_deg + 15, 100)
        intensity_curve = []
        for th_val in theta_scan_range:
            temp_normal = self.rotate_crystal(chi_deg, th_val, phi_deg)
            _, _, temp_intensity = self.calculate_diffraction(temp_normal)
            intensity_curve.append(temp_intensity)
        
        ax_bragg.plot(theta_scan_range, intensity_curve, color='steelblue',
                      label=f'Rocking curve at χ={chi_deg:.1f}°, φ={phi_deg:.1f}°')
        ax_bragg.axvline(theta_deg, color='purple', lw=2, label=f'Current θ={theta_deg:.1f}°')
        ax_bragg.set_xlabel('θ (Y-rotation) [degrees]')
        ax_bragg.set_ylabel('Calculated Intensity')
        ax_bragg.set_title('Calculated θ-scan (Rocking Curve)')
        ax_bragg.legend()
        ax_bragg.grid(True, alpha=0.4)
    
        # ------------------------------------------------------------
        # 4. Info Panel
        # ------------------------------------------------------------
        # ax_info = fig.add_subplot(224)
        # ax_info.axis('off')
        theta_bragg_deg = np.degrees(self.theta_bragg)
        txt = (
            f"CRYSTAL & BEAM\n"
            f"  a         = {self.a:.3f} Å\n"
            f"  λ         = {self.wavelength:.3f} Å\n"
            f"  d₁₁₁      = {self.d111:.3f} Å\n"
            f"  Bragg θ_B = {theta_bragg_deg:.2f}° (glancing)\n\n"
            #f"STATUS\n"
        )
        if intensity > 0.01:
            txt += (
                f"\N{WHITE HEAVY CHECK MARK}  Bragg Condition Met\n"
                f"Spot (Y,Z)  = ({det_y:.1f}, {det_z:.1f}) mm\n"
                f"Intensity   = {intensity:.3f}"
            )
            box_color = '#d6f5d6'  # light green background
            text_color = 'green'
        else:
            txt += "✗ Not in Bragg Condition"
            box_color = '#f5d6d6'  # light red background
            text_color = 'red'

        #ax_info.text(0.05, 0.95, txt, fontsize=11, family='monospace', verticalalignment='top')
        #ax_info.set_title('Info Panel')
        # Display in top right of ax_bragg
        ax_bragg.text(1.05, 0.95, txt,
              transform=ax_bragg.transAxes,
              fontsize=10,
              va='top', ha='left',
              family='monospace',
              color=text_color,
              bbox=dict(facecolor=box_color, alpha=0.6, edgecolor='gray', boxstyle='round,pad=0.4') )
        #fig.tight_layout(rect=[0, 0, 1, 0.96])
        return fig

# --------------------------------------------------------------
#  The rest of your UI code is perfect and needs no changes.
# --------------------------------------------------------------
sim = XRDSimulator()

sim = XRDSimulator()

def interactive_simulation():
    '''
    Interactive simulation for single-crystal X-ray diffraction from the (111) plane
    of an FCC crystal, using a Bragg-reflection geometry.

    -----------------------------------------------------------
    GEOMETRY SETUP:
    -----------------------------------------------------------

    • Incident X-ray beam travels along the +X direction:
        --> k_in = [1, 0, 0]  ← [100] in crystal/lab space

    • Initial (111) crystal plane has its normal aligned along +Z:
        --> n_111 = [0, 0, 1]  ← [111] is parallel to lab-frame [001]

    • Detector is a flat square positioned perpendicular to the beam:
        --> Located at X = L, spanning the YZ plane (200 mm × 200 mm)

    -----------------------------------------------------------
    ROTATION CONTROLS:
    -----------------------------------------------------------

    ▶️ χ (X-rotation): Rotation around beam direction (X)
        • Spot moves in a circle (Debye-Scherrer cone) if Bragg is already satisfied.
        • Radius = L · tan(2θ_B)
        • Intensity remains "on" throughout χ rotation.

    ▶️ θ (Y-rotation): Rocking scan
        • Controls whether Bragg condition is met.
        • Spot appears sharply only when θ ≈ Bragg angle (θ_B).
        • Results in narrow peak in rocking curve.

    ▶️ φ (Z-rotation): Rotation around (111) normal
        • Has no effect — the plane normal is unchanged.
        • No diffraction unless other angles are correct.

    -----------------------------------------------------------
    OBSERVATIONS:
    -----------------------------------------------------------

    ✔️ Y-Rotation (θ):
        • Most sensitive angle. A small deviation from θ_B suppresses intensity.
        • Produces sharp Gaussian-shaped rocking curve.

    ✔️ X-Rotation (χ):
        • Causes diffraction spot to move in a circle.
        • Useful for visualizing cone of diffracted rays.

    ✔️ Z-Rotation (φ):
        • No effect — the normal vector is unchanged.
        • Bragg condition cannot be satisfied alone.

    -----------------------------------------------------------
    UI DESCRIPTION:
    -----------------------------------------------------------

    This widget interface includes sliders and input boxes for:
        - χ, θ, φ: Rotation angles in degrees
        - a: Lattice constant (Å)
        - λ: X-ray wavelength (Å)

    Live updates include:
        - 3D geometry of beam and rotated plane
        - Detector view with diffraction spot (if any)
        - Rocking curve (θ-scan)
        - Info panel with simulation parameters and status
    '''

    def linked_pair(label, min_, max_, step_, init_):
        sld = widgets.FloatSlider(value=init_, min=min_, max=max_, step=step_, description=label, layout=widgets.Layout(width='70%'), continuous_update=False)
        txt = widgets.BoundedFloatText(value=init_, min=min_, max=max_, step=step_, layout=widgets.Layout(width='25%'))
        widgets.jslink((sld, 'value'), (txt, 'value'))
        return sld, txt

    # Rotation sliders
    chi_sld, chi_txt = linked_pair('χ (X-rot)', -180, 180, 0.1, 0.0)
    th_sld,  th_txt  = linked_pair('θ (Y-rot)', -90, 90, 0.1, 0.0)
    phi_sld, phi_txt = linked_pair('φ (Z-rot)', -90, 90, 0.1, 0.0)

    # Lattice & wavelength sliders
    a_sld = widgets.FloatSlider(value=sim.a, min=2.0, max=6.0, step=0.01, description='a (Å)', continuous_update=False)
    lam_sld = widgets.FloatSlider(value=sim.wavelength, min=0.5, max=3.0, step=0.01, description='λ (Å)', continuous_update=False)

    # Layout
    out = widgets.Output()
    row_chi = widgets.HBox([chi_sld, chi_txt])
    row_th  = widgets.HBox([th_sld,  th_txt])
    row_phi = widgets.HBox([phi_sld, phi_txt])
    row_mat = widgets.HBox([a_sld, lam_sld])
    ui = widgets.VBox([out, row_chi, row_th, row_phi, row_mat])

    # Plot updater
    def _draw(_=None):
        sim.a = a_sld.value
        sim.wavelength = lam_sld.value
        sim.d111 = sim.a / np.sqrt(3)
        sim.theta_bragg = np.arcsin(sim.wavelength / (2 * sim.d111))
        with out:
            clear_output(wait=True)
            fig = sim.visualize_setup(
                chi_deg=chi_sld.value,
                theta_deg=th_sld.value,
                phi_deg=phi_sld.value
            )
            plt.show(fig)

    # Attach callbacks
    for w in (chi_sld, th_sld, phi_sld, chi_txt, th_txt, phi_txt, a_sld, lam_sld):
        w.observe(_draw, names='value')

    display(ui)
    _draw()


# --- RUN THE INTERACTIVE UI ---
# (Make sure to run this in a Jupyter Notebook or JupyterLab environment)
#interactive_simulation()