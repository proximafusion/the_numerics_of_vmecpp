\chapter{VMEC Outputs}

Legacy Fortran VMEC produces four main output files:
\texttt{wout},
\texttt{jxbout},
\texttt{mercier} and
\texttt{threed1}.
%
These are merged into a single HDF5 file in VMEC++.
HDF5 is chosen, because it allows (in principle)
to adjust the other tools in the VMEC toolchain
to adopt the new VMEC++ output file format,
since a well-established HDF5 Fortran interface exists
and is maintained by the HDF group.
%
An intermediate conversion program from the single HDF5
into the legacy output files is available as part of the Python code vendored with VMEC++,
since most of the VMEC toolchain relies on the legacy output files.

The HDF5 file produced by VMEC++ has the following structure:
\begin{itemize}
 \item \texttt{/indata} - This is a 1:1 duplicate of the input data used to create this VMEC++ run.
 \item \texttt{/jxbout} - This is the equivalent of the \texttt{jxbout} output file produced by Fortran VMEC.
 \item \texttt{/mercier} - This is the equivalent of the \texttt{mercier} output file produced by Fortran VMEC.
 \item \texttt{/threed1} - This is the equivalent of the \texttt{threed1} output file produced by Fortran VMEC.
 \item \texttt{/wout} - This is the equivalent of the \texttt{wout} output file produced by Fortran VMEC.
\end{itemize}

\section{\texttt{indata}}
The \texttt{indata} members are specified and described above.
Their definitions are therefore not repeated here.

% input_extension
% mgrid_file
% pcurr_type
% pmass_type
% piota_type
% gamma
% nfp
% ns
% mpol
% ntor
% lasym__logical__
% lfreeb__logical__
% nextcur
% extcur
% curlabel

% am
% ac
% ai
% am_aux_s
% am_aux_f
% ac_aux_s
% ac_aux_f
% ai_aux_s
% ai_aux_f


% niter
% ftolv

% \section{\texttt{jxbout}}

% legend -> better move into docs and/or doc attributes
% mpol (already in INDATA)
% ntor (already in INDATA)
% radial_surfaces -> duplicate of ns
% poloidal_grid_points -> effective ntheta
% toroidal_grid_points -> effective nzeta

% radial profiles:
% phin
% avforce -> duplicate of equif?
% surf_av_jdotb
% bdotgradv
% pprime
% aminfor
% amaxfor

% full-volume:
% sqrt_g__bdotk
% sqrt_g_
% jsupu
% jsupv
% jsups
% bsupu
% bsupv
% jcrossb
% jxb_gradp
% bsubu
% bsubv
% bsubs

% \section{\texttt{mercier}}

% first table:
% s
% phi
% iota
% shear
% vp
% well
% itor
% itor'
% pres
% pres'

% result table:
% s
% DMerc
% DShear
% DCurr
% DWell
% DGeod

% \section{\texttt{threed1}}

% first table:
% s
% radial_force
% toroidal_flux
% iota
% avg_jsupu
% avg_jsupv
% d_volume_d_phi
% d_pressure_d_phi
% spectral_width
% pressure
% buco_full
% bvco_full
% j_dot_b
% b_dot_b

% geometric and magnetic quantities:
% toroidal_flux
% circum_p
% surf_area_p
% cross_area_p
% volume_p
% Rmajor_p
% Aminor_p
% aspect
% kappa_p
% rcen
% aminr1
% pavg
% factor
% b0
% rmax_surf
% rmin_surf
% zmax_surf
% bmin
% bmax
% waist
% height
% betapol
% betatot
% betator
% VolAvgB
% IonLarmor
% jpar_perp
% jparPS_perp
% toroidal_current
% rbtor
% rbtor0
% psi
% ygeo
% yinden
% yellip
% ytrian
% yshift
% loc_jpar_perp
% loc_jparPS_perp

% volumetrics:
% int_p
% avg_p
% int_bpol
% avg_bpol
% int_btor
% avg_btor
% int_modb
% avg_modb
% int_ekin
% avg_ekin

% axis geometry:
% raxis_symm
% zaxis_symm
% raxis_asym
% zaxis_asym

% betas:
% betatot
% betapol
% betator
% rbtor
% betaxis
% betstr

% Shafranov Integrals:
% scaling_ratio
% r_lao
% f_lao
% f_geo
% smaleli
% betai
% musubi
% lambda
% s11
% s12
% s13
% s2
% s3
% delta1
% delta2
% delta3

\section{\texttt{wout}}

The first set of members in the \texttt{wout} group
contains the following members:
\begin{itemize}
 \item \texttt{version} \\
       Type: \texttt{string} \\
       This is the version number of VMEC++, e.g.,
       \texttt{"8.52"} for the initial release based on serial VMEC 8.52 from STELLOPT's \texttt{v251} branch,
       on which \texttt{educational\_VMEC} is based, on which in turn VMEC++ is based.
 \item \texttt{sign\_of\_jacobian} \\
       Type: \texttt{int} \\
       This is the sign of the Jacobian of the coordintate transform between flux coordinates
       and cylindrical coordinates. It encodes the handedness of the coordinate system.
       As in Fortran VMEC (where this is called \texttt{signgs}), this is fixed to $-1$ as of now.
 \item \texttt{mgrid\_mode} \\
       Type: \texttt{string} \\
       In case of a free-boundary run,
       this indicates if the mgrid file was normalized to unit currents~(\texttt{"S"}) or not~(\texttt{"R"}).
       For a given mgrid file, this changes the interpretation of the coil currents.
 \item \texttt{ier\_flag} \\
       Type: \texttt{int} \\
       This is the status code indicating success of or problems during this VMEC++ run.
\end{itemize}
%
The next set of members in the \texttt{wout} group
indicates the convergence of this VMEC++ run:
\begin{itemize}
 \item \texttt{fsqt} \\
       Type: \texttt{double[]} \\
       Evolution of the total force residual along the run.
       Not implemented yet.
 \item \texttt{wdot} \\
       Type: \texttt{double[]} \\
       Evolution of the MHD energy decay along the run.
       Not implemented yet.
 \item \texttt{nfev} \\
       Type: \texttt{int} \\
       Number of evaluations of the ideal MHD forward model during this run.
 \item \texttt{fsqr} \\
       Type: \texttt{double} \\
       Invariant force residual of the force on $R$ at end of the run.
 \item \texttt{fsqz} \\
       Type: \texttt{double} \\
       Invariant force residual of the force on $Z$ at end of the run.
 \item \texttt{fsql} \\
       Type: \texttt{double} \\
       Invariant force residual of the force on $\lambda$ at end of the run.
\end{itemize}
%
A number of scalar plasma properties are available
in the next set of members in the \texttt{wout} group:
\begin{itemize}
 \item \texttt{wb} \\
       Type: \texttt{double} \\
       Magnetic energy: volume integral of $\vert\mathbf{B}\vert^2/(2 \mu_0)$.
 \item \texttt{wp} \\
       Type: \texttt{double} \\
       Kinetic energy: volume integral of $p$.
 \item \texttt{Rmajor\_p} \\
       Type: \texttt{double} \\
       Major radius of the plasma.
 \item \texttt{Aminor\_p} \\
       Type: \texttt{double} \\
       Minor radius of the plasma.
 \item \texttt{volume\_p} \\
       Type: \texttt{double} \\
       Plasma volume.
 \item \texttt{b0} \\
       Type: \texttt{double} \\
       Toroidal magnetic flux density from poloidal current and magnetic axis position at $\varphi=0$:
       \begin{equation}
        \texttt{b0} = F_\psi(0) / R(s=0, \theta=0, \varphi=0)
       \end{equation}
       with the poloidal current at the magnetic axis~$F_\psi(0)$
       (obtained by linear extrapolation from the innermost two half-grid points):
       \begin{equation}
        F_\psi(0) = \frac{3}{2} \langle B_\varphi \rangle(s_{0.5}) - \frac{1}{2} \langle B_\varphi \rangle(s_{1.5})
       \end{equation}
       and in there the covariant magnetic field component~$B_\varphi$ averaged over flux surfaces:
       \begin{equation}
        \langle B_\varphi \rangle (s)= \int\limits_0^{2 \pi} \int\limits_0^{2 \pi} B_\varphi(s ,\theta, \varphi) \,\mathrm{d}\theta \,\mathrm{d}\varphi
       \end{equation}
       (\texttt{bvcoH} in VMEC++, on half-grid).
       One can think of this as approximating the toroidal field as originating from a line current along the $z$ axis,
       which combines all poloidal currents in the coils and the plasma,
       and then asking what magnetic field at the radius of the magnetic axis at $\varphi=0$ this would create.
       Because the quantity is divided by $R$ at $\varphi = 0$,
       it is not a toroidal average and should not be used when, e.g., scaling equilibria based on magnetic field strength
       \footnote{Thanks to Damien Huet for stressing this important point.}.
 \item \texttt{rmax\_surf} \\
       Type: \texttt{double} \\
       Maximum of $R$ on the plasma boundary over all grid points.
 \item \texttt{rmin\_surf} \\
       Type: \texttt{double} \\
       Minimum of $R$ on the plasma boundary over all grid points.
 \item \texttt{zmax\_surf} \\
       Type: \texttt{double} \\
       Maximum of $Z$ on the plasma boundary over all grid points.
 \item \texttt{aspect} \\
       Type: \texttt{double} \\
       Aspect ratio, i.e., ratio of major radius over minor radius, of the plasma boundary.
 \item \texttt{betatotal} \\
       Type: \texttt{double} \\
       Total plasma beta.
 \item \texttt{betapol} \\
       Type: \texttt{double} \\
       Polodial plasma beta.
 \item \texttt{betator} \\
       Type: \texttt{double} \\
       Torodial plasma beta.
 \item \texttt{betaxis} \\
       Type: \texttt{double} \\
       Plasma beta on the magnetic axis.
 \item \texttt{volavgB} \\
       Type: \texttt{double} \\
       Volume-averaged magnetic field strength.
 \item \texttt{IonLarmor} \\
       Type: \texttt{double} \\
       Larmor radius of plasma ions.
 \item \texttt{ctor} \\
       Type: \texttt{double} \\
       Net toroidal plasma current.
 \item \texttt{rbtor} \\
       Type: \texttt{double} \\
       Poloidal ribbon current in the plasma at the boundary,
       i.e., sum of all (assumed, in case of a fixed-boundary run) coil currents
       plus the net poloidal plasma current.
 \item \texttt{rbtor0} \\
       Type: \texttt{double} \\
       Poloidal ribbon current in the plasma at the axis,
       i.e., sum of all (assumed, in case of a fixed-boundary run) coil currents,
       but without the poloidal plasma current.
\end{itemize}
%
A number of radial profiles of plasma parameters are available
in the next set of members in the \texttt{wout} group:
\begin{itemize}
 \item \texttt{iotaf} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Rotational transform~$\iota$ on the full-grid.
 \item \texttt{q\_factor} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Safety factor~$1/\iota$ on the full-grid.
 \item \texttt{presf} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Kinetic pressure~$p$ on the full-grid.
 \item \texttt{phi} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Enclosed toroidal magnetic flux~$\phi$ on the full-grid.
 \item \texttt{phipf} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Radial derivative of enclosed toroidal magnetic flux~$\phi'$ on the full-grid.
 \item \texttt{chi} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Enclosed poloidal magnetic flux~$\chi$ on the full-grid.
 \item \texttt{chipf} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Radial derivative of enclosed toroidal magnetic flux~$\chi'$ on the full-grid.
 \item \texttt{jcuru} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Radial derivative of enclosed poloidal current on full-grid.
 \item \texttt{jcurv} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Radial derivative of enclosed toroidal current on full-grid.
 \item \texttt{iotas} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Rotational transform~$\iota$ on the half-grid.
 \item \texttt{mass} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Plasma mass profile~$m$ on half-grid.
 \item \texttt{pres} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Kinetic pressure~$p$ on the half-grid.
 \item \texttt{beta\_vol} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Flux-surface averaged plasma beta on half-grid.
 \item \texttt{buco} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Profile of encosed toroidal current~$I$ on half-grid.
 \item \texttt{bvco} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Profile of encosed poloidal ribbon current~$G$ on half-grid.
 \item \texttt{vp} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Differential volume~$V'$ on half-grid.
 \item \texttt{equif} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Radial force balance residual on full-grid.
 \item \texttt{specw} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Spectral width~$M$ on full-grid.
 \item \texttt{phips} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       Radial derivative of enclosed toroidal magnetic flux~$\phi'$ on the half-grid.
 \item \texttt{over\_r} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns - 1]} \\
       $\langle \tau / R \rangle / V'$ on half-grid
 \item \texttt{jdotb} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       $\mathbf{j} \cdot \mathbf{B}$ on full-grid.
 \item \texttt{bdotb} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       $\mathbf{B} \cdot \mathbf{B}$ on full-grid.
 \item \texttt{bdotgradv} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Flux-surface-averaged toroidal magnetic field component $\mathbf{B} \cdot \nabla \zeta$ on full-grid.
 \item \texttt{DMerc} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Full Mercier stability criterion on the full-grid.
 \item \texttt{DShear} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Mercier stability criterion contribution due to magnetic shear on the full-grid.
 \item \texttt{DWell} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Mercier stability criterion contribution due to magnetic well on the full-grid.
 \item \texttt{DCurr} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Mercier stability criterion contribution due to plasma currents on the full-grid.
 \item \texttt{DGeod} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ns]} \\
       Mercier stability criterion contribution due to geodesic curvature on the full-grid.
\end{itemize}
%
The Fourier representation of the magnetic axis geometry
is the next set of members in the \texttt{wout} group:
\begin{itemize}
 \item \texttt{raxis\_cc}, \texttt{raxis\_cs} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ntor + 1]} \\
       Fourier coefficients of~$R(\zeta)$ of the magnetic axis geometry.
 \item \texttt{zaxis\_cs}, \texttt{zaxis\_cc} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[ntor + 1]} \\
       Fourier coefficients of~$Z(\zeta)$ of the magnetic axis geometry.
\end{itemize}
%
The state vector if VMEC++, i.e., the Fourier representation of the flux surface geometry ($R$ and $Z$)
and the stream function~$\lambda$
form the next set of members in the \texttt{wout} group:
\begin{itemize}
 \item \texttt{mnmax} \\
       Type: \texttt{int} \\
       Number of Fourier coefficients for the state vector.
 \item \texttt{xm} \\
       Type: \texttt{int[]} \\
       Dimensions: \texttt{[mnmax]} \\
       Poloidal mode numbers~$m$ for the Fourier coefficients in the state vector.
 \item \texttt{xn} \\
       Type: \texttt{int[]} \\
       Dimensions: \texttt{[mnmax]} \\
       Toroidal mode numbers times number of toroidal field periods~$n * n_\textrm{fp}$
       for the Fourier coefficients in the state vector.
 \item \texttt{rmnc}, \texttt{rmns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns][mnmax]} \\
       $R$ of the geometry of the flux surfaces on the full-grid.
 \item \texttt{zmns}, \texttt{zmnc} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns][mnmax]} \\
       $Z$ of the geometry of the flux surfaces on the full-grid.
 \item \texttt{lmns\_full}, \texttt{lmnc\_full} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns][mnmax]} \\
       $\lambda$ stream function on the full-grid.
 \item \texttt{lmns}, \texttt{lmnc} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns][mnmax]} \\
       $\lambda$ stream function on the half-grid (for backwards-compatibility with ancient versions of Fortran VMEC).
\end{itemize}
%
The Fourier representation of
the Jacobian and the magnetic field on the flux surfaces
form the next set of members in the \texttt{wout} group:
\begin{itemize}
 \item \texttt{mnyq} \\
       Type: \texttt{int} \\
       Number of poloidal Fourier modes used to describe the Nyquist-quantites,
       i.e., the Jacobian, the magnetic field strength, and the magnetic field components.
 \item \texttt{nnyq} \\
       Type: \texttt{int} \\
       Number of toroidal Fourier modes used to describe the Nyquist-quantites,
       i.e., the Jacobian, the magnetic field strength, and the magnetic field components.
 \item \texttt{mnmax\_nyq} \\
       Type: \texttt{int} \\
       Number of Fourier coefficients for the Nyquist-quantites,
       i.e., the Jacobian, the magnetic field strength, and the magnetic field components.
 \item \texttt{xm} \\
       Type: \texttt{int[]} \\
       Dimensions: \texttt{[mnmax\_nyq]} \\
       Poloidal mode numbers~$m$ for the Fourier coefficients in the Nyquist-quantites,
       i.e., the Jacobian, the magnetic field strength, and the magnetic field components.
 \item \texttt{xn} \\
       Type: \texttt{int[]} \\
       Dimensions: \texttt{[mnmax\_nyq]} \\
       Toroidal mode numbers times number of toroidal field periods~$n * n_\textrm{fp}$
       for the Fourier coefficients in the Nyquist-quantites,
       i.e., the Jacobian, the magnetic field strength, and the magnetic field components.
 \item \texttt{gmnc}, \texttt{gmns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the Jacobian~$\sqrt{g}$ of the coordinate transform between flux coordinates
       and cylindrical coordinates on the half-grid.
 \item \texttt{bmnc}, \texttt{bmns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the magnetic field strength~$\vert\mathbf{B}\vert$ on the half-grid.
 \item \texttt{bsubumnc}, \texttt{bsubumns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the covariant magnetic field component~$B_\theta$ on the half-grid.
 \item \texttt{bsubvmnc}, \texttt{bsubvmns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the covariant magnetic field component~$B_\zeta$ on the half-grid.
 \item \texttt{bsubsmns}, \texttt{bsubsmnc} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns][mnmax\_nyq]} \\
       Fourier coefficients of the covariant magnetic field component~$B_s$ on the full-grid.
 \item \texttt{bsupumnc}, \texttt{bsupumns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the contravariant magnetic field component~$B^\theta$ on the half-grid.
 \item \texttt{bsupvmnc}, \texttt{bsupvmns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the contravariant magnetic field component~$B^\zeta$ on the half-grid.
 \item \texttt{currumnc}, \texttt{currumns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns][mnmax\_nyq]} \\
       Fourier coefficients of the contravariant current density component~$j^\theta$ on the full-grid.
       Not implemented yet.
 \item \texttt{currvmnc}, \texttt{currvmns} \\
       Type: \texttt{double[][]} \\
       Dimensions: \texttt{[ns - 1][mnmax\_nyq]} \\
       Fourier coefficients of the contravariant current density component~$j^\zeta$ on the full-grid.
       Not implemented yet.
\end{itemize}
%
In case of a free-boundary run,
the Fourier representation of
the vacuum scalar magnetic potential on the plasma boundary
is the last set of members in the \texttt{wout} group:
\begin{itemize}
 \item \texttt{mnmaxpot} \\
       Type: \texttt{int} \\
       Number of Fourier modes used to represent
       the vacuum scalar magnetic potential on the boundary.
       Not implemented yet.
 \item \texttt{xmpot} \\
       Type: \texttt{int[]} \\
       Dimensions: \texttt{[mnmaxpot]} \\
       Poloidal mode numbers~$m$ for the Fourier coefficients
       of the vacuum scalar magnetic potential.
       Not implemented yet.
 \item \texttt{xnpot} \\
       Type: \texttt{int[]} \\
       Dimensions: \texttt{[mnmaxpot]} \\
       Toroidal mode numbers times number of toroidal field periods~$n * n_\textrm{fp}$
       for the Fourier coefficients
       of the vacuum scalar magnetic potential.
       Not implemented yet.
 \item \texttt{potsin}, \texttt{potcos} \\
       Type: \texttt{double[]} \\
       Dimensions: \texttt{[mnmaxpot]} \\
       Fourier coefficients of the vacuum scalar magnetic potential on the plasma boundary.
       Not implemented yet.
\end{itemize}
