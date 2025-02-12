\chapter{VMEC Inputs}

\section{Numerical Resolution and Symmetry Assumptions}

\begin{itemize}
 \item \texttt{lasym} \\
       Type: \texttt{bool} \\
       Default: \texttt{false} \\
       This flag controls the assumption of stellarator symmetry in the computation.
       If set to \texttt{true}, the additional Fourier coefficients in \eqn{r_full} and \eqn{z_full}
       are assumed throughout the computation instead of the reduced sets in \eqn{r_symm} and \eqn{z_symm}.

 \item \texttt{nfp} \\
       Type: \texttt{int} \\
       Range: $>0$ \\
       Default: 1 \\
       This is the number of field periods, $n_\textrm{fp}$, to make use of toroidal symmetry.

 \item \texttt{mpol} \\
       Type: \texttt{int} \\
       Range: $\geq 2$ \\
       Default: 6 \\
       This controls the poloidal resolution of the computation in Fourier space.
       The set of poloidal mode numbers taken into account is
       $\{0, 1, ..., (\texttt{mpol} - 1)\}$.
       The number of poloidal Fourier modes taken into account is thus \texttt{mpol}.

 \item \texttt{ntor} \\
       Type: \texttt{int} \\
       Range: $\geq 0$ \\
       Default: 0 \\
       This controls the toroidal resolution of the computation in Fourier space.
       The set of toroidal mode numbers taken into account is
       $\{-\texttt{ntor}, ..., -1, 0, 1, ..., \texttt{ntor}\}$.
       The number of toroidal Fourier modes taken into account is thus $2 \times \texttt{ntor} + 1$.
       If $0$ is specified, an axisymmetric (Tokamak) run is done.

 \item \texttt{ntheta} \\
       Type: \texttt{int} \\
       Range: $\geq 0$ \\
       Default: 0 \\
       This controls the poloidal resolution of the computation in real space,
       i.e., the number of grid points in the poloidal direction.
       A value of $0$ implies that it is automatically chosen as the minimally allowed number of poloidal grid points.
       Any given value is checked against the minimally allowed number of poloidal grid points,
       which is $2 \times \texttt{mpol} + 6$.

 \item \texttt{nzeta} \\
       Type: \texttt{int} \\
       Range: $\geq 0$ \\
       Default: 0 \\
       This controls the toroidal resolution of the computation in real space,
       i.e., the number of grid points per field period in the toroidal direction.
       A value of $0$ implies that it is automatically chosen as the minimally allowed number of toroidal grid points per period.
       Any given value is checked against the minimally allowed number of toroidal grid points,
       which is $2 \times \texttt{ntor} + 4$.
       In case of an axisymmetric calculation ($\texttt{ntor} = 0$),
       a value of $0$ is replaced with $1$, but a different positive number can also be taken into account.
       Think of this as being able to average over many poloidal cutplanes for an axisymmetric run.

\end{itemize}

\section{Multi-Grid Steps}

\begin{itemize}

 \item \texttt{ns\_array} \\
       Type: \texttt{vector<int>} \\
       Range: each entry $\ge 3$, each entry $\geq$ the previous one \\
       Default: $\{ 31 \}$ \\
       This controls the radial resolution, i.e., the number of flux surfaces,
       in a series of consecutive resolution steps.
       The equilibrium is solved for the first entry in \texttt{ns\_array},
       then interpolated onto the next radial resolution grid, solved there, and so on.

 \item \texttt{ftol\_array} \\
       Type: \texttt{vector<double>} \\
       Range: each entry $> 0$ \\
       Default: $\{ 10^{-10} \}$ \\
       This controls the convergence criterion for each multi-grid step.
       The sum of the raw, un-preconditioned, force residuals has to fall below
       the given entry in \texttt{ftol\_array} for that multi-grid step to be considered converged.

 \item \texttt{niter\_array} \\
       Type: \texttt{vector<int>} \\
       Range: each entry $> 0$ \\
       Default: $\{ 100 \}$ \\
       This controls the number of allowed iterations for each multi-grid step
       for trying to achieve the requested force tolerance levels in \texttt{ftol\_array}.

\end{itemize}

\section{Global Physics Parameters}

\begin{itemize}

 \item \texttt{phiedge} \\
       Type: \texttt{double} \\
       Default: $1.0$ \\
       This specifies the toroidal magnetic flux enclosed by the boundary.
       In a fixed-boundary computation, this implicitly defines the magnetic field strength.
       In a free-boundary computation, this mainly influences the area of the poloidal cross section of the boundary
       and thus sets the volume of the plasma.

 \item \texttt{ncurr} \\
       Type: \texttt{int} \\
       Range: 0 or 1 \\
       Default: 0 \\
       This flag selects whether the rotational transform profile ($\texttt{ncurr} = 0$)
       or the toroidal current profile ($\texttt{ncurr} = 1$)
       are used as the second radial profile (next to the mass or pressure profile)
       to define the equilibrium to compute.

\end{itemize}

\section{Profile of Mass or Pressure}

\begin{itemize}

 \item \texttt{pmass\_type} \\
       Type: \texttt{string} \\
       Range: one of $\{ \texttt{power\_series},
                         \texttt{gauss\_trunc},
                         \texttt{two\_lorentz},
                         \texttt{two\_power},
                         \texttt{two\_power\_gs}, \\
                         \texttt{akima\_spline},
                         \texttt{cubic\_spline},
                         \texttt{pedestal},
                         \texttt{rational},
                         \texttt{line\_segment} \}$ \\
       Default: \texttt{power\_series} \\
       This defines the parameterization used for specifying the mass or pressure profile.

 \item \texttt{am} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the expansion coefficients for a parametric mass or pressure profile.
       Their interpretation changes depending on the value of \texttt{pmass\_type}.
       The pressure profile would be specified in Pascals.

 \item \texttt{am\_aux\_s} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the knot locations for a spline-parameterized specification
       of the mass or pressure profile.
       The first entry should be $0.0$ and the last entry should be $1.0$.
       At least 4 entries need to be provided.
       Equal number of entries need to be specified here and in \texttt{am\_aux\_f}.

 \item \texttt{am\_aux\_f} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the knot values for a spline-parameterized specification
       of the mass or pressure profile.
       At least 4 entries need to be provided.
       Equal number of entries need to be specified here and in \texttt{am\_aux\_s}.
       The pressure profile would be specified in Pascals.

 \item \texttt{pres\_scale} \\
       Type: \texttt{double} \\
       Range: $\geq 0$ \\
       Default: $1.0$ \\
       This is a global scaling factor for the mass or pressure profile in Pascals.

 \item \texttt{gamma} \\
       Type: \texttt{double} \\
       Range: $\geq 0$ \\
       Default: $0.0$ \\
       This is the adiabatic index (also ratio of specific heats) $\gamma$.
       Specifying $0$ implies that the pressure profile is specified.
       For all other values, the mass profile is specified.

 \item \texttt{spres\_ped} \\
       Type: \texttt{double} \\
       Range: $0 < \texttt{spres\_ped} \leq 1$ \\
       Default: $1.0$ \\
       This specifies the radial position in $s$ of a pedestal in the pressure profile.
       Radially outside of this location, the pressure is assumed to be constant.

\end{itemize}

\section{(Initial Guess for) Rotational Transform Profile}

Note that the rotational transform profile can be specified in the input file
even in case of a constrained-toroidal-current computation
(where $\texttt{ncurr} = 1$ and the rotational transform profile is thus computed on-the-fly internally).
The rotational transform profile specified in the input file
is then only used as an initial guess for computing the radial profiles of the magnetic fluxes
in the first iteration.
The key takeaway here is that a weird rotational transform profile in the input file
can mess up a VMEC++ computation even if the input file specifies
to run VMEC++ in constrained-toroidal-current mode
(and the user thus expects not having to care about the rotational transform profile).
The suggested workflow here is to only specify the rotational transform profile in the input file
if a constrained-rotational-transform computation is supposed to be done (with $\texttt{ncurr} = 0$).
The option to specify an initial guess for the rotational transform profile
in case of a constrained-toroidal-current computation is then reserved for the case
where the default initial guess for the rotational transform profil ($\iota = 1$) does not work.

\begin{itemize}

 \item \texttt{piota\_type} \\
       Type: \texttt{string} \\
       Range: one of $\{ \texttt{power\_series},
                         \texttt{sum\_atan},
                         \texttt{akima\_spline},
                         \texttt{cubic\_spline}, \\
                         \texttt{rational},
                         \texttt{line\_segment},
                         \texttt{nice\_quadratic} \}$ \\
       Default: \texttt{power\_series} \\
       This defines the parameterization used for specifying the rotational transform profile.

 \item \texttt{ai} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the expansion coefficients for a parametric rotational transform profile.
       Their interpretation changes depending on the value of \texttt{piota\_type}.

 \item \texttt{ai\_aux\_s} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the knot locations for a spline-parameterized specification
       of the rotational transform profile.
       The first entry should be $0.0$ and the last entry should be $1.0$.
       At least 4 entries need to be provided.
       Equal number of entries need to be specified here and in \texttt{ai\_aux\_f}.

 \item \texttt{ai\_aux\_f} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the knot values for a spline-parameterized specification
       of the rotational transform profile.
       At least 4 entries need to be provided.
       Equal number of entries need to be specified here and in \texttt{ai\_aux\_s}.

\end{itemize}

\section{(Initial Guess for) Toroidal Current Profile}

\begin{itemize}

 \item \texttt{pcurr\_type} \\
       Type: \texttt{string} \\
       Range: one of $\{ \texttt{power\_series},
                         \texttt{power\_series\_i},
                         \texttt{gauss\_trunc},
                         \texttt{sum\_atan}, \\
                         \texttt{two\_power},
                         \texttt{two\_power\_gs},
                         \texttt{akima\_spline\_i},
                         \texttt{akima\_spline\_ip}, \\
                         \texttt{cubic\_spline\_i},
                         \texttt{cubic\_spline\_ip},
                         \texttt{pedestal},
                         \texttt{rational}, \\
                         \texttt{line\_segment\_i},
                         \texttt{line\_segment\_ip},
                         \texttt{sum\_cossq\_s},
                         \texttt{sum\_cossq\_sqrts},
                         \texttt{sum\_cossq\_s\_free} \}$ \\
       Default: \texttt{power\_series} \\
       This defines the parameterization used for specifying the toroidal current profile.
       Traditionally, the profile of the radial derivative of the enclosed toroidal current
       was specified. This is referred to as the ``I-prime'' profile.
       The radial profile of the enclosed toroidal current was then computed by numerical quadrature.
       Later, additional profile types were implemented that allow to specify the radial profile
       of the enclosed toroidal current directly.
       These are distinguished between with the suffixes \texttt{\_i} (enclosed current)
       and \texttt{\_ip} (radial derivative of enclosed current)
       in the names of the profile parameterizations.

 \item \texttt{ac} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the expansion coefficients for a parametric toroidal current profile.
       Their interpretation changes depending on the value of \texttt{pcurr\_type}.

 \item \texttt{ac\_aux\_s} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the knot locations for a spline-parameterized specification
       of the toroidal current profile.
       The first entry should be $0.0$ and the last entry should be $1.0$.
       At least 4 entries need to be provided.
       Equal number of entries need to be specified here and in \texttt{ac\_aux\_f}.

 \item \texttt{ac\_aux\_f} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       These are the knot values for a spline-parameterized specification
       of the toroidal current profile.
       At least 4 entries need to be provided.
       Equal number of entries need to be specified here and in \texttt{ac\_aux\_s}.

 \item \texttt{curtor} \\
       Type: \texttt{double} \\
       Default: $0.0$ \\
       This is the net toroidal current in Amperes.
       The toroidal current profile is scaled to yield a total net toroidal current as given here.

 \item \texttt{bloat} \\
       Type: \texttt{double} \\
       Range: $\texttt{bloat} > 1$ \\
       Default: $1.0$ \\
       This is a scalar factor to re-scale the radial scale for the toroidal current profile.
       This factor can only deviate from $1.0$ in case of a constrained-current computation,
       i.e., when $\texttt{ncurr} = 1$.

\end{itemize}

\section{Free-Boundary Parameters}

\begin{itemize}

 \item \texttt{lfreeb} \\
       Type: \texttt{bool} \\
       Default: \texttt{false} \\
       This flag controls whether VMEC++ is run in fixed-boundary mode (\texttt{false}) or in free-boundary mode (\texttt{true}).

 \item \texttt{mgrid\_file} \\
       Type: \texttt{string} \\
       Default: \texttt{"NONE"} \\
       This is the filename to the MGRID file (in NetCDF format) that contains the magnetic field response factors for the external coils.

 \item \texttt{extcur} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ \}$ \\
       This is the array of coil currents that determines the net external magnetic field as produced by the coil contributions in the MGRID file.

 \item \texttt{nvacskip} \\
       Type: \texttt{int} \\
       Range: $1 \leq \texttt{nvacskip} \leq 10$ \\
       Default: 1 \\
       This is the interval for performing a full update of the free-boundary computation.
       Values greater than 1 mean that the matrix in the Nestor module is re-used for the next \texttt{nvacskip} iterations
       after the one that computed it (such as the first iteration that the free-boundary computation was activated).

 % \item \texttt{free\_boundary\_method} \\
 %       Type: \texttt{string} \\
 %       Default: \texttt{"nestor"} \\
 %       This enum controls which method is used in VMEC++
 %       to compute the free-boundary force contribution.
 %       Options are \texttt{"nestor"} for the standard
 %       NEumann Solver for TORoidal systems
 %       and \texttt{"biest"} for the Boundary Integral Equation Solver
 %       for Toroidal systems.
 %       Note that BIEST currently only works
 %       with zero net toroidal current (\texttt{curtor = 0}).

\end{itemize}

\section{Tweaking Parameters}

\begin{itemize}

 \item \texttt{nstep} \\
       Type: \texttt{int} \\
       Range: $\geq 1$ \\
       Default: 10 \\
       This is the interval of iterations at which the convergence progress (force residuals, etc.) are printed to the screen.

 \item \texttt{aphi} \\
       Type: \texttt{vector<double>} \\
       Default: $\{ 1.0 \}$ \\
       This is a power series expansion, without a constant offset,
       of the remapping profile for changing the radial location of the flux surfaces
       as a function of the normalized toroidal flux.

 \item \texttt{delt} \\
       Type: \texttt{double} \\
       Range: $0 < \texttt{delt} \leq 1$ \\
       Default: $1.0$ \\
       This is the initial value for the artificial time step used in the iterative solver in VMEC++.
       The time step used along the computation is automatically reduced on-the-fly
       if a too-large time step lead to crossing of the flux surfaces.

 \item \texttt{tcon0} \\
       Type: \texttt{double} \\
       Range: $0 \leq \texttt{tcon0} \leq 1$ \\
       Default: $1.0$ \\
       This is a scaling factor that can be used to tune down the spectral condensation force contribution.

 \item \texttt{lforbal} \\
       Type: \texttt{bool} \\
       Default: \texttt{false} \\
       This is a flag to use a non-variational form of the MHD forces
       for the innermost flux surface (the first one outside the magnetic axis).
       It can sometimes help convergence to enable this option.

\end{itemize}

\section{Initial Guess for Magnetic Axis Geometry}

\begin{itemize}

 \item \texttt{raxis\_c} \\
       Type: \texttt{vector<double>} \\
       Range: At least 1 value needs to be given, at most entries up to $n = \texttt{ntor}$ are considered. \\
       Default: $\{ \}$ \\
       These are the stellarator-symmetric (cosine) expansion coefficients $\hat{R}_{n}^{\cos}$ for the $R$ coordinate of the magnetic axis.

 \item \texttt{zaxis\_s} \\
       Type: \texttt{vector<double>} \\
       Range: At most entries up to $n = \texttt{ntor}$ are considered. The first entry ($n = 0$) is ignored. \\
       Default: $\{ \}$ \\
       These are the stellarator-symmetric (sine) expansion coefficients $\hat{Z}_{n}^{\sin}$ for the $Z$ coordinate of the magnetic axis.

 \item \texttt{raxis\_s} \\
       Type: \texttt{vector<double>} \\
       Range: At most entries up to $n = \texttt{ntor}$ are considered. The first entry ($n = 0$) is ignored. \\
       Default: $\{ \}$ \\
       These are the non-stellarator-symmetric (sine) expansion coefficients $\hat{R}_{n}^{\sin}$ for the $R$ coordinate of the magnetic axis.

 \item \texttt{zaxis\_c} \\
       Type: \texttt{vector<double>} \\
       Range: At least 1 value needs to be given, at most entries up to $n = \texttt{ntor}$ are considered. \\
       Default: $\{ \}$ \\
       These are the non-stellarator-symmetric (cosine) expansion coefficients $\hat{Z}_{n}^{\cos}$ for the $Z$ coordinate of the magnetic axis.

\end{itemize}

\section{(Initial Guess for) Boundary Geometry}

\begin{itemize}

 \item \texttt{rbc} \\
       Type: \texttt{vector<double>[]} \\
       Default: $\{ \}$ \\
       These are the stellarator-symmetric (cosine) expansion coefficients $\hat{R}_{m, n}^{\cos}$ for the $R$ coordinate of the (initial guess for the) boundary.

 \item \texttt{zbs} \\
       Type: \texttt{vector<double>[]} \\
       Default: $\{ \}$ \\
       These are the stellarator-symmetric (sine) expansion coefficients $\hat{Z}_{m, n}^{\sin}$ for the $Z$ coordinate of the (initial guess for the) boundary.

 \item \texttt{rbs} \\
       Type: \texttt{vector<double>[]} \\
       Default: $\{ \}$ \\
       These are the non-stellarator-symmetric (sine) expansion coefficients $\hat{R}_{m, n}^{\sin}$ for the $R$ coordinate of the (initial guess for the) boundary.

 \item \texttt{zbc} \\
       Type: \texttt{vector<double>[]} \\
       Default: $\{ \}$ \\
       These are the non-stellarator-symmetric (cosine) expansion coefficients $\hat{Z}_{m, n}^{\cos}$ for the $Z$ coordinate of the (initial guess for the) boundary.

\end{itemize}
