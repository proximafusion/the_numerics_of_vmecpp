\chapter{BNORM}

\FloatBarrier
\section{Introduction}
Every current is surrounded by an accompanying magnetic field.
A toroidal plasma cannot confine itself by its own currents.
Thus, external confinement coils are required to produce the (mostly toroidal) confinement magnetic field.
In a Stellarator, these confinement coils are required (in VMEC: assumed)
to produce good flux surfaces in the confimenent region.
The diamagnetic current~$\mathbf{j}_\perp$ inside the plasma arises from a finite plasma pressure:
\begin{equation}
 \nabla p = \mathbf{j}_\perp \times \mathbf{B} \, .
\end{equation}
The current paths inside the plasma are closed: $\nabla \cdot \mathbf{j} = 0$.
This imples $\nabla \cdot \mathbf{j}_\perp = -\nabla \cdot \mathbf{j}_{||}$,
where $\mathbf{j} = \mathbf{j}_\perp + \mathbf{j}_{||}$ is the total current density in the plasma
and $\mathbf{j}_{||}$ is the Pfirsch-Schlüter current density.
The diamagnetic current is a current flowing inside the plasma itself.
It has its own magnetic field associated with it,
which is superposed to the confinement magnetic field,
thereby contributing to the total magnetic field inside the plasma domain.
VMEC can be used to compute the total magnetic field inside the plasma domain
such that flux surfaces are in force balance.

It is a non-trivial problem to separate the individual contributions (plasma current and coil currents)
to the total magnetic field from each other.
The magnetic field due to the confinement coils (``external'' magnetic field)
is typically computed using the Biot-Savart law
for a given geometry of the coils and given coil currents along those current paths.
It can be computed in all of space; both inside the plasma domain
as well as outside of it and also on the boundary.

The total magnetic field at the plasma boundary is assumed to be purely tangential to the boundary.
This implies that the normal component of the magnetic field contribution from the coils
is cancelled by the normal component of the magnetic field due to plasma currents.
A Stellarator coil configuration must have flux surfaces in vacuum alone;
otherwise, it is not a suitable Stellarator magnetic configuration.
The magnetic field on each of these vacuum flux surfaces
is purely tangential to the flux surfaces by definition.
However, if a plasma is created now inside the flux surface region
and pressure builds up inside the plasma,
the magnetic field due to coils inside the plasma region is expelled from the plasma region.
This is the so-called diamagnetic behavior of the plasma.
Good flux surfaces are assumed to persist even in case of a finite plasma pressure in VMEC.
In fact, stochastization of the magnetic field inside the plasma (breaking up of flux surfaces)
poses a limit on the achievable plasma pressure. SPEC can be used to compute these limits,
but this is not the topic of this particular note.
The buildup of plasma pressure implies that currents start to flow inside the plasma.
These currents necessarily change the magnetic field structure inside the plasma as well as outside of it.
Assuming that good flux surfaces still exist,
the flux surfaces change shape and shift around in space (e.g. the Shafranov shift is one of these outcomes).
Suddenly, the vacuum flux surfaces (on which the external magnetic field is purely tangential)
do not coincide with the plasma flux surfaces anymore.
This implies that a non-zero normal component of the plasma magnetic field is present
and that (because of the assumption on the continued existence of flux surfaces in the plasma)
only the new total magnetic field is purely tangential to the new plasma boundary.

Inside the plasma domain, the magnetic field due to plasma currents
can be obtained by subtracting the coil field (available from Biot-Savart)
from the total magnetic field (available from VMEC).
The currents inside the plasma can be computed using Ampere's law:
\begin{equation}
 \mu_0 \mathbf{j} = \nabla \times \mathbf{B} \, .
\end{equation}
A covariant representation of the total magnetic field inside the plasma is available in VMEC
and this can be used to compute the contravariant representation of the current density inside the plasma.
This can then be used to compute the magnetic field due to plasma currents anywhere in space.
Of particular interest is the computation of the magnetic field due to plasma currents
outside the plasma domain, where it is not available trivially.
Here, a Biot-Savart (three-dimensional) volume integral over the plasma current density
can be used to compute the magnetic field outside the plasma domain.
This approach is commonly used to predict magnetic diagnostic signals.

The virtual casing principle can be used to reduce the volume integral
to a surface integral over the plasma boundary
by replacing the volume current density inside the plasma domain
with a surface current density flowing inside the plasma boundary surface,
which generates the same magnetic field as the volume plasma currents
outside the plasma domain.
This method is available in, e.g., DIAGNO.
The surface integral is nearly-singular when evaluating the external field due to plasma currents
close to the surface. This is the reason why an adaptive integration routine
is used in DIAGNO to compute the predictions for the magnetic diagnostics
(which are typically located quite close to the LCFS).

BNORM computes the normal component of the external magnetic field via the virtual casing integral \emph{on} the LCFS.
This is the normal component of the magnetic field due to plasma currents alone.
The integral is singular in this case, since the evaluation location is located
in the same surface as where the surface current density (the source term of the Biot-Savart integral) is defined.
A regularization method is therefore required.

\FloatBarrier
\section{Theoretical Foundation}
Assume an arbitrary magnetic field $\mathbf{B}_0$ is given in a toroidal domain.
In case of a VMEC computation, this is the total magnetic field
due to both external confinement coils and the plasma currents.
Furthermore, $\mathbf{B}_0$ shall be purely tangential to the boundary of the domain,
such that $\mathbf{B}_0 \cdot \mathbf{n} = 0$, where $\mathbf{n}$ is the exterior normal of the boundary surface.
The virtual-casing principle allows to define a surface current $\mathbf{j}$ flowing on the boundary
that is given by:
\begin{equation}
  \mu_0 \mathbf{j} = \mathbf{B}_0 \times \mathbf{n}
\end{equation}
which leads to the same magnetic field outside the toroidal domain
as the current source terms inside the volume of the domain.
When computing the magnetic field outside the toroidal domain,
one can replace the Biot-Savart volume integral over the current within the domain
by a surface integral over above surface current density.
The magnetic vector potential $\mathbf{A}(\mathbf{x})$ outside the domain
due to the surface current is given by the Biot-Savart law,
which takes on the form of a two-dimensional surface integral in this case:
\begin{equation}
 \mathbf{A}(\mathbf{x}) = \frac{\mu_0}{4 \pi} \iint \frac{\mathbf{j}(\mathbf{x}')}{|\mathbf{x} - \mathbf{x}'|} \,\mathrm{d}^2\mathbf{x}' \, .
\end{equation}
The magnetic field outside the toroidal domain
is then given by $\mathbf{B} = \nabla \times \mathbf{A}$.
The normal component of this magnetic field~$B_\mathrm{n}$ is then given by
\begin{equation}
 B_\mathrm{n} = -\mathbf{n} \cdot \mathbf{B} = -\mathbf{n} \cdot (\nabla \times \mathbf{A}) \label{eqn:b_norm_def} \, .
\end{equation}
This is the main quantity of interest in this computation.
The outward-directed surface normal~$\mathbf{n}$ is given by:
\begin{equation}
 \mathbf{n} = - \frac{\mathbf{x}_u \times \mathbf{x}_v}{|\mathbf{x}_u \times \mathbf{x}_v|} \, .
\end{equation}
Inserting this into Eqn.~(\ref{eqn:b_norm_def}), it follows:
\begin{align}
 B_\mathrm{n}
 =&\ \frac{1}{|\mathbf{x}_u \times \mathbf{x}_v|}
     \left( \mathbf{x}_u \times \mathbf{x}_v \right) \cdot (\nabla \times \mathbf{A}) \nonumber \\
 =&\ \frac{1}{|\mathbf{x}_u \times \mathbf{x}_v|}
     \left(  \mathbf{x}_v \cdot \frac{\partial \mathbf{A}}{\partial u}
           - \mathbf{x}_u \cdot \frac{\partial \mathbf{A}}{\partial v} \right) \nonumber \\
 =&\ \frac{1}{|\mathbf{x}_u \times \mathbf{x}_v|}
     \left[  \frac{\partial}{\partial u} \left(\mathbf{x}_v \cdot \mathbf{A} \right)
           - \frac{\partial}{\partial v} \left(\mathbf{x}_u \cdot \mathbf{A} \right) \right] \, .
\end{align}
In above formula, the following trick was used:
\begin{align}
 ~&\,   \frac{\partial}{\partial u} \left(\mathbf{x}_v \cdot \mathbf{A} \right)
      - \frac{\partial}{\partial v} \left(\mathbf{x}_u \cdot \mathbf{A} \right) \nonumber \\
 =&\,   \mathbf{x}_v \cdot \frac{\partial \mathbf{A}}{\partial u} \bcancel{+ \mathbf{x}_{uv} \cdot \mathbf{A}}
      - \mathbf{x}_u \cdot \frac{\partial \mathbf{A}}{\partial v} \bcancel{- \mathbf{x}_{uv} \cdot \mathbf{A}} \nonumber \\
 =&\,   \mathbf{x}_v \cdot \frac{\partial \mathbf{A}}{\partial u}
      - \mathbf{x}_u \cdot \frac{\partial \mathbf{A}}{\partial v} \, .
\end{align}
The derivatives are conveniently computed by Fourier-transforming
$\mathbf{x}_v \cdot \mathbf{A}$ and $\mathbf{x}_u \cdot \mathbf{A}$
and then computing the tangential derivatives analytically.
The remainder of this note deals with the computation of the vector potential.
Have a closer look at the current density:
\begin{align}
 \mu_0 |\mathbf{x}_u \times \mathbf{x}_v| \,\mathbf{j}
 =&\, \mathbf{B}_0 \times \left( \mathbf{x}_u \times \mathbf{x}_v \right) \nonumber \\
 =&\,   \mathbf{x}_u \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_v \right)
      - \mathbf{x}_v \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_u \right)
\end{align}
according to the BAC-CAB rule.
The Jacobian from realspace to normalized angle-like coordinates on the surface is:
\begin{equation}
 \mathrm{d}^2\mathbf{x} = |\mathbf{x}_u \times \mathbf{x}_v| \,\mathrm{d}u \,\mathrm{d}v \, .
\end{equation}
It follows for the vector potential~$\mathbf{A}$:
\begin{equation}
 \mathbf{A}(u, v) = \frac{1}{4 \pi} \iint
    \frac{  \mathbf{x}_{u'} \cdot \left( \mathbf{B}_0' \cdot \mathbf{x}_{v'} \right)
          - \mathbf{x}_{v'} \cdot \left( \mathbf{B}_0' \cdot \mathbf{x}_{u'} \right) }
         {|\mathbf{x}(u,v) - \mathbf{x}(u',v')|}
    \,\mathrm{d}u' \,\mathrm{d}v' \, .
\end{equation}
This function is evaluated on a grid in $u$ and $v$,
which is the same grid on which the integral over $u'$ and $v'$ is computed,
since this allows to re-use certain geometric quantities.
However, this approach introduces a singularity in the integrand.
The distance $|\mathbf{x}(u,v) - \mathbf{x}(u',v')|$ in the denominator of the integrand
vanishes when $u = v'$ simulatenous with $v = v'$, leading to a $1/x$-type singularity in the integrand.
A subtraction method is applied to solve this issue.
An function with the same singularity, which can be Fourier-transformed analytically,
is subtracted and its analytical Fourier transform is added back in afterwards.
The regularizing function approximates the actual integrand:
\begin{equation}
  \frac{  \mathbf{x}_{u'} \cdot \left( \mathbf{B}_0' \cdot \mathbf{x}_{v'} \right)
        - \mathbf{x}_{v'} \cdot \left( \mathbf{B}_0' \cdot \mathbf{x}_{u'} \right) }
       {|\mathbf{x}(u,v) - \mathbf{x}(u',v')|}
  \approx
  \frac{  \mathbf{x}_{u} \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_{v} \right)
        - \mathbf{x}_{v} \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_{u} \right)}
       {\left[ g_{uu} \delta u^2 + 2 g_{uv} \delta u \delta v + g_{vv} \delta v^2 \right]^{1/2} }
\end{equation}
where
\begin{align}
 g_{uu} =&\, \mathbf{x}_{u} \cdot \mathbf{x}_{u} \\
 g_{uv} =&\, \mathbf{x}_{u} \cdot \mathbf{x}_{v} \\
 g_{vv} =&\, \mathbf{x}_{v} \cdot \mathbf{x}_{v} \\
 \delta u =&\, \frac{1}{\pi} \tan\left(\pi (u - u') \right) \\
 \delta v =&\, \frac{1}{\pi} \tan\left(\pi (v - v') \right) \, .
\end{align}
The analytical Fourier transform of this equivalently-singular function is:
\begin{equation}
 I(a,b,c) = \pi \int\limits_0^1 \int\limits_0^1
   \frac{\mathrm{d}u \,\mathrm{d}v}
        {\left[ a \tan^2(\pi u) + 2 b \tan(\pi u) \tan(\pi v) + c \tan^2(\pi v) \right]^{1/2} } \, .
\end{equation}
The analytical Fourier transform results in:
\begin{equation}
 I(a,b,c) = T_{0}^{+} + T_{0}^{-}
\end{equation}
with
\begin{equation}
 T_{0}^{\pm} = \frac{1}{\sqrt{a \pm 2b + c}} \log\left(
   \frac{\sqrt{c \left(a \pm 2b + c\right)} + c \pm b}
        {\sqrt{a \left(a \pm 2b + c\right)} - a \mp b} \right) \, .
\end{equation}
The vector potential is then computed as follows:
\begin{equation}
 \mathbf{A}(u,v) = \mathbf{A}_\textrm{reg}(u,v) + \mathbf{A}_\textrm{sing}(u,v)
\end{equation}
with
\begin{align}
 \mathbf{A}_\textrm{reg}(u,v)
 =&\, \frac{1}{4 \pi} \iint \,\mathrm{d}u' \,\mathrm{d}v' \Biggl\{
      \frac{  \mathbf{x}_{u'} \cdot \left( \mathbf{B}_0' \cdot \mathbf{x}_{v'} \right)
            - \mathbf{x}_{v'} \cdot \left( \mathbf{B}_0' \cdot \mathbf{x}_{u'} \right) }
           {|\mathbf{x}(u,v) - \mathbf{x}(u',v')|} \nonumber \\
 ~&\,   - \frac{\pi \left[  \mathbf{x}_{u} \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_{v} \right)
                      - \mathbf{x}_{v} \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_{u} \right) \right] }
           {\left[ g_{uu} \delta u^2 + 2 g_{uv} \delta u \delta v + g_{vv} \delta v^2 \right]^{1/2} }
         \Biggr\} \label{eqn:a_sing}
\end{align}
and
\begin{equation}
 \mathbf{A}_\textrm{sing}(u,v)
 = \left[   \mathbf{x}_{u} \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_{v} \right)
          - \mathbf{x}_{v} \cdot \left( \mathbf{B}_0 \cdot \mathbf{x}_{u} \right) \right]
   I(g_{uu}, g_{uv}, g_{vv}) \, .
\end{equation}
The discrete Fourier transform in Eqn.~(\ref{eqn:a_sing}) is performed
numerically, where the singular point is not considered in the summation.


% \FloatBarrier
% \section{Numerical Implementation}



\FloatBarrier
\section{Further Considerations}
The source term of the normal component of the magnetic field computed by BNORM
is the diamagnetic current density, which arises due to a finite plasma pressure.
Thus, for a vacuum computation, the resulting normal field should vanish,
since the source term is missing in case of a zero-pressure/vacuum computation.

% TODO: run VMEC for a vanishing-pressure case and run BNORM on the output; expected outcome: B.n vanishes as well.

BNORM computes B.n due to plasma currents;
NESTOR computes $|B|^2/2$ due to plasma currents.

BNORM plays a central role in the classical two-stage Stellarator optimization approach.
The plasma is optimized subject to physics objective functions using fixed-boundary VMEC.
BNORM is then used to compute the normal component of the magnetic field due to plasma currents
on the LCFS. This components needs to be cancelled by the external magnetic field.
A jump can be present in the tangential magnetic field components across the LCFS
between vacuum and plasma, as long as the total magnitude of the magnetic field
is conserved between plasma and vacuum.
The normal component of the magnetic field due to plasma currents
is therefore (along with the surface geometry) input to the coil design code
and defines the target to match/cancel.

% Workflow:
% \begin{enumerate}
%  \item fixed-boundary optimization using fixed-bdy VMEC
%  \item BNORM: get normal component of magnetic field on LCFS
%  \item NESCOIL/... to design coils for given BNORM output and LCFS geometry from VMEC
% \end{enumerate}
