\section{NESTOR for Tokamaks}
This is a variant of NESTOR, which makes explicit use of axisymmetry
to carry out the toroidal integrals analytically.

The grid steps are:
\begin{align}
 \Delta_u      =&\, \frac{1}{n_{\theta,1}} \\
 \Delta_\theta =&\, 2 \pi \Delta_u \, .
\end{align}
%
The Fourier basis functions can be evaluated
for a fixed grid in advance.
Here are the regular basis functions:
\begin{align}
 c^{lm} =&\, \cos(\Delta_\theta (m+1) l) \Delta_u \\
 s^{lm} =&\, \sin(\Delta_\theta (m+1) l) \Delta_u \, .
\end{align}
Here are the poloidal derivative basis functions:
\begin{align}
 c_m^{lm} =&\, 2 \pi (m+1) \cos(\Delta_\theta (m+1) l) \Delta_u \\
 s_m^{lm} =&\, 2 \pi (m+1) \sin(\Delta_\theta (m+1) l) \Delta_u \, .
\end{align}
%
The Fourier basis functions are used
to compute the realspace geometry
of the last closed flux surface
as provided by VMEC:
\begin{align}
 R(\theta_l) = \sum\limits_{m = 0}^{M}
   \hat{R}^\textrm{cos}_m \cos(m l \Delta_\theta) \\
 Z(\theta_l) = \sum\limits_{m = 0}^{M}
   \hat{Z}^\textrm{cos}_m \sin(m l \Delta_\theta)
\end{align}
as well as the poloidal derivatives:
\begin{align}
 R_u(\theta_l) = \sum\limits_{m = 0}^{M}
   (- 2 \pi m) \hat{R}^\textrm{cos}_m \sin(m l \Delta_\theta) \\
 Z_u(\theta_l) = \sum\limits_{m = 0}^{M}
   (  2 \pi m) \hat{Z}^\textrm{cos}_m \cos(m l \Delta_\theta)
\end{align}
plus, in case of a non-stellarator-symmetric
(respectively non-up/down-symmetric) case,
the other sin/cos terms.
%
The metric element~$g_{\theta \theta}$ is given by:
\begin{equation}
 g_{\theta \theta} = R_u^2 + Z_u^2 \, .
\end{equation}
%
The cylindrical components of the external magnetic field
$B^R(\theta_l)$, $B^\varphi(\theta_l)$ and $B^Z(\theta_l)$
are interpolated from the \texttt{mgrid} file
at the circumference of the plasma boundary.
%
The normal component of the magnetic field is given by:
\begin{equation}
 (\mathbf{B} \cdot \mathbf{n})(\theta_l) =
 2 \pi R(\theta_l) \left[
   Z_u(\theta_l) B^R(\theta_l) - R_u(\theta_l) B^Z(\theta_l)
 \right] \, .
\end{equation}
%
The total poloidal current~$G$
is furthermore computed as:
\begin{equation}
 G = \frac{2 \pi}{n_\theta} \sum\limits_{l = 0}^{n_\theta - 1}
   R(\theta_l) B^\varphi(\theta_l)
\end{equation}
which is simply an average over the plasma boundary:
\begin{equation}
 G = 2 \pi \langle R B^\varphi \rangle \, .
\end{equation}
%
The normal magnetic field component
is Fourier-transformed:
\begin{align}
 \hat{B}_{n,\textrm{ext}}^s =&\,
   \sum\limits_{l = 0}^{n_\theta - 1}
     (\mathbf{B} \cdot \mathbf{n})(\theta_l) s^{l m} \\
 \hat{B}_{n,\textrm{ext}}^c =&\,
   \sum\limits_{l = 0}^{n_\theta - 1}
     (\mathbf{B} \cdot \mathbf{n})(\theta_l) c^{l m} \, .
\end{align}
%
Define the following shorthand notation:
\begin{align}
 \sigma^R_{l l'} =&\, R(\theta_l) + R(\theta_{l'}) \\
 \delta^R_{l l'} =&\, R(\theta_l) - R(\theta_{l'}) \\
 \delta^R_{l l'} =&\, Z(\theta_l) - Z(\theta_{l'})
\end{align}
as well as:
\begin{align}
 a^{+}_{l l'} =&\, \sqrt{(\sigma^R_{l l'})^2 + (\delta^Z_{l l'})^2} \\
 a^{-}_{l l'} =&\, \sqrt{(\delta^R_{l l'})^2 + (\delta^Z_{l l'})^2} \, .
\end{align}
Then define
\begin{equation}
 w_{l l'} = \frac{a^{+}_{l l'} - a^{-}_{l l'}}{a^{+}_{l l'} + a^{-}_{l l'}}
\end{equation}
which leads to
\begin{align}
 x_{l l'} =&\, 1 - \frac{4 w_{l l'}}{(1 + w_{l l'})^2} \\
          =&\, \frac{(1 - w_{l l'})^2}{(1 + w_{l l'})^2} \, .
\end{align}
This is used as argument to complete elliptic integrals to compute:
\begin{equation}
 E_{l l'} = \frac{1}{\pi w_{l l'}} \left[
   \frac{1 + w_{l l'}^2}{1 + w_{l l'}} \mathcal{K}(x_{l l'})
   - (1 + w_{l l'}) \mathcal{E}(x_{l l'})
 \right] \, .
\end{equation}
%
Moreover, it is computed:
\begin{equation}
 a_\textrm{sing}(l - l')
 = \log\left( \frac{1}{\pi} \sin\left( \frac{\pi}{n_\theta} (l - l') \right) \right) \, .
\end{equation}
%
Furthermore, we compute:
\begin{equation}
 A = \frac{1}{2} \left( a^{+}_{l l'} + a^{-}_{l l'} \right) \, .
\end{equation}
%
The off-diagonal elements of the interaction matrix in realspace
can now be computed as:
\begin{equation}
 G_{l,l',\textrm{reg}}
 =
 \pi R(\theta_l) R(\theta_{l'}) \frac{E_{l l'}}{A_{l l'}}
 + \frac{1}{2} \left[ R(\theta_l) + R(\theta_{l'}) \right] a_\textrm{sing}(l - l') \, .
\end{equation}
%
The diagonal elements of the interaction matrix in realspace are computed as:
\begin{equation}
 G_{l,l,\textrm{reg}}
 =
 R(\theta_l) \left(
   3 \log(2) - 2 - \log(\frac{\sqrt{g_{\theta \theta}}}{R(\theta_l)})
 \right) \, .
\end{equation}
%
The matrix is symmetric
and therefore, the upper triangle computed above
is mirrored into the lower triangle.
%
The Fourier transforms for forming the matrix blocks are computed now:
\begin{align}
  f^\texttt{cc}_{m,m'} =&\, f^\texttt{cc}_{m,m',\textrm{reg}} -  f^\texttt{cc}_{m,m',\textrm{sing}} \\
  f^\texttt{cs}_{m,m'} =&\, f^\texttt{cs}_{m,m',\textrm{reg}} -  f^\texttt{cs}_{m,m',\textrm{sing}} \\
  f^\texttt{ss}_{m,m'} =&\, f^\texttt{ss}_{m,m',\textrm{reg}} -  f^\texttt{ss}_{m,m',\textrm{sing}} \, .
\end{align}
%
The regular matrix block contributions are:
\begin{align}
 f^\texttt{cc}_{m,m',\textrm{reg}}
 =&\, \sum\limits_{l=0}^{n_\theta-1} \sum\limits_{l'=0}^{n_\theta-1}
   c_m^{lm} G_{l',l,\textrm{reg}} c_m^{l' m'} \\
 f^\texttt{cs}_{m,m',\textrm{reg}}
 =&\, \sum\limits_{l=0}^{n_\theta-1} \sum\limits_{l'=0}^{n_\theta-1}
   c_m^{lm} G_{l',l,\textrm{reg}} s_m^{l' m'} \\
 f^\texttt{ss}_{m,m',\textrm{reg}}
 =&\, \sum\limits_{l=0}^{n_\theta-1} \sum\limits_{l'=0}^{n_\theta-1}
   s_m^{lm} G_{l',l,\textrm{reg}} s_m^{l' m'} \, .
\end{align}
%
The singular matrix kernel is
\begin{equation}
 g^s_m = - \frac{n_\theta}{2 (m + 1)}
\end{equation}
and this is used to form:
\begin{equation}
 \hat{g}^s_{m m'} = \frac{1}{2} \left(g^s_m + g^s_{m'} \right) \, .
\end{equation}
%
This is used to form the singular matrix blocks:
\begin{align}
 f^\texttt{cc}_{m,m',\textrm{sing}}
 =&\, \sum\limits_{l=0}^{n_\theta-1}
   c_m^{lm} R(\theta_l) \hat{g}^s_{m m'} c_m^{l' m'} \\
 f^\texttt{cs}_{m,m',\textrm{sing}}
 =&\, \sum\limits_{l=0}^{n_\theta-1}
   c_m^{lm} R(\theta_l) \hat{g}^s_{m m'} s_m^{l' m'} \\
 f^\texttt{ss}_{m,m',\textrm{sing}}
 =&\, \sum\limits_{l=0}^{n_\theta-1}
   s_m^{lm} R(\theta_l) \hat{g}^s_{m m'} s_m^{l' m'} \, .
\end{align}
%
The right-hand-side vector of the linear system~$\mathbf{g}$
is computed as follows:
\begin{equation}
 \mathbf{g} =
 \begin{pmatrix}
  \hat{g}^s \\
  \hat{g}^c
 \end{pmatrix}
\end{equation}
with
\begin{align}
 \hat{g}^s =&\, \hat{g}_\textrm{reg}^s + \hat{g}_\textrm{sing}^s \\
 \hat{g}^c =&\, \hat{g}_\textrm{reg}^c - \hat{g}_\textrm{sing}^c \, .
\end{align}
In there, we use:
\begin{align}
 \hat{g}_{m,\textrm{reg}}^s
 =&\,
 \sum\limits_{l = 0}^{n_\theta}
   \sum\limits_{l' = 0}^{n_\theta}
     G_{l',l,\textrm{reg}} c_m^{lm} \Delta_u \\
 \hat{g}_{m,\textrm{reg}}^c
 =&\, -
 \sum\limits_{l = 0}^{n_\theta}
   \sum\limits_{l' = 0}^{n_\theta}
     G_{l',l,\textrm{reg}} s_m^{lm} \Delta_u \, .
\end{align}
%
The singular part is compute as follows:
\begin{align}
 \hat{g}_{m,\textrm{sing}}^s
 =&\,
 \sum\limits_{l = 0}^{n_\theta}
   \alpha_m R(\theta_l) s^{lm} \\
 \hat{g}_{m,\textrm{sing}}^c
 =&\,
 \sum\limits_{l = 0}^{n_\theta}
   \alpha_m R(\theta_l) c^{lm}
\end{align}
with
\begin{equation}
 \alpha_m =
   \frac{1}{2}
   \left[
     \pi + 2 \pi (m+1) \log(2 \pi)
   \right] \, .
\end{equation}
%
The matrix of the linear system~$\mathcal{H}$ is
now constructed from the blocks:
\begin{equation}
 \mathcal{H} =
 \begin{pmatrix}
   f^\texttt{cc}_{m,m'} & -f^\texttt{cs}_{m,m'} \\
  -f^\textrm{cs}_{m',m} &  f^\textrm{ss}_{m,m'}
 \end{pmatrix}
\end{equation}
%
The inverse of the matrix is obtained
by computing the Choleskey decomposition of the matrix
and then solving for a set of right-hand-side vectors
forming the equally-sized identity matrix.
%
The solution to the linear system is then
computed as follows:
\begin{equation}
 \begin{pmatrix}
  \hat{\Phi}^s \\
  \hat{\Phi}^c
 \end{pmatrix}
 =
 \mathcal{H}^{-1}
 \left[
 I_\textrm{tor}
 \begin{pmatrix}
  \hat{g}^s \\
  \hat{g}^c
 \end{pmatrix}
 +
 \begin{pmatrix}
  \hat{B}_{n,\textrm{ext}}^s \\
  \hat{B}_{n,\textrm{ext}}^c
 \end{pmatrix}
 \right]
\end{equation}
%
The potential derivative~$\Phi_u$ is computed
directly from the Fourier coefficients of the potential:
\begin{equation}
 \Phi_u(\theta_l)
 =
 I_\textrm{tor}
 - n_{\theta,1} \sum\limits_{m = 0}^{m_f}
     \left[
       \hat{\Phi}^s_m c_m^{lm} - \hat{\Phi}^c_m s_m^{lm}
     \right] \, .
\end{equation}
%
The vacuum magnetic pressure~$|\mathbf{B}_\textrm{vac}|^2/2$
is computed as follows:
\begin{equation}
 \frac{\vert\mathbf{B}_\textrm{vac}\vert^2}{2}
 =
 \frac{1}{2} \left[
   \frac{G^2}{(2 \pi R)^2}
    +
    \frac{\Phi_u^2}{g_{\theta \theta}^2}
 \right] \, .
\end{equation}
