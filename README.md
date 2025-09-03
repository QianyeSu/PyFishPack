# PyFishPack - A python wrapper for FishPack

Here we need a logo or banner

Here are some badges (pipy, conda, doi, unit tests, docs, downloads, etc.)

## 1. Introduction
 
**F**ast **I**terative **S**olver for **H**elmholtz equations, a.k.a **FishPack**, is a collection of Fortran programs and subroutines that solve 2nd- and 4th-order finite difference approximations to separable elliptic Partial Differential Equations (PDEs) developed by NCAR.  Here is [a link to the FishPack](https://github.com/NCAR/NCAR-Classic-Libraries-for-Geophysics).

These legacy codes, although robust and efficient, is less accessible to mordern programming language like **Python** and its ecosystem.  Here we develop this project to

- make these codes much easier to use for pythoners, and
- provide some examples in solving the PDE problems in atmospheric science and physical oceanography.

---

## 2. How to install

**Requirements**
`PyFishPack` is developed under the environment with `xarray` (=version 0.15.0), `dask` (=version 2.11.0), `numpy` (=version 1.15.4).  Older versions of these packages are not well tested.

**Install via pip**
```
pip install PyFishPack (not yet finished)
```

**Install via conda**
```
conda install -c conda-forge PyFishPack (not yet finished)
```

**Install from github**
```
git clone https://github.com/QianyeSu/PyFishPack.git
cd PyFishPack
python setup.py install (not yet finished)
```

---
## 3. Equations

Naming conventions for a subroutine contains four parts, for example:
> [t][h][ws][crt]

- [t] or None, means test or the main program for the subroutine **hwscrt**;
- [h] means Helmholtz;
- [ws | st | w3] means [Arakawa-A | Arakawa-C staggered | third-order] grids;
- [crt | plr | cyl | ssp | csp] means [cartesian | polor | cylindric | spherical | spherical] coordinates;


hwscrt and hstcrt

>$$
\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \lambda u = f(x,y)
$$

hwsplr and hstplr

>$$
\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial u}{\partial r}\right) + \frac{1}{r^2}\frac{\partial^2 u}{\partial \theta^2} + \lambda u = f(r,\theta)
$$

hwscyl and hstcyl

>$$
\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial u}{\partial r}\right) + \frac{\partial^2 u}{\partial z^2} + \frac{\lambda}{r^2} u = f(r,z)
$$


hwsssp

>$$
\frac{1}{\sin\theta}\frac{\partial}{\partial \theta}\left(\sin\theta\frac{\partial u}{\partial \theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2 u}{\partial \phi^2} + \lambda u = f(\theta,\phi)
$$


hwscsp

>$$
\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial u}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial \theta}\left(\sin\theta\frac{\partial u}{\partial \theta} \right)+ \frac{\lambda}{r\sin^2\theta} u = f(\theta,r)
$$


hstcsp

>$$
\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial u}{\partial r}\right) + \frac{1}{\sin\theta}\frac{\partial}{r\partial\theta}\left(\sin\theta\frac{\partial u}{r\partial \theta}\right) + \frac{\lambda}{r\sin^2\theta} u = f(\theta,r)
$$


hw3crt

>$$
\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2} + \lambda u = f(x,y,z)
$$


sepx4

>$$
a(x)\frac{\partial^2 u}{\partial x^2} + b(x)\frac{\partial u}{\partial x} + c(x)u + \frac{\partial^2 u}{\partial y^2} = g(x,y)
$$

sepeli

>$$
a(x)\frac{\partial^2 u}{\partial x^2} + b(x)\frac{\partial u}{\partial x} + c(x)u + d(y)\frac{\partial^2 u}{\partial y^2} + e(y)\frac{\partial u}{\partial y} + f(y)u = g(x,y)
$$


These include Helmholtz equations in cartesian, polar, cylindrical, and spherical coordinates, as well as more general separable elliptic equations. The solvers use the cyclic reduction algorithm. When the problem is singular, a least-squares solution is computed. Singularities induced by the coordinate system are handled, including at the origin *r=0* in cylindrical coordinates, and at the poles in spherical coordinates.

Test programs are provided for the 19 solvers in the `examples` folder. Each serves two purposes: as a template to guide you in writing your own codes utilizing the fishpack solvers, and as a demonstration that you can correctly produce the executables. 

---

## 4. Quick examples


