"""Utilities for drawing two-dimensional flow boundaries."""

import matplotlib.pyplot as plt
import numpy as np

import aerokit.aero.degree as deg
from aerokit.aero.plot.defaultstyle import sty_wall, sty_wall_fill


class Wall:
    """A polyline or single-valued curved boundary.

    One coordinate may be a callable and the other a two-value range or an
    array of sample points.  ``location`` identifies the solid side of the
    boundary and may be ``"top"``, ``"bottom"``, ``"left"`` or ``"right"``.
    """

    _locations = {"top", "bottom", "left", "right"}

    def __init__(self, x, y, location=None, npts=101, fill_style=None, **style):
        if callable(x) and callable(y):
            raise ValueError("only one of x and y may be callable")
        if location is not None and location not in self._locations:
            raise ValueError("location must be one of {}".format(self._locations))
        self._xdef = x
        self._ydef = y
        self.location = location
        self.npts = npts
        self.style = {**sty_wall, **style}
        self.fill_style = {**sty_wall_fill, **({} if fill_style is None else fill_style)}

    def _samples(self, values):
        values = np.asarray(values)
        if values.ndim != 1:
            raise ValueError("boundary coordinates must be one-dimensional")
        if values.size == 2:
            return np.linspace(values[0], values[1], self.npts)
        return values

    @property
    def x(self):
        if callable(self._xdef):
            return np.asarray(self._xdef(self.y))
        return self._samples(self._xdef)

    @property
    def y(self):
        if callable(self._ydef):
            return np.asarray(self._ydef(self.x))
        return self._samples(self._ydef)

    def angle(self, x=None):
        """Return the local boundary angle, in degrees."""
        if x is None:
            return deg.atan(np.gradient(self.y, self.x))
        x = np.asarray(x)
        isegment = np.searchsorted(self.x, x, side="right") - 1
        isegment = np.clip(isegment, 0, self.x.size - 2)
        slope = np.diff(self.y)[isegment] / np.diff(self.x)[isegment]
        return deg.atan(slope)

    def plot(self, ax=None, fill=True):
        """Draw the boundary and, when requested, its solid side."""
        ax = plt.gca() if ax is None else ax
        if fill and self.location is not None:
            self._fill(ax)
        style = {'zorder': 3, **self.style}
        return ax.plot(self.x, self.y, **style)

    def _fill(self, ax):
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        if self.location == "bottom":
            x = np.r_[self.x, self.x[-1], self.x[0]]
            y = np.r_[self.y, ymin, ymin]
        elif self.location == "top":
            x = np.r_[self.x, self.x[-1], self.x[0]]
            y = np.r_[self.y, ymax, ymax]
        elif self.location == "left":
            x = np.r_[self.x, xmin, xmin]
            y = np.r_[self.y, self.y[-1], self.y[0]]
        else:  # right
            x = np.r_[self.x, xmax, xmax]
            y = np.r_[self.y, self.y[-1], self.y[0]]
        return ax.fill(x, y, **self.fill_style)


class Geom:
    """A collection of :class:`Wall` objects defining a flow geometry."""

    def __init__(self, walls=None):
        self.walls = [] if walls is None else list(walls)

    def add_wall(self, x, y, location=None, npts=101, fill_style=None, **style):
        """Create, store, and return a boundary wall."""
        wall = Wall(x, y, location=location, npts=npts, fill_style=fill_style, **style)
        self.walls.append(wall)
        return wall

    def subplots(self, figsize=(14, 8), xlim=None, ylim=None, ax=None):
        """Create or configure an equal-aspect, axis-free geometry plot."""
        if ax is None:
            self.fig, self.ax = plt.subplots(figsize=figsize, facecolor="white")
        else:
            self.ax = ax
            self.fig = ax.figure
        self.ax.set(aspect="equal")
        self.ax.axis("off")
        if xlim is not None:
            self.ax.set_xlim(xlim)
        if ylim is not None:
            self.ax.set_ylim(ylim)
        return self.fig, self.ax

    def plot(self, ax=None, fill=True):
        """Draw all boundaries and return their Matplotlib line artists."""
        ax = getattr(self, "ax", None) if ax is None else ax
        if ax is None:
            _, ax = self.subplots()
        return [line for wall in self.walls for line in wall.plot(ax=ax, fill=fill)]
