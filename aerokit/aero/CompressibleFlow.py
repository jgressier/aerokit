"""Removed legacy compatibility module."""

raise ImportError(
    "aerokit.aero.CompressibleFlow is no longer available. "
    "Import the required functions from aerokit.aero.Isentropic, "
    "aerokit.aero.Supersonic, or aerokit.aero.MassFlow."
)
