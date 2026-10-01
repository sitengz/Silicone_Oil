#!/usr/bin/env python3
"""Rasterize the three SVG plot sources as baseline 300-dpi RGB TIFF files.

This optional preview exporter uses Linux librsvg/cairo via ctypes; MATLAB's
plot_pdms_series.m independently writes native TIFF and editable FIG files.
"""

from __future__ import annotations

import ctypes
import struct
import xml.etree.ElementTree as ET
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCALE = 2
DPI = 300


def render_rgb(path: Path) -> tuple[int, int, bytes]:
    svg_root = ET.parse(path).getroot()
    width = int(svg_root.attrib["width"]) * SCALE
    height = int(svg_root.attrib["height"]) * SCALE
    rsvg = ctypes.CDLL("librsvg-2.so.2")
    cairo = ctypes.CDLL("libcairo.so.2")
    gobject = ctypes.CDLL("libgobject-2.0.so.0")

    rsvg.rsvg_handle_new_from_file.argtypes = (ctypes.c_char_p, ctypes.c_void_p)
    rsvg.rsvg_handle_new_from_file.restype = ctypes.c_void_p
    rsvg.rsvg_handle_render_cairo.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
    rsvg.rsvg_handle_render_cairo.restype = ctypes.c_int
    cairo.cairo_image_surface_create.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int)
    cairo.cairo_image_surface_create.restype = ctypes.c_void_p
    cairo.cairo_create.argtypes = (ctypes.c_void_p,)
    cairo.cairo_create.restype = ctypes.c_void_p
    cairo.cairo_scale.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
    cairo.cairo_surface_flush.argtypes = (ctypes.c_void_p,)
    cairo.cairo_image_surface_get_stride.argtypes = (ctypes.c_void_p,)
    cairo.cairo_image_surface_get_stride.restype = ctypes.c_int
    cairo.cairo_image_surface_get_data.argtypes = (ctypes.c_void_p,)
    cairo.cairo_image_surface_get_data.restype = ctypes.c_void_p
    cairo.cairo_destroy.argtypes = (ctypes.c_void_p,)
    cairo.cairo_surface_destroy.argtypes = (ctypes.c_void_p,)
    gobject.g_object_unref.argtypes = (ctypes.c_void_p,)

    handle = rsvg.rsvg_handle_new_from_file(str(path).encode(), None)
    if not handle:
        raise RuntimeError(f"librsvg could not open {path}")
    surface = cairo.cairo_image_surface_create(0, width, height)  # ARGB32
    context = cairo.cairo_create(surface)
    try:
        cairo.cairo_scale(context, float(SCALE), float(SCALE))
        if rsvg.rsvg_handle_render_cairo(handle, context) != 1:
            raise RuntimeError(f"librsvg could not render {path}")
        cairo.cairo_surface_flush(surface)
        stride = cairo.cairo_image_surface_get_stride(surface)
        raw = ctypes.string_at(cairo.cairo_image_surface_get_data(surface),
                               stride * height)
    finally:
        cairo.cairo_destroy(context)
        cairo.cairo_surface_destroy(surface)
        gobject.g_object_unref(handle)

    rgb = bytearray(width * height * 3)
    for row in range(height):
        bgra = raw[row * stride:row * stride + width * 4]
        start = row * width * 3
        end = start + width * 3
        rgb[start:end:3] = bgra[2::4]
        rgb[start + 1:end:3] = bgra[1::4]
        rgb[start + 2:end:3] = bgra[0::4]
    return width, height, bytes(rgb)


def write_tiff(path: Path, width: int, height: int, rgb: bytes) -> None:
    if len(rgb) != width * height * 3:
        raise ValueError("RGB byte count does not match TIFF dimensions")
    entries = 12
    bits_offset = 8 + 2 + entries * 12 + 4
    xres_offset = bits_offset + 6
    yres_offset = xres_offset + 8
    pixels_offset = yres_offset + 8
    tags = (
        (256, 4, 1, width),             # ImageWidth
        (257, 4, 1, height),            # ImageLength
        (258, 3, 3, bits_offset),       # BitsPerSample = 8, 8, 8
        (259, 3, 1, 1),                 # Compression = none
        (262, 3, 1, 2),                 # PhotometricInterpretation = RGB
        (273, 4, 1, pixels_offset),     # StripOffsets
        (277, 3, 1, 3),                 # SamplesPerPixel
        (278, 4, 1, height),            # RowsPerStrip
        (279, 4, 1, len(rgb)),         # StripByteCounts
        (282, 5, 1, xres_offset),       # XResolution
        (283, 5, 1, yres_offset),       # YResolution
        (296, 3, 1, 2),                 # ResolutionUnit = inch
    )
    with path.open("wb") as handle:
        handle.write(b"II")
        handle.write(struct.pack("<HIH", 42, 8, entries))
        for tag in tags:
            handle.write(struct.pack("<HHII", *tag))
        handle.write(struct.pack("<I", 0))
        handle.write(struct.pack("<HHH", 8, 8, 8))
        handle.write(struct.pack("<II", DPI, 1))
        handle.write(struct.pack("<II", DPI, 1))
        handle.write(rgb)


def main() -> None:
    for mean in (16, 32, 64):
        source = HERE / f"SZ_N{mean}.svg"
        destination = HERE / f"SZ_N{mean}.tif"
        width, height, rgb = render_rgb(source)
        write_tiff(destination, width, height, rgb)
        print(f"{destination.name}: {width} x {height}, {DPI} dpi")


if __name__ == "__main__":
    main()
