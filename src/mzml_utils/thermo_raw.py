"""Thermo RAW access through Thermo's RawFileReader, run under Mono.

Only what the spectrum cache needs: every MS1 scan's centroids, read from the RAW file itself.

Why it exists. A conversion preset that keeps only the N most intense peaks per scan (msconvert's
``threshold count N most-intense``) truncates every busy MS1 survey, while MS2 scans rarely reach
the cap. Everything that reads MS1 peaks from such an mzML -- precursor isotope checks,
co-isolation, XICs, MS1 quantification -- then reads censored data. The RAW file still holds every
centroid, so the fix is to read MS1 from it and patch the spectrum cache
(``spectrum_cache.patch_ms1``) instead of re-converting.

The centroids are the instrument's own label peaks, *including* reference and exception peaks
(``GetCentroidStream(scan, true)``). That is what ProteoWizard's "Thermo/Xcalibur peak picking"
writes: on an Orbitrap Eclipse run, every capped mzML MS1 peak was found in this export at 0 ppm
and identical intensity, and every capped scan held exactly the export's N most intense peaks.
Leaving reference/exception peaks out drops about 4% of them, mostly shoulders beside intense ions.

The Thermo DLLs are not redistributed here. ``find_pwiz_dir`` looks for
``ThermoFisher.CommonCore.{Data,RawFileReader}.dll`` in, in order: an explicit directory,
``$THERMO_PWIZ_DIR`` (one directory, or several separated by ``os.pathsep``), a local copy made on
first use (``~/.cache/mzml_utils/thermo/dlls``), and a ThermoRawFileParser install under
``~/.local/opt/ThermoRawFileParser`` (it bundles the same two DLLs). A ProteoWizard install has
them too. The first hit off the local copy is copied into it, so later runs never read the DLLs
from a network share. ``mono`` and ``mcs`` (or ``csc``) must be on PATH.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import struct
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np

DLLS = ("ThermoFisher.CommonCore.Data.dll", "ThermoFisher.CommonCore.RawFileReader.dll")

#: Local installs that bundle the DLLs, searched last.
LOCAL_DLL_DIRS = (str(Path.home() / ".local" / "opt" / "ThermoRawFileParser"),)


def _cache_root() -> Path:
    return Path(os.environ.get("MZML_UTILS_CACHE", Path.home() / ".cache" / "mzml_utils"))


def _has_dlls(d: Path) -> bool:
    return all((d / name).is_file() for name in DLLS)


def find_pwiz_dir(user_value: Optional[str] = None, required: bool = True) -> Optional[Path]:
    """Directory holding the two RawFileReader DLLs, or None (``required=False``).

    Search order: ``user_value``, ``$THERMO_PWIZ_DIR``, the local copy, ``LOCAL_DLL_DIRS``. A hit
    outside the local copy is copied into it so the next call does not touch the network.
    """
    local = _cache_root() / "thermo" / "dlls"
    env = [d for d in os.environ.get("THERMO_PWIZ_DIR", "").split(os.pathsep) if d]
    candidates = [user_value, *env, str(local), *LOCAL_DLL_DIRS]
    for c in candidates:
        if not c:
            continue
        d = Path(c)
        try:
            if not _has_dlls(d):
                continue
        except OSError:          # an unreachable share
            continue
        if d != local:
            try:
                local.mkdir(parents=True, exist_ok=True)
                for name in DLLS:
                    if not (local / name).is_file() or (local / name).stat().st_size != (d / name).stat().st_size:
                        shutil.copy2(d / name, local / name)
            except OSError:
                return d
        return local if _has_dlls(local) else d
    if required:
        raise FileNotFoundError(
            "Thermo RawFileReader DLLs not found (ThermoFisher.CommonCore.Data.dll, "
            "ThermoFisher.CommonCore.RawFileReader.dll). Pass pwiz_dir= / --pwiz-dir or set "
            "THERMO_PWIZ_DIR.")
    return None


def raw_reader_available(pwiz_dir: Optional[str] = None) -> bool:
    """True when mono, a C# compiler and the DLLs are all present."""
    return bool(shutil.which("mono") and (shutil.which("mcs") or shutil.which("csc"))
                and find_pwiz_dir(pwiz_dir, required=False))


# The exporter. Binary output, little-endian:
#   b"MS1X", int32 version
#   per MS1 scan: int32 scan, float64 rt_min, int32 source, int32 n, float64[n] mz, float64[n] intensity
#   int32 -1
# source: 1 = centroid stream (label peaks incl. reference/exception), 2 = centroid scan data,
#         3 = profile data (no centroid stream; NOT what msconvert's vendor peak picking writes).
_MS1_EXPORT_CS = r'''
using System;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

public static class Ms1Export
{
    public static int Main(string[] args)
    {
        if (args.Length < 3) { Console.Error.WriteLine("usage: Ms1Export <dll_dir> <raw> <out.bin>"); return 2; }
        var dllDir = args[0];
        AppDomain.CurrentDomain.AssemblyResolve += (s, e) =>
        {
            var p = Path.Combine(dllDir, new AssemblyName(e.Name).Name + ".dll");
            return File.Exists(p) ? Assembly.LoadFrom(p) : null;
        };
        return Run(args[1], args[2]);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    static int Run(string raw, string outPath)
    {
        using (var rf = RawFileReaderAdapter.FileFactory(raw))
        {
            if (rf == null || !rf.IsOpen || rf.IsError) { Console.Error.WriteLine("cannot open " + raw); return 3; }
            rf.SelectInstrument(Device.MS, 1);
            int first = rf.RunHeaderEx.FirstSpectrum, last = rf.RunHeaderEx.LastSpectrum;
            int n = 0;
            using (var bw = new BinaryWriter(File.Create(outPath)))
            {
                bw.Write(new byte[] { (byte)'M', (byte)'S', (byte)'1', (byte)'X' });
                bw.Write(1);
                for (int s = first; s <= last; s++)
                {
                    var f = rf.GetFilterForScanNumber(s);
                    if (f.MSOrder != MSOrderType.Ms) continue;
                    double rt = rf.RetentionTimeFromScanNumber(s);
                    double[] mz, it; int src;
                    var cs = rf.GetCentroidStream(s, true);
                    if (cs != null && cs.Length > 0) { mz = cs.Masses; it = cs.Intensities; src = 1; }
                    else
                    {
                        var st = rf.GetScanStatsForScanNumber(s);
                        var seg = rf.GetSegmentedScanFromScanNumber(s, st);
                        mz = seg.Positions; it = seg.Intensities; src = st.IsCentroidScan ? 2 : 3;
                    }
                    bw.Write(s); bw.Write(rt); bw.Write(src); bw.Write(mz.Length);
                    foreach (var v in mz) bw.Write(v);
                    foreach (var v in it) bw.Write(v);
                    n++;
                }
                bw.Write(-1);
            }
            Console.WriteLine("ms1_scans\t" + n + "\tfirst\t" + first + "\tlast\t" + last);
        }
        return 0;
    }
}
'''

EXPORTER_VERSION = hashlib.sha256(_MS1_EXPORT_CS.encode()).hexdigest()[:12]


@dataclass
class RawMs1Scan:
    scan: int
    rt: float              # minutes
    source: int            # see the exporter comment above
    mz: np.ndarray         # float64
    intensity: np.ndarray  # float64


def compile_ms1_exporter(pwiz_dir: Optional[str] = None) -> Path:
    """Compile the exporter once per source version; returns the .exe (cached)."""
    dlls = find_pwiz_dir(pwiz_dir)
    compiler = shutil.which("mcs") or shutil.which("csc")
    if not shutil.which("mono") or not compiler:
        raise RuntimeError("mono and mcs (or csc) are required to read Thermo RAW files")
    out_dir = _cache_root() / "thermo" / f"ms1export_{EXPORTER_VERSION}"
    exe = out_dir / "Ms1Export.exe"
    if exe.is_file():
        return exe
    out_dir.mkdir(parents=True, exist_ok=True)
    src = out_dir / "Ms1Export.cs"
    src.write_text(_MS1_EXPORT_CS)
    refs = [f"-r:{dlls / name}" for name in DLLS]
    tmp = out_dir / "Ms1Export.exe.tmp"
    proc = subprocess.run([compiler, "-nologo", f"-out:{tmp}", *refs, str(src)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"compiling the MS1 exporter failed:\n{proc.stdout}\n{proc.stderr}")
    os.replace(tmp, exe)
    return exe


def read_ms1x(path: str) -> Dict[int, RawMs1Scan]:
    """Read the exporter's binary file into {scan: RawMs1Scan}."""
    out: Dict[int, RawMs1Scan] = {}
    with open(path, "rb") as fh:
        if fh.read(4) != b"MS1X":
            raise ValueError(f"{path}: not an MS1 export")
        (version,) = struct.unpack("<i", fh.read(4))
        if version != 1:
            raise ValueError(f"{path}: unsupported export version {version}")
        while True:
            head = fh.read(4)
            if len(head) < 4:
                raise ValueError(f"{path}: truncated (no terminator)")
            (scan,) = struct.unpack("<i", head)
            if scan < 0:
                break
            rt, src, n = struct.unpack("<dii", fh.read(16))
            mz = np.frombuffer(fh.read(8 * n), dtype="<f8")
            it = np.frombuffer(fh.read(8 * n), dtype="<f8")
            if len(mz) != n or len(it) != n:
                raise ValueError(f"{path}: truncated in scan {scan}")
            out[scan] = RawMs1Scan(scan, rt, src, mz, it)
    return out


def write_ms1x(path: str, scans: Iterable[RawMs1Scan]) -> None:
    """Write scans in the exporter's format (tests and hand-built fixtures)."""
    with open(path, "wb") as fh:
        fh.write(b"MS1X" + struct.pack("<i", 1))
        for s in scans:
            mz = np.asarray(s.mz, dtype="<f8")
            it = np.asarray(s.intensity, dtype="<f8")
            fh.write(struct.pack("<idii", int(s.scan), float(s.rt), int(s.source), len(mz)))
            fh.write(mz.tobytes())
            fh.write(it.tobytes())
        fh.write(struct.pack("<i", -1))


def sha256_file(path: str, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def export_ms1(raw_path: str, pwiz_dir: Optional[str] = None, out_path: Optional[str] = None,
               timeout: Optional[float] = None) -> Dict[int, RawMs1Scan]:
    """Every MS1 scan's centroids from a Thermo RAW file, as {scan: RawMs1Scan}.

    The RAW path is resolved first: RawFileReader can open a symlink yet fail to select the MS
    device. ``out_path`` keeps the binary export; by default a temp file is used and removed.
    """
    raw = os.path.realpath(raw_path)
    if not os.path.isfile(raw):
        raise FileNotFoundError(raw)
    exe = compile_ms1_exporter(pwiz_dir)
    dlls = find_pwiz_dir(pwiz_dir)
    keep = out_path is not None
    if out_path is None:
        fd, out_path = tempfile.mkstemp(suffix=".ms1x")
        os.close(fd)
    try:
        proc = subprocess.run(["mono", str(exe), str(dlls), raw, out_path],
                              capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            raise RuntimeError(f"MS1 export failed for {raw} (exit {proc.returncode}):\n{proc.stderr}")
        return read_ms1x(out_path)
    finally:
        if not keep and os.path.exists(out_path):
            os.remove(out_path)
