"""recorded_precursor_ms1: the MS1 the instrument recorded for an MSn scan, on both readers.

A hand-written five-scan mzML: MS1 1 and 2, then MS2 scans whose precursor element records MS1 1
(not the nearest MS1, 2), nothing, a non-MS1 scan, and a scan that does not exist."""
import base64
import os

import numpy as np

from mzml_utils import spectrum_cache as sc
from mzml_utils.reader import MzMLReader


def _arr(values):
    return base64.b64encode(np.asarray(values, dtype="<f8").tobytes()).decode()


def _spectrum(idx, scan, level, ref=None):
    nid = f"controllerType=0 controllerNumber=1 scan={scan}"
    mz, it = [400.0, 500.0], [10.0, 20.0]
    prec = ""
    if level == 2:
        ref_attr = f' spectrumRef="controllerType=0 controllerNumber=1 scan={ref}"' if ref else ""
        prec = (f'<precursorList count="1"><precursor{ref_attr}><selectedIonList count="1"><selectedIon>'
                '<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="500.0"/>'
                '<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2"/>'
                '</selectedIon></selectedIonList><activation>'
                '<cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>'
                '</activation></precursor></precursorList>')
    filt = "FTMS + p NSI Full ms" if level == 1 else "FTMS + c NSI d Full ms2 500.0@hcd30.00"
    arrays = "".join(
        f'<binaryDataArray encodedLength="{len(_arr(v))}">'
        '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>'
        '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>'
        f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value=""/><binary>{_arr(v)}</binary></binaryDataArray>'
        for v, acc, name in ((mz, "MS:1000514", "m/z array"), (it, "MS:1000515", "intensity array")))
    return (f'<spectrum index="{idx}" id="{nid}" defaultArrayLength="2">'
            f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>'
            f'<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" '
            f'value="{scan / 10:.2f}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>'
            f'<cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="{filt}"/></scan></scanList>'
            f'{prec}<binaryDataArrayList count="2">{arrays}</binaryDataArrayList></spectrum>')


def _write(path):
    specs = [(1, 1, None), (2, 1, None), (3, 2, 1), (4, 2, None), (5, 2, 3), (6, 2, 99)]
    body = "".join(_spectrum(i, s, lv, ref) for i, (s, lv, ref) in enumerate(specs))
    path.write_text('<?xml version="1.0" encoding="utf-8"?><mzML xmlns="http://psi.hupo.org/ms/mzml" '
                    f'version="1.1.0"><run id="r"><spectrumList count="{len(specs)}">{body}'
                    '</spectrumList></run></mzML>')


EXPECTED = {3: (1, "recorded"), 4: (2, "nearest_preceding:no_recorded_ref"),
            5: (2, "nearest_preceding:recorded_ref_not_ms1"), 6: (2, "nearest_preceding:recorded_ref_missing"),
            42: (None, "scan_missing")}


def test_mzml_reader(tmp_path):
    p = tmp_path / "run.mzML"
    _write(p)
    with MzMLReader(str(p)) as r:
        assert {s: r.recorded_precursor_ms1(s) for s in EXPECTED} == EXPECTED


def test_cache_agrees_with_the_mzml(tmp_path):
    p = tmp_path / "run.mzML"
    _write(p)
    db = sc.cache_path_for(str(p))
    os.makedirs(os.path.dirname(db), exist_ok=True)
    sc.build_cache(str(p), db, local_copy=False)
    with sc.SpectrumCache(db, mirror_local=False) as c:
        assert {s: c.recorded_precursor_ms1(s) for s in EXPECTED} == EXPECTED
