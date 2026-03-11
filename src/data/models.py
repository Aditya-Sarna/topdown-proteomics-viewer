from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np


@dataclass
class Spectrum:
    scan_id: str
    mz_array: np.ndarray
    intensity_array: np.ndarray
    retention_time: float = 0.0
    precursor_mz: float = 0.0
    precursor_charge: int = 0
    precursor_mass: float = 0.0
    ms_level: int = 2
    file_name: str = ""

    def normalize(self) -> np.ndarray:
        max_i = self.intensity_array.max()
        if max_i > 0:
            return self.intensity_array / max_i * 100.0
        return self.intensity_array.copy()

    def to_dict(self) -> dict:
        return {
            'scan_id': self.scan_id,
            'mz': self.mz_array.tolist(),
            'intensity': self.intensity_array.tolist(),
            'retention_time': self.retention_time,
            'precursor_mz': self.precursor_mz,
            'precursor_charge': self.precursor_charge,
            'precursor_mass': self.precursor_mass,
            'ms_level': self.ms_level,
            'file_name': self.file_name,
        }

    @staticmethod
    def from_dict(d: dict) -> 'Spectrum':
        return Spectrum(
            scan_id=d['scan_id'],
            mz_array=np.array(d['mz']),
            intensity_array=np.array(d['intensity']),
            retention_time=d.get('retention_time', 0.0),
            precursor_mz=d.get('precursor_mz', 0.0),
            precursor_charge=d.get('precursor_charge', 0),
            precursor_mass=d.get('precursor_mass', 0.0),
            ms_level=d.get('ms_level', 2),
            file_name=d.get('file_name', ''),
        )


@dataclass
class Modification:
    position: int        # 1-based position in the form sequence
    name: str
    mass_shift: float
    residue: str = ""

    def to_dict(self) -> dict:
        return {'position': self.position, 'name': self.name,
                'mass_shift': self.mass_shift, 'residue': self.residue}

    @staticmethod
    def from_dict(d: dict) -> 'Modification':
        return Modification(position=d['position'], name=d['name'],
                            mass_shift=d['mass_shift'], residue=d.get('residue', ''))


@dataclass
class Proteoform:
    sequence: str
    protein_name: str = ""
    start_pos: int = 1
    end_pos: int = -1
    modifications: List[Modification] = field(default_factory=list)
    theoretical_mass: float = 0.0
    observed_mass: float = 0.0
    mass_error_da: float = 0.0
    mass_error_ppm: float = 0.0
    score: float = 0.0
    matched_ions: int = 0
    total_ions: int = 0

    def __post_init__(self):
        if self.end_pos < 0:
            self.end_pos = len(self.sequence)

    def to_dict(self) -> dict:
        return {
            'sequence': self.sequence,
            'protein_name': self.protein_name,
            'start_pos': self.start_pos,
            'end_pos': self.end_pos,
            'modifications': [m.to_dict() for m in self.modifications],
            'theoretical_mass': self.theoretical_mass,
            'observed_mass': self.observed_mass,
            'mass_error_da': self.mass_error_da,
            'mass_error_ppm': self.mass_error_ppm,
            'score': self.score,
            'matched_ions': self.matched_ions,
            'total_ions': self.total_ions,
        }

    @staticmethod
    def from_dict(d: dict) -> 'Proteoform':
        pf = Proteoform(
            sequence=d['sequence'],
            protein_name=d.get('protein_name', ''),
            start_pos=d.get('start_pos', 1),
            end_pos=d.get('end_pos', -1),
            modifications=[Modification.from_dict(m) for m in d.get('modifications', [])],
            theoretical_mass=d.get('theoretical_mass', 0.0),
            observed_mass=d.get('observed_mass', 0.0),
            mass_error_da=d.get('mass_error_da', 0.0),
            mass_error_ppm=d.get('mass_error_ppm', 0.0),
            score=d.get('score', 0.0),
            matched_ions=d.get('matched_ions', 0),
            total_ions=d.get('total_ions', 0),
        )
        return pf


@dataclass
class FragmentIon:
    ion_type: str      # b, y, c, z, a
    position: int
    charge: int
    mz: float
    mass: float
    sequence: str = ""
    matched: bool = False
    observed_mz: float = 0.0
    mass_error_da: float = 0.0
    mass_error_ppm: float = 0.0

    def label(self) -> str:
        sup = {2: '²', 3: '³', 4: '⁴', 5: '⁵'}
        return f"{self.ion_type}{self.position}{sup.get(self.charge, '')}"

    def to_dict(self) -> dict:
        return {
            'ion_type': self.ion_type, 'position': self.position,
            'charge': self.charge, 'mz': self.mz, 'mass': self.mass,
            'sequence': self.sequence, 'matched': self.matched,
            'observed_mz': self.observed_mz,
            'mass_error_da': self.mass_error_da,
            'mass_error_ppm': self.mass_error_ppm,
        }

    @staticmethod
    def from_dict(d: dict) -> 'FragmentIon':
        return FragmentIon(**d)


@dataclass
class Feature:
    feature_id: str
    mz_apex: float
    rt_apex: float
    rt_start: float
    rt_end: float
    mz_start: float
    mz_end: float
    intensity: float
    charge: int
    monoisotopic_mass: float = 0.0
    sequence: str = ""
    proteoform_id: str = ""

    def to_dict(self) -> dict:
        return self.__dict__

    @staticmethod
    def from_dict(d: dict) -> 'Feature':
        return Feature(**d)


@dataclass
class SearchResult:
    proteoform: Proteoform
    fragment_ions: List[FragmentIon] = field(default_factory=list)
    matched_b: int = 0
    matched_y: int = 0
    matched_c: int = 0
    matched_z: int = 0
    sequence_coverage: float = 0.0

    def to_dict(self) -> dict:
        return {
            'proteoform': self.proteoform.to_dict(),
            'fragment_ions': [ion.to_dict() for ion in self.fragment_ions],
            'matched_b': self.matched_b,
            'matched_y': self.matched_y,
            'matched_c': self.matched_c,
            'matched_z': self.matched_z,
            'sequence_coverage': self.sequence_coverage,
        }

    @staticmethod
    def from_dict(d: dict) -> 'SearchResult':
        return SearchResult(
            proteoform=Proteoform.from_dict(d['proteoform']),
            fragment_ions=[FragmentIon.from_dict(i) for i in d.get('fragment_ions', [])],
            matched_b=d.get('matched_b', 0),
            matched_y=d.get('matched_y', 0),
            matched_c=d.get('matched_c', 0),
            matched_z=d.get('matched_z', 0),
            sequence_coverage=d.get('sequence_coverage', 0.0),
        )
