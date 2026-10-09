"""Physical fields from an explicitly selected independent contact eigenmode."""
import json
from pathlib import Path
import numpy as np
from .continuum_contact_mode import ContinuumContactMode
from .finite_mass_contact_reference import FiniteMassContact
from .implicit_plasma_records import sha


class ArchivedContactMode(ContinuumContactMode):
    def __init__(self, directory='outputs/radial16_reference', index=1):
        directory=Path(directory);summary=directory/'summary.json'
        record=json.loads(summary.read_text())['continuum'][index]
        path=directory/record['arrays_file']
        if sha(path)!=record['arrays_sha256']:raise ValueError('Unchanged independent continuum arrays required')
        a=dict(np.load(path));vector=a['leading_vector']
        if record['leading_imaginary']!=0 or np.max(abs(vector.imag))!=0:
            raise ValueError('A real independent growing mode is required')
        self.input_sha256=sha(path);self.input_summary_sha256=sha(summary)
        self.directory=str(directory);self.index=index;self.radial_mode=record['radial_mode']
        self.ref=FiniteMassContact(velocity=0.);self.q=float(a['q'])
        self.growth=record['maximum_real'];self.points=len(a['z'])
        vector=vector.real/np.max(abs(vector.real))
        if vector[np.argmax(abs(vector))]<0:vector=-vector
        physical=(vector*a['scale']).reshape(9,self.points)
        self.coefficients=np.fft.fft(physical,axis=1)/self.points
        self.frequency=2*np.pi*np.fft.fftfreq(self.points,self.ref.length/self.points)
        self.frequency[self.points//2]=0.
        self.derivative=a['derivative'];self.curl=a['curl'];self.physical_nodes=physical
