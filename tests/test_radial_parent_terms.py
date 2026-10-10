"""Independent cylindrical curl and original continuum cancellation checks."""
import numpy as np
from examples.archived_contact_mode import ArchivedContactMode
from examples.radial_parent_terms import terms


def test_component_curls_against_independent_fourth_order_differences():
    m=ArchivedContactMode();r=np.array([.009,.017,.029]);z=np.array([-.083,-.014,.071])
    result=terms(m,r,z);hr=1e-6;hz=1e-5
    for name,value in result['electric'].items():
        dr=sum(c*terms(m,r+i*hr,z)['electric'][name] for i,c in [(-2,1),(-1,-8),(1,8),(2,-1)])/(12*hr)
        dz=sum(c*terms(m,r,z+i*hz)['electric'][name] for i,c in [(-2,1),(-1,-8),(1,8),(2,-1)])/(12*hz)
        negative_curl=np.stack((dz[...,1],dr[...,2]-dz[...,0],-dr[...,1]-value[...,1]/r),axis=-1)
        expected=result['induction'][name]
        assert np.linalg.norm(negative_curl-expected)/np.linalg.norm(expected)<2e-7


def test_full_parent_mode_fields_and_rates():
    m=ArchivedContactMode();r=np.array([.009,.017,.029]);z=np.array([-.083,-.014,.071])
    a=terms(m,r,z);f=m.fields(r,z)
    assert np.linalg.norm(sum(a['electric'].values())-f['electric'])/np.linalg.norm(f['electric'])<2e-14
    for actual,target in [(sum(a['momentum'].values()),m.growth*f['physical'][...,1:4]),(sum(a['induction'].values()),m.growth*f['physical'][...,6:9])]:
        assert np.linalg.norm(actual-target)/np.linalg.norm(target)<2e-7
