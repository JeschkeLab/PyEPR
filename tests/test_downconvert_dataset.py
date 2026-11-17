from pyepr.dataset import create_dataset_from_axes, downconvert_dataset
import numpy as np
import pytest

DC_IF = 0.15          # GHz
DC_FS = 1.0           # GHz, detector sampling rate
DC_N_TX = 256
DC_ECHO_POS = 128     # samples
DC_ECHO_WIDTH = 10    # samples
DC_PHASE = 0.3        # rad
DC_T_AMP = np.exp(-np.linspace(0, 1, 8))
DC_SCAN_AMP = np.array([1.0, 0.5, 2.0])


def build_echo_dataset(complex_data, with_scan):
    """Gaussian echo on an IF carrier, with dims ('t','tx') or ('scan','t','tx').

    Returns the dataset and the expected echo amplitude for every non-tx point.
    """
    tx = np.arange(DC_N_TX)
    t = np.linspace(0, 1, DC_T_AMP.shape[0])
    envelope = np.exp(-0.5 * ((tx - DC_ECHO_POS) / DC_ECHO_WIDTH)**2)
    carrier_phase = 2*np.pi*DC_IF*tx/DC_FS + DC_PHASE
    if complex_data:
        carrier = np.exp(1j*carrier_phase)
    else:
        carrier = np.cos(carrier_phase)

    amp = DC_T_AMP
    axes = [t, tx]
    labels = ['t', 'tx']
    if with_scan:
        amp = DC_SCAN_AMP[:, None] * amp[None, :]
        axes = [np.arange(DC_SCAN_AMP.shape[0])] + axes
        labels = ['scan'] + labels

    data = amp[..., None] * envelope * carrier
    params = {'IF_freq': DC_IF, 'det_rate': DC_FS}
    dset = create_dataset_from_axes(data, axes, params=params, axes_labels=labels)
    return dset, amp


complex_param = pytest.mark.parametrize(
    'complex_data', [True, False], ids=['complex', 'real'])
scan_param = pytest.mark.parametrize(
    'with_scan', [False, True], ids=['t_tx', 'scan_t_tx'])
filter_param = pytest.mark.parametrize('filter_type', ['boxcar', 'cheby'])


@complex_param
@scan_param
def test_downconvert_dataset_mixing_only(complex_data, with_scan):
    dset, amp = build_echo_dataset(complex_data, with_scan)
    dc = downconvert_dataset(dset, filter_type=None)

    assert dc.dims == dset.dims
    assert dc.shape == dset.shape
    assert np.iscomplexobj(dc)

    envelope = np.exp(-0.5 * ((dset.tx.values - DC_ECHO_POS) / DC_ECHO_WIDTH)**2)
    expected = amp[..., None] * envelope * np.exp(1j*DC_PHASE)
    np.testing.assert_allclose(dc.values, expected, atol=1e-6)


@complex_param
@scan_param
@filter_param
def test_downconvert_dataset_reduced(complex_data, with_scan, filter_type):
    dset, amp = build_echo_dataset(complex_data, with_scan)
    dc = downconvert_dataset(dset, filter_type=filter_type, reduce=True)

    assert dc.dims == tuple(d for d in dset.dims if d != 'tx')
    assert dc.shape == amp.shape
    assert np.iscomplexobj(dc)
    for dim in dc.dims:
        np.testing.assert_array_equal(dc[dim].values, dset[dim].values)

    assert dc.attrs['filter_type'] == filter_type
    assert abs(dc.attrs['max_echo_pos'] - DC_ECHO_POS) <= 1

    # Every trace is the same echo scaled, so relative amplitudes are kept
    # and the downconverted signal has the same constant phase everywhere.
    mag = np.abs(dc.values)
    np.testing.assert_allclose(mag / mag.max(), amp / amp.max(), rtol=1e-6)
    np.testing.assert_allclose(np.angle(dc.values), DC_PHASE, atol=1e-2)


@complex_param
@scan_param
@filter_param
def test_downconvert_dataset_unreduced(complex_data, with_scan, filter_type):
    dset, amp = build_echo_dataset(complex_data, with_scan)
    dc_full = downconvert_dataset(dset, filter_type=filter_type, reduce=False)

    assert dc_full.dims == dset.dims
    assert dc_full.shape == dset.shape
    assert np.iscomplexobj(dc_full)

    peak_pos = np.abs(dc_full.values).argmax(axis=-1)
    assert np.all(np.abs(peak_pos - DC_ECHO_POS) <= 1)

    dc = downconvert_dataset(dset, filter_type=filter_type, reduce=True)
    np.testing.assert_allclose(
        dc.values, dc_full.values[..., dc.attrs['max_echo_pos']])


@scan_param
@filter_param
def test_downconvert_dataset_real_matches_complex(with_scan, filter_type):
    dset_complex, _ = build_echo_dataset(True, with_scan)
    dset_real, _ = build_echo_dataset(False, with_scan)

    dc_complex = downconvert_dataset(dset_complex, filter_type=filter_type, reduce=False)
    dc_real = downconvert_dataset(dset_real, filter_type=filter_type, reduce=False)

    scale = np.abs(dc_complex.values).max()
    np.testing.assert_allclose(dc_real.values, dc_complex.values, atol=1e-6*scale)


@complex_param
@filter_param
def test_downconvert_dataset_scans_match_summed(complex_data, filter_type):
    dset, _ = build_echo_dataset(complex_data, with_scan=True)
    dset_summed = dset.sum('scan', keep_attrs=True)

    dc = downconvert_dataset(dset, filter_type=filter_type)
    dc_summed = downconvert_dataset(dset_summed, filter_type=filter_type)

    assert dc.attrs['max_echo_pos'] == dc_summed.attrs['max_echo_pos']
    np.testing.assert_allclose(dc.sum('scan').values, dc_summed.values)

    # Each scan is downconverted independently of the others
    for i in range(dset.sizes['scan']):
        dc_scan = downconvert_dataset(
            dset.isel(scan=i), filter_type=filter_type,
            max_echo_pos=dc.attrs['max_echo_pos'])
        np.testing.assert_allclose(dc.isel(scan=i).values, dc_scan.values)

