"""Headless rendering through the CPU backend (no GPU needed)."""

import numpy as np
import pytest
import vizviz
from conftest import FIXTURES

pytestmark = pytest.mark.skipif(not vizviz.HAS_RENDER, reason="built without the render feature")

DARK_BACKGROUND = [23, 23, 28, 255]  # StylePreset::DarkPresentation.background() as RGBA8


def drawn(img, background):
    return (img != np.asarray(background, dtype=np.uint8)).any(axis=2)


def test_render_returns_an_rgba_image(hemoglobin):
    img = hemoglobin.render(64, 48)
    assert img.shape == (48, 64, 4)
    assert img.dtype == np.uint8
    assert img.flags.c_contiguous
    assert img[0, 0].tolist() == DARK_BACKGROUND, "corner is background"
    assert img[24, 32].tolist() != DARK_BACKGROUND, "center hits the molecule"
    mask = drawn(img, DARK_BACKGROUND)
    assert 0.1 < mask.mean() < 0.9


def test_styles_and_background(hemoglobin):
    white = hemoglobin.render(32, 32, style="publication_white")
    assert white[0, 0].tolist() == [255, 255, 255, 255]
    assert np.array_equal(white, hemoglobin.render(32, 32, style="white"))
    custom = hemoglobin.render(32, 32, background=(1, 2, 3, 4))
    assert custom[0, 0].tolist() == [1, 2, 3, 4]
    for style in ["dark_presentation", "glossy", "flat_cel"]:
        img = hemoglobin.render(32, 32, style=style)
        assert img.shape == (32, 32, 4)
    with pytest.raises(ValueError):
        hemoglobin.render(32, 32, style="neon")


def test_coloring_changes_pixels(hemoglobin):
    by_element = hemoglobin.render(48, 48, coloring="element")
    by_chain = hemoglobin.render(48, 48, coloring="chain")
    by_b = hemoglobin.render(48, 48, coloring="b_factor")
    assert not np.array_equal(by_element, by_chain)
    assert not np.array_equal(by_element, by_b)
    with pytest.raises(ValueError):
        hemoglobin.render(32, 32, coloring="plaid")


def test_camera_controls(hemoglobin):
    base = hemoglobin.render(48, 48)
    turned = hemoglobin.render(48, 48, yaw=1.0, pitch=0.4)
    assert not np.array_equal(base, turned)
    far = hemoglobin.render(48, 48, zoom=3.0)
    assert drawn(far, DARK_BACKGROUND).sum() < drawn(base, DARK_BACKGROUND).sum()
    with pytest.raises(ValueError):
        hemoglobin.render(0, 48)
    with pytest.raises(ValueError):
        hemoglobin.render(48, 48, zoom=0.0)


def test_write_png(hemoglobin, tmp_path):
    img = hemoglobin.render(40, 30)
    out = tmp_path / "hemoglobin.png"
    vizviz.write_png(out, img)
    data = out.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    # A non-contiguous slice is accepted (copied internally).
    vizviz.write_png(tmp_path / "crop.png", img[5:25, 10:30])
    with pytest.raises(ValueError):
        vizviz.write_png(tmp_path / "bad.png", img[:, :, :3])
    with pytest.raises(OSError):
        vizviz.write_png(tmp_path / "no-such-dir" / "x.png", img)


def test_session_render_uses_the_scene_state():
    session = vizviz.Session()
    with pytest.raises(ValueError):
        session.render(32, 32)  # nothing loaded
    sid = session.load(FIXTURES / "4HHB.cif")
    default = session.render(48, 48)
    assert default.shape == (48, 48, 4)
    assert np.array_equal(default, session.render(48, 48, id=sid))
    session.set_coloring(sid, "chain")
    assert not np.array_equal(default, session.render(48, 48))
    with pytest.raises(KeyError):
        session.render(32, 32, id=sid + 50)
    session.set_representation(sid, "ball_and_stick")
    with pytest.warns(UserWarning, match="spacefill"):
        session.render(32, 32)
