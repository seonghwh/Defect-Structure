from __future__ import annotations

import math

import numpy as np
import pytest

from batio3_defects.solver import B_SITE_DENSITY_CM3, solve_equilibrium, solve_quenched


TK = 1150 + 273
TQK = 500 + 273
AIR_PO2 = np.array([0.21])
MN_TOTAL = 7.7e19
Y_TOTAL = 7.7e19
BACKGROUND_ACCEPTOR = 100.0e-6 * B_SITE_DENSITY_CM3


SAMPLES = {
    "undoped_ab0999": {
        "ratio_ab": 0.999,
        "mn_total": 0.0,
        "y_total": 0.0,
        "acc_cm3": BACKGROUND_ACCEPTOR,
        "kr_model": "undoped_effective",
        "expected": {
            "n": 4.1882342888108934e6,
            "p": 2.1747654037556146e19,
            "VO2": 5.401172981224034e18,
            "VBa2": 1.5500000000000014e19,
            "VTi4": 0.0,
            "Mn0": 0.0,
            "Mn1": 0.0,
            "Mn2": 0.0,
            "YBa": 0.0,
            "Acc": BACKGROUND_ACCEPTOR,
        },
    },
    "mn_ab0994": {
        "ratio_ab": 0.994,
        "mn_total": MN_TOTAL,
        "y_total": 0.0,
        "acc_cm3": 0.0,
        "kr_model": "mn_effective",
        "expected": {
            "n": 1.4342377485588606e7,
            "p": 6.350709318089828e18,
            "VO2": 8.984174986102456e19,
            "VBa2": 9.300000000000008e19,
            "VTi4": 0.0,
            "Mn0": 7.696579096853248e19,
            "Mn1": 3.4209022810248908e16,
            "Mn2": 8.657268903257547e9,
            "YBa": 0.0,
            "Acc": 0.0,
        },
    },
    "mn_ab0999": {
        "ratio_ab": 0.999,
        "mn_total": MN_TOTAL,
        "y_total": 0.0,
        "acc_cm3": 0.0,
        "kr_model": "mn_effective",
        "expected": {
            "n": 3.5978321909798995e7,
            "p": 2.5316430980201375e18,
            "VO2": 1.4277056993220735e19,
            "VBa2": 1.5500000000000014e19,
            "VTi4": 0.0,
            "Mn0": 7.691424297001583e19,
            "Mn1": 8.575697554274026e16,
            "Mn2": 5.444142060080585e10,
            "YBa": 0.0,
            "Acc": 0.0,
        },
    },
    "mny_ab0999": {
        "ratio_ab": 0.999,
        "mn_total": MN_TOTAL,
        "y_total": Y_TOTAL,
        "acc_cm3": 0.0,
        "kr_model": "mn_effective",
        "expected": {
            "n": 4.7747052703154045e10,
            "p": 1.907641732518341e15,
            "VO2": 8.106387535340188e12,
            "VBa2": 1.5500000000000014e19,
            "VTi4": 0.0,
            "Mn0": 3.1036767240938635e19,
            "Mn1": 4.592454171136226e19,
            "Mn2": 3.869104769910242e16,
            "YBa": Y_TOTAL,
            "Acc": 0.0,
        },
    },
}


def solve_air_value(sample: dict[str, object]):
    eq = solve_equilibrium(
        AIR_PO2,
        TK,
        sample["ratio_ab"],
        Mn_total=sample["mn_total"],
        Y_total=sample["y_total"],
        acc_cm3=sample["acc_cm3"],
        kr_model=sample["kr_model"],
    )
    quench = solve_quenched(
        AIR_PO2,
        TQK,
        sample["ratio_ab"],
        frozen_eq=eq,
        Mn_total_quench=sample["mn_total"],
        Y_total_quench=sample["y_total"],
        acc_cm3=sample["acc_cm3"],
        vo_equilibrates=True,
        kr_model=sample["kr_model"],
    )
    return quench.iloc[0]


@pytest.mark.parametrize("sample_id", SAMPLES)
def test_current_air_values_are_authoritative(sample_id: str) -> None:
    sample = SAMPLES[sample_id]
    row = solve_air_value(sample)

    assert row.KR_model == sample["kr_model"]
    assert row.ratio_AB == pytest.approx(sample["ratio_ab"])
    assert row.TQK == TQK
    assert row.pO2 == pytest.approx(float(AIR_PO2[0]))

    for column, expected in sample["expected"].items():
        assert row[column] == pytest.approx(expected, rel=2e-10, abs=1e-12)


@pytest.mark.parametrize("sample_id", SAMPLES)
def test_mn_and_y_site_balances(sample_id: str) -> None:
    sample = SAMPLES[sample_id]
    row = solve_air_value(sample)

    mn_total = row.Mn0 + row.Mn1 + row.Mn2
    assert mn_total == pytest.approx(sample["mn_total"], rel=1e-12, abs=1e-6)
    assert row.YBa == pytest.approx(sample["y_total"])


def test_background_acceptor_conversion() -> None:
    assert math.isclose(BACKGROUND_ACCEPTOR, 1.55e18, rel_tol=1e-15)
