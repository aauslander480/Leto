import pandas as pd

from leto.donor_weights import filter_weights


def test_filtered_weights_sum_to_one():
    frame = pd.DataFrame({
        "SEG_ID": [1, 1, 1, 2, 2],
        "STAND_ID": ["MU_1", "MU_1", "MU_1", "MU_2", "MU_2"],
        "PLT_CN": ["A", "B", "C", "D", "E"],
        "WEIGHT": [0.90, 0.06, 0.04, 0.03, 0.02],
    })
    result = filter_weights(frame, 0.05)
    sums = result.groupby("SEG_ID")["WEIGHT"].sum()
    assert all(abs(value - 1.0) < 1e-12 for value in sums)
    assert set(result.loc[result["SEG_ID"] == 2, "PLT_CN"]) == {"D"}
