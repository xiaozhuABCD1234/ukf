import csv
from pathlib import Path

import matplotlib.pyplot as plt

csv_path = Path("/home/ubuntu/code/ukf/output/ukf_results.csv")
if not csv_path.exists():
    raise FileNotFoundError(f"CSV not found: {csv_path}")

steps = []
true_x = []
true_y = []
true_z = []
est_x = []
est_y = []
est_z = []
error_mm = []

with csv_path.open() as f:
    reader = csv.DictReader(f)
    for row in reader:
        steps.append(int(row["step"]))
        true_x.append(float(row["true_x"]))
        true_y.append(float(row["true_y"]))
        true_z.append(float(row["true_z"]))
        est_x.append(float(row["est_x"]))
        est_y.append(float(row["est_y"]))
        est_z.append(float(row["est_z"]))
        error_mm.append(float(row["error_mm"]))

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

axes[0].plot(steps, true_x, label="True X")
axes[0].plot(steps, est_x, label="Est X", linestyle="--")
axes[0].plot(steps, true_y, label="True Y")
axes[0].plot(steps, est_y, label="Est Y", linestyle="--")
axes[0].plot(steps, true_z, label="True Z")
axes[0].plot(steps, est_z, label="Est Z", linestyle="--")
axes[0].set_ylabel("Position (m)")
axes[0].legend(ncol=3)
axes[0].grid(True, alpha=0.3)

axes[1].plot(steps, error_mm, color="tab:red", label="Position Error (mm)")
axes[1].set_xlabel("Step")
axes[1].set_ylabel("Error (mm)")
axes[1].grid(True, alpha=0.3)
axes[1].legend()

fig.suptitle("UKF Position Tracking")
output_path = Path("/home/ubuntu/code/ukf/output/ukf_results.png")
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(output_path, dpi=150)
print(f"Saved plot: {output_path}")
