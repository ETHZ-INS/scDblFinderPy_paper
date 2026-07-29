#!/usr/bin/env python3
"""Poll `nvidia-smi` at a fixed interval, logging GPU utilization/memory to CSV.

Runs until killed (SIGTERM/SIGINT), which is how the Snakefile uses it: start
this in the background before a GPU run, then kill it once the run finishes.
"""
import csv
import subprocess
import sys
import time


def main():
    out_path = sys.argv[1]
    interval = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["timestamp", "gpu_util_pct", "mem_used_mib", "mem_total_mib"])
        f.flush()
        while True:
            out = subprocess.run(
                ["nvidia-smi", "-i", "0",
                 "--query-gpu=utilization.gpu,memory.used,memory.total",
                 "--format=csv,noheader,nounits"],
                capture_output=True, text=True, check=True,
            ).stdout.strip()
            util, used, total = (x.strip() for x in out.split(","))
            writer.writerow([time.time(), util, used, total])
            f.flush()
            time.sleep(interval)


if __name__ == "__main__":
    main()
