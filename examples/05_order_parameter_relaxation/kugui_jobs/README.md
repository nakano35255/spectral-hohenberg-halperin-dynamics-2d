# Kugui F1cpu replicas for order-parameter relaxation

This directory contains PBS scripts for long `sine2d` relaxation runs on Kugui.

Default layout:

- queue: `F1cpu`
- allocation: `1` node, `128` MPI ranks
- replicas: `16` samples
- parallel layout: `8` MPI ranks per sample
- walltime: `12:00:00`
- grid: `256 256`
- length: `8192 8192`
- timestep: `dt = 4.0`
- run time: `2.4e6`, i.e. `600000` steps
- snapshot interval: `48000`, i.e. `12000` steps
- snapshots per sample: `50`
- transport: `eta = 0.12`, `M[0,0] = 0.12`

The run length is chosen from the previous i8cpu timing: `50000` steps took
about `27` minutes, so `600000` steps should be roughly `5.5` hours with the
same `8` ranks per sample. The `12` hour PBS walltime leaves margin.

Submit from the repository root on Kugui:

```bash
qsub examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_F1cpu_16samples_T2400000.pbs
```

Build before running:

```bash
qsub -v BUILD_BEFORE_RUN=1 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_F1cpu_16samples_T2400000.pbs
```

Useful overrides:

```bash
qsub -v RUN_TIME=4800000.0,SNAPSHOT_DTOUT=96000.0,RUN_NAME=sine2d_grid256_L8192_dt4_T4800000 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_F1cpu_16samples_T2400000.pbs
```

To collect more samples without overwriting existing sample IDs:

```bash
qsub -v SAMPLE_OFFSET=16,RUN_NAME=sine2d_grid256_L8192_dt4_T2400000_batch2 \
  examples/05_order_parameter_relaxation/kugui_jobs/job_kugui_F1cpu_16samples_T2400000.pbs
```

Default output root:

```text
/work/k0565/k056500/spectral-hohenberg-halperin-dynamics-2d/examples/05_order_parameter_relaxation/raw_data/sine2d_grid256_L8192_dt4_T2400000
```

Override `WORK_BASE` if your Kugui work directory is different.
