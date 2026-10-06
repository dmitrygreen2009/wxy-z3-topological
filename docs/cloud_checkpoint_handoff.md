# Minimal cloud checkpoint handoff

Prepared 2026-10-05. No VM has been provisioned. The active Mac calculation
and its queue have not been interrupted. AWS CLI/profile access is absent.

## Machine and budget

Start with one AWS EC2 r7g.4xlarge Linux ARM64 Spot instance in us-east-1:
16 vCPU, 128 GiB RAM, Ubuntu 24.04 ARM64, 100 GB gp3 persistent storage
(baseline IOPS/throughput). Set the data volume DeleteOnTermination=false.
Keep off-instance checkpoint backups. No GPU, managed cluster, or NAT gateway
is required. Check current per-AZ Spot price and interruption history before
launch; use $0.40/hour as a purchasing ceiling, not an asserted live quote.
Public estimates inspected today were approximately $0.305-$0.390/hour:
$7.32-$9.36 per 24 compute hours. Storage, public IPv4, backup and transfer
charges are additional; provisionally budget about $8-$11/day all-in, verify
account/region pricing. No completion-time estimate is yet justified.

Sources:
- https://aws.amazon.com/ec2/instance-types/r7g/
- https://www.doit.com/compute/compute/aws/us-east-1/r7g.4xlarge
- https://sparecores.com/server/aws/r7g.4xlarge
- https://aws.amazon.com/ebs/pricing/
- https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/spot-interruptions.html

## Handoff sequence

1. Provide authenticated cloud access or an existing VM SSH endpoint. Launch
   only one worker. Retain local checkpoint originals throughout validation.
2. At a completed-sector boundary in the current finite number scan, take a
   stable copy of its results, completed checkpoints and input seeds. Avoid
   copying a live checkpoint payload and metadata from different generations.
   The nine immutable resource_restart_snapshot pairs can be copied now.
   Transfer the exact working source tree (including required uncommitted
   changes), Project.toml, Manifest.toml, small results and checkpoint pairs.
   Exclude caches, compiled packages, scratch and .git internals; clone Git
   separately if needed. About 4.43 GiB of .jls files were present on inspection.
   Verify payload SHA256 against JSON before and after transfer.
3. Install Linux aarch64 Julia 1.10.10, not the macOS Julia executable or depot.
   Instantiate the existing manifest without Pkg.update. Keep Julia and BLAS
   threads at one initially, matching the original runs. Extra vCPUs are not
   evidence of single-job speedup; RAM eliminates the Mac's memory bottleneck.
   Keep the same pinned package source revisions, seed/reset policy, solver
   algorithms, tolerances, cutoffs, stage budgets, physical/QN indices and targets.
4. On Linux, load each transferred checkpoint and construct the exact original
   Hamiltonian using validate_only below; check saved norms/QNs and measurements.
   Julia Serialization portability is not assumed from ARM64 alone. Require
   successful target-machine validation before releasing any live Mac job.
   Platform-dependent floating-point roundoff is possible, even with unchanged
   numerical controls. Run inexpensive existing regression tests for the new
   platform; do not rerun completed production calculations.
5. Resume one job at a time in existing scientific priority order. Use the
   target_cap from results/julia_checkpoint_restart_plans.json. Keep each job's
   outputs in its own working copy if two independent runs could write the same
   result name. Do not transplant Mac PIDs, launch timestamps or /private/tmp
   executable paths into a Linux process watcher. No scheduler redesign needed.
6. Back up completed atomic checkpoint pairs off-instance. Spot termination
   notice is normally only two minutes, so existing periodic checkpoints are
   the recovery mechanism; do not depend on an emergency final solver iteration.
   Restart from the newest valid pair after interruption, preserving the
   original future-stage budget. Release idle compute, retain needed data.

Example (run in the copied repository, with Julia 1.10.10 on PATH):

```sh
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Read/setup-only checkpoint validation, followed by the existing resume entrypoint:

```sh
julia --project=. -e 'include("scripts/resume_checkpoint.jl"); println(resume_checkpoint(ARGS[1],parse(Int,ARGS[2]);validate_only=true,preserve_stage_budget=true))' results/checkpoints/infinite_armchair_w2_cell2_N36_chi32_Sz0_seed7103_v3_qn_star_active_latest_resource_restart_snapshot.jls 64
julia --project=. -e 'include("scripts/resume_checkpoint.jl"); resume_checkpoint(ARGS[1],parse(Int,ARGS[2]);preserve_stage_budget=true)' results/checkpoints/infinite_armchair_w2_cell2_N36_chi32_Sz0_seed7103_v3_qn_star_active_latest_resource_restart_snapshot.jls 64
```

## Finite workload exceptions: do not silently change settings

The active scripts/cylinder_number_search.jl skips completed, load-verified
number-sector points. Transfer it at a completed-sector boundary and re-run
that driver with all saved point files present to continue pending sectors.
The generic finite resume_checkpoint currently sets noise=0, whereas the
number scan uses [1e-7,1e-8,0]. It is therefore NOT an exact-settings replacement
for a mid-sweep transfer of that scan. Do not use it for that purpose. If a
mid-sector migration is necessary, implement/validate remaining-sweep control
with the original noise and initialization/partner logic before stopping Mac.

The protected paused finite batch PID39835 is NOT fully restartable yet.
Its readable chi256 seed does not reproduce the live chi512 partial sweep and
batch control state. Leave it paused; transferring that seed is not permission
to terminate its only live state.

## Independence and second VM

From their saved snapshots, the armchair width2 U1 chi32->64 continuation and
zigzag width1 unrestricted warm-start chi64->128 continuation are independent
optimizations. Other infinite width/filling/initialization/cell-size runs are
also independently restartable; same-geometry runs may overwrite filenames
unless working directories are isolated. Finite number sectors generally do
not exchange optimizer state, but the existing driver may use an earlier
particle-hole partner as initialization: preserve that dependency/order.

A second VM could overlap independent jobs, reducing elapsed queue time. There
is no timing basis for a quantified overall speedup, and clean-sector production
and analysis remain gates. Start with one VM; add a second only if measured
remaining time in two useful independent branches makes it cost-effective.
Retain existing scientific priority classification: archived unrestricted
warm starts do not automatically become clean winding-sector evidence.
