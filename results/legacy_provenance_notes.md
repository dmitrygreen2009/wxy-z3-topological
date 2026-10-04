# Launch provenance interpretation

The original parent commit did not contain the calculation source. Historical
results therefore have explicit provenance gaps. Preserve their audits and logs;
do not replace unknown values with plausible launch settings.

Version-3 processes initially read source fingerprints and Git state at each
solver-stage start. A process already running when files were edited can thus
record the disk version at that time rather than the code already loaded in
memory. The initial stage fingerprints, saved scripts, first source commit, and
recorded change history provide the available reconstruction evidence. Treat
stage-time disk hashes as such, rather than claiming they certify loaded code.

The updated audit implementation captures Git state and all source fingerprints
once at process launch. Later stage records retain that immutable snapshot and
identify it with `provenance_snapshot: process launch`. This applies to newly
started jobs; ongoing jobs continue without interruption.
